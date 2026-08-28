"""Quantization and inverse geometry: the core address <-> position mapping.

Satisfies OGC DGGS Core requirements 11 (simple cell geometry), 12 (direct
position) and 16 (quantization).

Origem: itacart_core/cells.py + sinusoidal_coordinates_to_dggs do notebook.

The shear
---------

ITACaRT cells are 45-degree parallelograms on the parallels plane, not
axis-aligned squares. For a lower-left vertex ``(x, y)`` and side ``l``
the vertices are ``(x, y)``, ``(x + l, y)``, ``(x, y + l)`` and
``(x - l, y + l)`` (paper, section 3, Figure 1).

A unit shear maps that tiling onto a square lattice::

    u = |x| + |y|      v = |y|        forward, Jacobian 1
    |x| = u - v        |y| = v        inverse

All indexing arithmetic -- the 2x2 and 5x5 subdivisions -- happens in
``(u, v)``, where it is ordinary integer flooring. This is the single
idea the whole module rests on, and it is what keeps the descent free of
trigonometry below resolution 1.

The base index reported at resolution 1 is ``X = u_col - Y`` and
``Y = j``, so ``X`` is the distance from the cell's representative point
to the prime meridian along the parallel, in 10 km units, and the anchor
recovers as ``x = X * l``. That is what makes the paper's neighbour rules
hold (section 3.1): the northern neighbour keeps ``X``, the eastern one
keeps ``Y``.

Boundary scope
--------------

Prime-meridian triangles, antemeridian trapezoids and the extension zones
are F4's. Until then a position resolving onto one of them raises rather
than returning a parallelogram that does not exist
(:class:`~itacart.exceptions.NonExistentCellError`,
:class:`~itacart.exceptions.AntemeridianError`).
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

from .constants import (
    ANTEMERIDIAN_LON,
    DESCENT_CLOSE,
    DESCENT_OPEN,
    MAX_RESOLUTION,
    QUINARY_GRID_SIZE,
    RES1_DIGITS,
    RES1_SEPARATOR,
    refinement_alphabet,
)
from .exceptions import DomainError, NonExistentCellError, ResolutionError
from .geodesy import geodetic_to_sinusoidal
from .index import is_atomic, is_valid_index, iter_cells, split_components
from .resolutions import cell_size, linear_refinement_ratio

if TYPE_CHECKING:
    from shapely.geometry import Polygon

__all__ = [
    "geo_to_cell",
    "cell_to_anchor",
    "cell_to_centroid",
    "cell_to_boundary",
    "cell_to_polygon",
    "sinusoidal_to_cell",
    "cell_to_sinusoidal",
    "is_quadrant_boundary_cell",
]

_L1: float = cell_size(1)
"""Resolution-1 cell side, in metres."""

FLOOR_EPSILON_M: float = 1e-6
"""Floor tolerance, in metres, applied wherever a coordinate meets a grid line.

One micrometre: four orders of magnitude below the finest cell (1 cm) and
well above the ~1e-9 m residual that a projection round trip leaves at
continental distances. Without it a position recovered from its own cell
anchor floors back into the previous cell, and ``geo_to_cell`` stops being
idempotent on its own output -- which is acceptance criterion 1.

Ported verbatim from ``itacart_core/cells.py`` (``_EPS_M``). It is applied
at three places and only three: the sheared column, the sheared row, and
each child selection during descent.
"""

_ROW_LETTERS: tuple[str, ...] = tuple(
    code[0] for code in refinement_alphabet(3)[::QUINARY_GRID_SIZE]
)
"""Row letters A..E of the quinary alphabet, south to north.

Sliced out of :func:`~itacart.constants.refinement_alphabet` rather than
written out, so the north-south orientation has exactly one definition in
the package. ``B-0.1`` was a second, inverted copy of that rule.
"""


# --------------------------------------------------------------------------
# Index assembly
# --------------------------------------------------------------------------


def _from_path(components: list[str]) -> str:
    """Assemble an atomic index from its per-level components.

    Inverse of :func:`itacart.index.split_components`. The origin calls
    this ``compositional_index.from_path``; the package has no public
    counterpart, so it lives here privately (``P-3.3``). It duplicates
    ``index._render_path``, and promoting one of the two is the chat
    ponte's call, not this phase's.
    """
    quadrant, *rest = components
    if not rest:
        return quadrant
    body = DESCENT_OPEN.join(rest)
    return f"{quadrant}{DESCENT_OPEN}{body}{DESCENT_CLOSE * len(rest)}"


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------


def _check_resolution(resolution: int) -> None:
    """Reject a resolution outside 1..13.

    Resolution 0 is a whole quadrant and is not addressable by
    quantization: no position is quantized *to* a hemisphere.
    """
    if not isinstance(resolution, int) or isinstance(resolution, bool):
        raise ResolutionError(
            f"resolution must be an int, got {type(resolution).__name__}"
        )
    if not 1 <= resolution <= MAX_RESOLUTION:
        raise ResolutionError(f"resolution {resolution} outside 1..{MAX_RESOLUTION}")


def _check_position(lon: float, lat: float) -> None:
    """Reject non-finite or out-of-range geodetic input."""
    for value, name in ((lon, "lon"), (lat, "lat")):
        if not math.isfinite(value):
            raise DomainError(f"{name} must be finite, got {value!r}")
    if not -90.0 <= lat <= 90.0:
        raise DomainError(f"lat {lat} outside [-90, 90]")
    if not -180.0 <= lon <= 180.0:
        raise DomainError(f"lon {lon} outside [-180, 180]")


def _resolve_extension(lon: float, lat: float) -> tuple[float, str]:
    """Longitude and quadrant a position addresses, honouring the extensions.

    Two hemispheres meet at the antemeridian, and one of them may reach
    across it. A position west of 180 degrees inside Fiji's or
    Chukotka's latitude band belongs to the *eastern* quadrant, addressed
    by a longitude past 180: the projection is linear in longitude and
    carries straight on.

    The prime meridian and the antemeridian are both single lines with a
    quadrant on either side; both are awarded to the east, so that the
    cell owning a boundary point is the one whose index encodes it.

    Returns:
        ``(lon, quadrant)`` ready for projection and descent.
    """
    from .boundary import extension_zone_for_point
    from .constants import EXTENSION_ZONES

    if lon == -ANTEMERIDIAN_LON:
        lon = ANTEMERIDIAN_LON
    zone = extension_zone_for_point(lon, lat)
    if zone is not None and lon < 0.0 and lon <= EXTENSION_ZONES[zone].lon_limit:
        lon += 2.0 * ANTEMERIDIAN_LON
    return lon, ("N" if lat >= 0.0 else "S") + ("E" if lon >= 0.0 else "W")


def _descend_triangle(
    px: float, py: float, quadrant: str, row: int, resolution: int
) -> list[str]:
    """Refinement codes of the cell containing a position on the meridian.

    The resolution-1 cell of the prime-meridian column is an isosceles
    triangle spanning both sides of the line, so the ordinary descent
    does not apply to it. Its children do not all inherit the shape:
    exactly one per sub-row straddles the meridian and stays a triangle,
    pointing toward the pole with its base toward the equator, and the
    rest are ordinary parallelograms. Once the descent steps off the
    line it never returns, and the remaining levels are the ordinary
    subdivision.
    """
    from .boundary import child_code, meridian_child_grid

    y_sign = -1.0 if quadrant[0] == "S" else 1.0
    side = _L1
    base = y_sign * row * _L1
    codes: list[str] = []

    for level in range(2, resolution + 1):
        size = linear_refinement_ratio(level)
        sub = side / size
        height = (py - base) * y_sign
        sub_row = min(max(int((height + FLOOR_EPSILON_M) // sub), 0), size - 1)
        within = min(max((height + FLOOR_EPSILON_M) / sub - sub_row, 0.0), 1.0)

        # The triangle of this sub-row narrows from a full side at its
        # base to nothing at its apex; the parallelograms beside it lean
        # toward the meridian at the same rate.
        across = px / sub
        span = 1.0 - within + FLOOR_EPSILON_M / sub
        if abs(across) <= span:
            offset = 0
        else:
            steps = int(abs(across) + within + FLOOR_EPSILON_M / sub)
            steps = min(max(steps, 1), size - sub_row - 1)
            offset = steps if across > 0.0 else -steps

        codes.append(child_code(*meridian_child_grid(sub_row, offset), level))
        base += y_sign * sub_row * sub
        side = sub
        if offset != 0:
            anchor_x = abs(offset) * sub
            anchor_y = abs(base)
            codes.extend(
                _descend_parallelogram(
                    abs(px) - anchor_x, abs(py) - anchor_y, sub, level + 1, resolution
                )
            )
            break

    return codes


def _descend_parallelogram(
    residual_x: float, residual_y: float, side: float, first: int, resolution: int
) -> list[str]:
    """Refinement codes below a parallelogram, from residuals at its anchor.

    Shared by the ordinary descent and by the tail of a meridian descent
    that has stepped off the line, so the subdivision rule is written
    once.
    """
    residual_u = residual_x + residual_y
    residual_v = residual_y
    codes: list[str] = []
    for level in range(first, resolution + 1):
        divisor = linear_refinement_ratio(level)
        sub = side / divisor
        if divisor == 2:
            child_u = 1 if residual_u + FLOOR_EPSILON_M >= sub else 0
            child_v = 1 if residual_v + FLOOR_EPSILON_M >= sub else 0
            codes.append(str(child_v * 2 + child_u + 1))
        else:
            child_u = _bucket(residual_u, sub, divisor)
            child_v = _bucket(residual_v, sub, divisor)
            codes.append(f"{_ROW_LETTERS[child_v]}{child_u + 1}")
        residual_u -= child_u * sub
        residual_v -= child_v * sub
        side = sub
    return codes


# --------------------------------------------------------------------------
# Quantization (OGC Req 16)
# --------------------------------------------------------------------------


def geo_to_cell(lon: float, lat: float, resolution: int) -> str:
    """Address the cell containing a geodetic position.

    Projects onto the ellipsoidal parallels plane, mirrors into the
    quadrant's positive octant, shears onto the square lattice, then
    descends: the resolution-1 Cartesian pair, then alternating 1-to-4 and
    1-to-25 refinements. The descent stays in plane metres throughout; no
    trigonometry runs below resolution 1.

    Args:
        lon: Longitude in decimal degrees.
        lat: Latitude in decimal degrees.
        resolution: Target resolution level, 1 to 13.

    Returns:
        The atomic compositional index of the containing cell.

    Raises:
        ResolutionError: If ``resolution`` is out of range.
        DomainError: If the position is outside the addressable domain.
        NonExistentCellError: If the position falls on a prime-meridian
            triangular cell.
        AntemeridianError: If the position lies near the antemeridian.
    """
    _check_resolution(resolution)
    _check_position(lon, lat)
    lon, quadrant = _resolve_extension(lon, lat)
    x, y = geodetic_to_sinusoidal(lon, lat)
    return _quantize(x, y, quadrant, resolution)


def sinusoidal_to_cell(x: float, y: float, resolution: int) -> str:
    """Address a cell from projection-plane coordinates.

    Skips the forward projection; useful when coordinates are already on
    the plane, as in bulk cell filling. The boundary screen still runs, so
    it inverts the position back to geodetic first -- the antemeridian is
    a meridian, not a plane coordinate, and cannot be tested on ``x``
    alone once the parallel has been scaled by ``cos(phi)``.

    Args:
        x: Easting in metres on the sinusoidal plane.
        y: Northing in metres on the sinusoidal plane.
        resolution: Target resolution level, 1 to 13.

    Returns:
        The atomic compositional index of the containing cell.

    Raises:
        ResolutionError: If ``resolution`` is out of range.
        DomainError: If the position is outside the addressable domain.
    """
    _check_resolution(resolution)
    for value, name in ((x, "x"), (y, "y")):
        if not math.isfinite(value):
            raise DomainError(f"{name} must be finite, got {value!r}")
    quadrant = ("N" if y >= 0.0 else "S") + ("E" if x >= 0.0 else "W")
    return _quantize(x, y, quadrant, resolution)


def _quantize(x: float, y: float, quadrant: str, resolution: int) -> str:
    """Descend the hierarchy from signed plane coordinates.

    Mirrors into the quadrant's positive octant, shears onto the square
    lattice and floors. Column zero is the prime-meridian column, and it
    is not a parallelogram: the descent hands over to the triangular
    lattice and the quadrant is forced east, since the triangle spans
    both sides of the meridian and the western column does not exist.

    Between the last existing column of the target resolution and the
    domain border lies a strip the lattice cannot name. The cell that does
    exist there absorbs the border and covers it, so a position falling in
    that strip is carried back into that cell before the descent starts.
    Without this the quantizer would answer a column that
    :func:`itacart.boundary.is_valid_cell` refuses, and the strip would
    belong to no cell at all.
    """
    from .boundary import last_lattice_column

    absolute_x, absolute_y = abs(x), abs(y)
    u = absolute_x + absolute_y

    side = cell_size(resolution)
    fine_row = int(math.floor((absolute_y + FLOOR_EPSILON_M) / side))
    fine_column = int(math.floor((u + FLOOR_EPSILON_M) / side)) - fine_row
    last = last_lattice_column(quadrant, fine_row, side)
    if fine_column > last:
        u = (last + fine_row + 0.5) * side

    row = int(math.floor((absolute_y + FLOOR_EPSILON_M) / _L1))
    column = int(math.floor((u + FLOOR_EPSILON_M) / _L1)) - row

    if column <= 0:
        quadrant = quadrant[0] + "E"
        codes = _descend_triangle(x, y, quadrant, row, resolution)
        base = f"{0:0{RES1_DIGITS}d}{RES1_SEPARATOR}{row:0{RES1_DIGITS}d}"
        return _from_path([quadrant, base, *codes])

    codes = _descend_parallelogram(
        u - (column + row) * _L1 - (absolute_y - row * _L1),
        absolute_y - row * _L1,
        _L1,
        2,
        resolution,
    )

    base = f"{column:0{RES1_DIGITS}d}{RES1_SEPARATOR}{row:0{RES1_DIGITS}d}"
    return _from_path([quadrant, base, *codes])


def _bucket(residual: float, sub: float, divisor: int) -> int:
    """Floor a residual into one of ``divisor`` buckets, clamped.

    The clamp is a guard, not arithmetic: the residual is bounded by the
    parent side by construction, so the only way to land outside is
    floating-point slack at the top edge. The origin clamps here too.
    """
    return min(max(int(math.floor((residual + FLOOR_EPSILON_M) / sub)), 0), divisor - 1)


# --------------------------------------------------------------------------
# Inverse geometry (OGC Req 11, 12)
# --------------------------------------------------------------------------


def _anchor_on_plane(cell: str) -> tuple[float, float, float, str]:
    """Signed plane anchor, cell side and quadrant of one atomic cell.

    Returns ``(x, y, side, quadrant)``. The shared core of every inverse
    operation: ascends the component path accumulating ``(u, v)``, inverts
    the shear, then re-applies the quadrant signs.

    Resolution 0 has no anchor, and the reason is structural rather than
    a gap in the implementation. A quadrant's lower-left vertex is the
    origin of the projection, so all four quadrants would report the same
    position, ``(0, 0)``. An anchor that four distinct cells share is not
    a direct position, which is precisely the property OGC Core
    requirement 12 asks the anchor to have.

    Raising :class:`~itacart.exceptions.ResolutionError` rather than
    :class:`~itacart.exceptions.InvalidIndexError` (``D-3.7``): by
    ``D-2.3`` a bare quadrant code is a *valid* index, so claiming
    malformation would be false. What fails is the operation at this
    resolution, which is the case ``ResolutionError`` names -- "outside
    the valid range or invalid for the operation". Same distinction as
    ``D-2.8``, and it keeps this module consistent with
    :mod:`itacart.resolutions`, which already answers resolution 0 with
    ``ResolutionError`` in ``cell_size`` and ``nominal_cell_area``.

    Rejected alternatives: :class:`~itacart.exceptions.NonExistentCellError`,
    which would assert the quadrant cell does not exist -- it does, only
    its anchor is undefined; and
    :class:`~itacart.exceptions.MinResolutionError`, which is reserved
    for an attempt to ascend above resolution 0 and would be overloaded
    here, since nothing is ascending.
    """
    components = split_components(cell)
    if len(components) < 2:
        raise ResolutionError(
            f"{cell!r} is a whole quadrant (resolution 0), which has no "
            "anchor: all four quadrants meet at the origin of the "
            "projection, so the position would not identify the cell"
        )
    quadrant, base, *codes = components
    column_str, _, row_str = base.partition(RES1_SEPARATOR)
    column, row = int(column_str), int(row_str)

    u = (column + row) * _L1
    v = row * _L1
    side = _L1
    for offset, code in enumerate(codes):
        divisor = linear_refinement_ratio(2 + offset)
        sub = side / divisor
        if divisor == 2:
            digit = int(code) - 1
            child_u, child_v = digit % 2, digit // 2
        else:
            child_v = _ROW_LETTERS.index(code[0])
            child_u = int(code[1]) - 1
        u += child_u * sub
        v += child_v * sub
        side = sub

    ax = u - v
    ay = v
    x = -ax if quadrant[1] == "W" else ax
    y = -ay if quadrant[0] == "S" else ay
    return x, y, side, quadrant


def _reflections(quadrant: str) -> int:
    """How many axis reflections map the NE octant onto ``quadrant``.

    Zero for NE, one for NW and SE, two for SW. Each reflection reverses
    ring orientation, so an odd count means the mirrored vertex sequence
    must be reversed to stay counter-clockwise. Getting this wrong is
    ``B-0.3``: the origin re-signs the coordinates and stops there,
    yielding clockwise rings in exactly the two single-reflection
    quadrants.
    """
    return (1 if quadrant[1] == "W" else 0) + (1 if quadrant[0] == "S" else 0)


def _vertices_on_plane(cell: str) -> list[tuple[float, float]]:
    """The four plane vertices of a cell, counter-clockwise.

    Built in the mirrored octant where the parallelogram leans one way,
    then re-signed and, when the reflection count is odd, reversed.
    """
    x, y, side, quadrant = _anchor_on_plane(cell)
    sign_x = -1.0 if quadrant[1] == "W" else 1.0
    sign_y = -1.0 if quadrant[0] == "S" else 1.0
    ring = [
        (x, y),
        (x + sign_x * side, y),
        (x, y + sign_y * side),
        (x - sign_x * side, y + sign_y * side),
    ]
    if _reflections(quadrant) % 2 == 1:
        ring.reverse()
        ring.insert(0, ring.pop())  # keep the anchor first
    return ring


def _per_cell(index: str, values: list[object]) -> object:
    """Return a scalar for an atomic index and a list for a composed one.

    ``D-0.2``: the vectorised return is positionally aligned with
    :func:`itacart.index.decompose`, and a single-cell index still gets a
    bare value rather than a list of one.
    """
    return values[0] if is_atomic(index) else values


def _anchor_on_plane_any(cell: str) -> tuple[float, float]:
    """Signed plane anchor of one atomic cell, whatever its shape.

    A parallelogram reports its lower-left vertex, the vertex its index
    literally encodes. A triangle reports the midpoint of its base,
    because that is what the paper says its index denotes -- and the
    midpoint of a base is not a vertex of the cell at all.

    Cells in the prime-meridian column may be either: only the child that
    straddles the line stays a triangle, and its siblings are ordinary
    parallelograms whose anchor rule is the ordinary one.
    """
    from .boundary import meridian_geometry

    components = split_components(cell)
    if len(components) >= 2 and int(components[1].partition(RES1_SEPARATOR)[0]) == 0:
        _, x, y, _, _, _ = meridian_geometry(cell)
        return x, y
    return _anchor_on_plane(cell)[:2]


def is_quadrant_boundary_cell(cell: str) -> bool | list[bool]:
    """Whether a cell's anchor is shared with a cell in another quadrant.

    A cell's anchor is a vertex, and a vertex lies on the border of every
    cell that meets there. The half-open convention awards each border
    point to exactly one cell: a cell owns the ordinates from its own row
    up to, but not including, the next, and likewise across. That makes
    the anchor identify its own cell everywhere the convention and the
    quadrant signs agree.

    They disagree in one place. The equator is awarded to the north, so a
    southern cell in row ``0000`` has its anchor on the edge it does not
    own, and requantizing that anchor answers the northern cell above it.
    The prime meridian is the same situation resolved differently: it is
    awarded to the east, and the paper removes the western column
    outright rather than leave it holding an anchor it does not own.

    Removal is available at the meridian because the eastern triangle
    already covers both sides of it. It is not available at the equator,
    where the two rows cover disjoint ground and deleting one would erase
    ten kilometres of the southern hemisphere. The property is therefore
    unsatisfiable for one of the two rows, and this predicate names which.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    values = []
    for atom in iter_cells(cell):
        quadrant = split_components(atom)[0]
        _, y = _anchor_on_plane_any(atom)
        values.append(y == 0.0 and quadrant[0] == "S")
    return _per_cell(cell, list(values))  # type: ignore[return-value]


def cell_to_anchor(cell: str) -> tuple[float, float] | list[tuple[float, float]]:
    """Representative position of a cell, as defined by ITACaRT.

    ITACaRT designates the lower-left vertex, not the centroid, so that
    addressing behaves like a Cartesian system for surveyors. This is the
    position that satisfies OGC Core requirement 12, and it is also why
    EAERS requirement 27 is only partially met.

    "Lower-left" is meant in the cell's own quadrant: the vertex nearest
    the equator and nearest the prime meridian. It is the vertex the
    index literally encodes, so ``cell_to_anchor`` inverts the descent
    exactly, with no averaging.

    **Triangular cells differ.** A prime-meridian cell reports the
    midpoint of its base, which the paper names directly: the cell index
    "intersects with the prime meridian and functions as the midpoint of
    the base of an isosceles triangle". That point is not a vertex of the
    cell, and for a triangle refined into inverted children the base is
    the side away from the equator rather than toward it.

    Use :func:`cell_to_centroid` when a centre point is wanted.

    **Restriction.** The anchor identifies its cell uniquely except where
    :func:`is_quadrant_boundary_cell` is true, which is the southern row
    ``0000`` and nowhere else. See that function for why the exception
    cannot be removed.

    Args:
        cell: Compositional index string.

    Returns:
        ``(lon, lat)`` for a single cell, or a positionally aligned list.

    Raises:
        ResolutionError: If the index addresses a whole quadrant, which
            has no anchor.
    """
    from .boundary import to_geodetic

    values = [to_geodetic(*_anchor_on_plane_any(atom)) for atom in iter_cells(cell)]
    return _per_cell(cell, list(values))  # type: ignore[return-value]


def cell_to_sinusoidal(
    cell: str,
) -> tuple[float, float] | list[tuple[float, float]]:
    """Anchor of a cell on the projection plane.

    Args:
        cell: Compositional index string.

    Returns:
        ``(x, y)`` in metres for a single cell, or a positionally aligned
        list.
    """
    values = [_anchor_on_plane_any(atom) for atom in iter_cells(cell)]
    return _per_cell(cell, list(values))  # type: ignore[return-value]


def cell_to_centroid(cell: str) -> tuple[float, float] | list[tuple[float, float]]:
    """Geodetic centroid of a cell.

    Computed on the projection plane, where the cell's edges are straight
    and its centroid is the area-weighted one, then inverted back to
    geodetic. The plane is equal-area, so the weighting the centroid
    depends on is the one the plane preserves.

    Averaging the four *geodetic* vertices instead would be wrong: it
    weights by longitude, which the projection compresses by
    ``cos(phi)``, and it breaks outright for any cell spanning the
    antemeridian.

    All three shapes are handled. For a parallelogram and for a triangle
    the area centroid coincides with the mean of the vertices; for a
    clipped trapezoid it does not, and the mean would sit off centre.

    Args:
        cell: Compositional index string.

    Returns:
        ``(lon, lat)`` for a single cell, or a positionally aligned list.

    Raises:
        NonExistentCellError: If the cell has no area in the domain.
    """
    from .boundary import plane_ring, ring_centroid, to_geodetic

    values = []
    for atom in iter_cells(cell):
        ring = _require_ring(atom, plane_ring(atom)[1])
        values.append(to_geodetic(*ring_centroid(ring)))
    return _per_cell(cell, list(values))  # type: ignore[return-value]


def _require_ring(
    cell: str, ring: list[tuple[float, float]]
) -> list[tuple[float, float]]:
    """Reject a cell whose ring the domain clipped away to nothing."""
    if not ring:
        raise NonExistentCellError(
            f"{cell!r} has no area inside the ITACaRT domain: the column "
            "lies beyond the meridian bounding its quadrant in that row"
        )
    return ring


def cell_to_boundary(
    cell: str, close: bool = False
) -> list[tuple[float, float]] | list[list[tuple[float, float]]]:
    """Geodetic vertices bounding a cell, counter-clockwise.

    Four vertices for a parallelogram, three for a prime-meridian
    triangle, four for a trapezoid clipped on its meridian side. A
    trapezoid clipped where the border line leans faster than the cell's
    own side -- which happens above roughly 18.5 degrees of latitude,
    where the border crosses more than one column per row -- keeps three
    or five instead, since a line through a convex quadrilateral can
    leave any of those.

    The ring is counter-clockwise in every quadrant. That is not free:
    mirroring the NE octant into NW or SE reverses orientation, so the
    sequence is reversed back.

    Args:
        cell: Compositional index string.
        close: Repeat the first vertex at the end, as ring conventions
            such as GeoJSON require.

    Returns:
        A vertex list for a single cell, or a positionally aligned list
        of vertex lists.

    Raises:
        NonExistentCellError: If the cell has no area in the domain.
    """
    from .boundary import plane_ring, to_geodetic

    values = []
    for atom in iter_cells(cell):
        plane = _require_ring(atom, plane_ring(atom)[1])
        ring = [to_geodetic(x, y) for x, y in plane]
        if close:
            ring.append(ring[0])
        values.append(ring)
    return _per_cell(cell, list(values))  # type: ignore[return-value]


def cell_to_polygon(cell: str) -> "Polygon | list[Polygon]":
    """Cell boundary as a Shapely polygon in EPSG:4326.

    Valid and closed for all three shapes. A cell inside an extension
    zone carries longitudes past 180 degrees rather than wrapping, so the
    polygon stays simple and its area stays measurable.

    Args:
        cell: Compositional index string.

    Returns:
        A polygon for a single cell, or a positionally aligned list.
    """
    from shapely.geometry import Polygon

    rings = cell_to_boundary(cell, close=True)
    if is_atomic(cell):
        return Polygon(rings)
    return [Polygon(ring) for ring in rings]


def _is_valid_atomic(cell: str) -> bool:
    """Guard used by the tests to assert the assembler cannot emit garbage."""
    return is_valid_index(cell) and is_atomic(cell)
