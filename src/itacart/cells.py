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
    EXTENSION_LON_PRECISION,
    MAX_RESOLUTION,
    PRIME_MERIDIAN_LON,
    QUINARY_GRID_SIZE,
    RES1_DIGITS,
    RES1_SEPARATOR,
    refinement_alphabet,
)
from .exceptions import (
    AntemeridianError,
    DomainError,
    NonExistentCellError,
    ResolutionError,
)
from .geodesy import geodetic_to_sinusoidal, sinusoidal_to_geodetic
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


def _check_boundary(lon: float, lat: float, x_index: int) -> None:
    """Reject positions resolving onto a cell F3 does not model.

    Two conditions, both from section 3.2 of the paper and both deferred
    to F4:

    - within :data:`~itacart.constants.EXTENSION_LON_PRECISION` degrees of
      the antemeridian, where cells are trapezoids or fall in an extension
      zone;
    - ``X < 1``, where the cell touches the prime meridian and is an
      isosceles triangle rather than a parallelogram.

    Raising is the point. Returning a parallelogram index for a triangular
    cell would give a cadastral caller a plausible area for a parcel that
    is a different shape (``D-3.5``).
    """
    if abs(abs(lon) - ANTEMERIDIAN_LON) < EXTENSION_LON_PRECISION:
        raise AntemeridianError(
            f"lon {lon} lies within {EXTENSION_LON_PRECISION} degrees of the "
            "antemeridian, where cells are trapezoidal or fall in an "
            "extension zone; delivered in F4"
        )
    if x_index < 1:
        raise NonExistentCellError(
            f"position lon={lon}, lat={lat} resolves onto a prime-meridian "
            f"cell (X={x_index}), which is an isosceles triangle rather than "
            "a parallelogram; delivered in F4"
        )


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
    if lon == PRIME_MERIDIAN_LON:
        raise NonExistentCellError(
            "the prime meridian itself is covered by triangular cells; "
            "delivered in F4"
        )
    quadrant = ("N" if lat >= 0.0 else "S") + ("E" if lon >= 0.0 else "W")
    x, y = geodetic_to_sinusoidal(lon, lat)
    return _quantize(abs(x), abs(y), quadrant, resolution, lon, lat)


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
    lon, lat = sinusoidal_to_geodetic(x, y)
    quadrant = ("N" if y >= 0.0 else "S") + ("E" if x >= 0.0 else "W")
    return _quantize(abs(x), abs(y), quadrant, resolution, lon, lat)


def _quantize(
    ax: float, ay: float, quadrant: str, resolution: int, lon: float, lat: float
) -> str:
    """Descend the hierarchy in the quadrant's mirrored positive octant.

    ``ax`` and ``ay`` are absolute plane coordinates; ``lon`` and ``lat``
    are carried only so the boundary screen and its messages can name the
    position the caller actually gave.
    """
    u = ax + ay
    v = ay
    u_col = int(math.floor((u + FLOOR_EPSILON_M) / _L1))
    row = int(math.floor((v + FLOOR_EPSILON_M) / _L1))
    column = u_col - row
    _check_boundary(lon, lat, column)

    residual_u = u - u_col * _L1
    residual_v = v - row * _L1

    codes: list[str] = []
    side = _L1
    for level in range(2, resolution + 1):
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


def is_quadrant_boundary_cell(cell: str) -> bool | list[bool]:
    """Whether a cell's anchor lies on a quadrant axis.

    True when the anchor sits on the prime meridian or on the equator, so
    that the mirror cell in the adjacent quadrant reports the *same*
    geodetic position. The condition is exact, not a tolerance: the
    anchor's plane coordinates are sums of exact multiples of the cell
    side, so a cell either lands on the axis or does not.

    This is the predicate that scopes acceptance criterion 1 (``D-3.91``).
    Requantizing an anchor is a statement about a tie-breaking
    convention, not about geometry, and on an axis no convention can
    satisfy both claimants -- see ``B-3.1``, which is open and belongs to
    F4.

    Distinct from :func:`itacart.boundary.is_boundary_cell`, which asks
    whether a cell touches a *grid discontinuity* -- the meridian
    triangles, the antemeridian, the extension zones. A cell can be a
    quadrant-boundary cell here and perfectly ordinary there: a southern
    ``Y = 0000`` cell is a plain parallelogram covering plain territory.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    values = []
    for atom in iter_cells(cell):
        x, y, _, _ = _anchor_on_plane(atom)
        values.append(x == 0.0 or y == 0.0)
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

    Use :func:`cell_to_centroid` when a centre point is wanted.

    **Restriction.** The anchor identifies its cell uniquely only when
    :func:`is_quadrant_boundary_cell` is false. A vertex always lies on
    the cell's border, and on a quadrant axis the mirror cell in the
    adjacent quadrant owns the very same vertex: reflection in the x axis
    turns a southern cell's lower-left corner into its northern
    extremity, so ``NE(X/0000)`` and ``SE(X/0000)`` report identical
    positions. Two cells claim one point and no tie-breaking rule
    satisfies both. The paper resolves the analogous case on the prime
    meridian with triangular cells and the non-existence of ``X = 0`` in
    the west, and says nothing about the equator; resolving it is F4's
    (``B-3.1``, ``D-3.91``).

    Args:
        cell: Compositional index string.

    Returns:
        ``(lon, lat)`` for a single cell, or a positionally aligned list.

    Raises:
        ResolutionError: If the index addresses a whole quadrant, which
            has no anchor.
    """
    values = [
        sinusoidal_to_geodetic(*_anchor_on_plane(atom)[:2]) for atom in iter_cells(cell)
    ]
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
    values = [_anchor_on_plane(atom)[:2] for atom in iter_cells(cell)]
    return _per_cell(cell, list(values))  # type: ignore[return-value]


def cell_to_centroid(cell: str) -> tuple[float, float] | list[tuple[float, float]]:
    """Geodetic centroid of a cell.

    Computed on the projection plane, where the cell is an exact
    parallelogram and its centroid is the mean of its four vertices, then
    inverted back to geodetic (``D-3.6``). The plane is equal-area, so the
    area weighting the centroid depends on is the one the plane preserves.

    Averaging the four *geodetic* vertices instead would be wrong: it
    weights by longitude, which the projection compresses by
    ``cos(phi)``, and it breaks outright for any cell spanning the
    antemeridian.

    For a parallelogram the mean simplifies to ``(x, y + side/2)`` in the
    mirrored octant -- the anchor's abscissa, half a side up. Triangular
    and trapezoidal cells report their own centroid from F4 onward.

    Args:
        cell: Compositional index string.

    Returns:
        ``(lon, lat)`` for a single cell, or a positionally aligned list.
    """
    values = []
    for atom in iter_cells(cell):
        ring = _vertices_on_plane(atom)
        mean_x = math.fsum(vertex[0] for vertex in ring) / len(ring)
        mean_y = math.fsum(vertex[1] for vertex in ring) / len(ring)
        values.append(sinusoidal_to_geodetic(mean_x, mean_y))
    return _per_cell(cell, list(values))  # type: ignore[return-value]


def cell_to_boundary(
    cell: str, close: bool = False
) -> list[tuple[float, float]] | list[list[tuple[float, float]]]:
    """Geodetic vertices bounding a cell, counter-clockwise.

    Four vertices for a parallelogram, starting at the anchor. Vertex
    count will depend on
    :func:`itacart.boundary.cell_shape` once F4 lands: three for a
    prime-meridian triangle, four for a clipped trapezoid.

    The ring is counter-clockwise in every quadrant. That is not free:
    mirroring the NE octant into NW or SE reverses orientation, so the
    sequence is reversed back (``B-0.3``).

    Args:
        cell: Compositional index string.
        close: Repeat the first vertex at the end, as ring conventions
            such as GeoJSON require.

    Returns:
        A vertex list for a single cell, or a positionally aligned list
        of vertex lists.
    """
    values = []
    for atom in iter_cells(cell):
        ring = [sinusoidal_to_geodetic(x, y) for x, y in _vertices_on_plane(atom)]
        if close:
            ring.append(ring[0])
        values.append(ring)
    return _per_cell(cell, list(values))  # type: ignore[return-value]


def cell_to_polygon(cell: str) -> "Polygon | list[Polygon]":
    """Cell boundary as a Shapely polygon in EPSG:4326.

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
