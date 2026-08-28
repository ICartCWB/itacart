"""Boundary behaviour: prime meridian, antemeridian, extension zones.

This module is what makes the domain global, complete and unique, and so
carries OGC DGGS Core requirements 8, 9 and 10 on its own.

Three departures from the uniform parallelogram grid are specified in
section 3.2 of the paper:

**Prime meridian.** Cells straddling 0 degrees are isosceles triangles
whose base is twice their height, mirrored about the meridian. This
applies at resolution 1; finer resolutions subdivide the adjacent eastern
cells rather than creating separate western ones. Consequently
resolution-1 cells with ``X = 0`` in the western quadrants do not exist.

**Antemeridian.** Uniform areas cannot be held everywhere at 180 degrees,
so the eastern quadrants are extended over the inhabited landmasses that
the meridian crosses: Fiji and the Chukotka/Wrangel group. Antarctica is
excluded. Cells in the oceans and in Antarctica along that boundary
therefore have unequal areas.

**Trapezoids.** Where a parallelogram vertex on the side opposite the
prime meridian would cross a boundary line, that side is clipped to the
line, extending the longer base into a trapezoid. The rule applies only
to the cells in that condition, not to their subdivisions.

The model
---------

Every cell has a *nominal* ring on the parallels plane -- a parallelogram
in the ordinary case, an isosceles triangle in the prime-meridian column
-- and an *effective* ring, which is the nominal one clipped to its
quadrant's domain. An eastern quadrant is bounded by the meridian its
extension reaches: 180 degrees outside a zone, 182 or 190.5 degrees
inside one. A western quadrant is bounded by the same line seen from the
other side, so it gives up exactly what the eastern one gains. Both are
bounded by the pole as well.

Clipping is the single mechanism. A cell exists when its clipped ring
still has area; it is a trapezoid when the clip changed the ring; it
honours the nominal area when the clip left it alone. Every predicate
here is a reading of that one construction.

Triangle refinement
-------------------

The children of a meridian triangle are triangles, at every resolution.
Figure 4(c) labels all 25 children of one triangle, and the arrangement
it draws is the ordinary grid folded onto itself about its diagonal:
child ``(row i, column j)`` sits in triangle row ``min(i, j)``, at
``|i - j|`` steps east of the meridian when ``j >= i`` and the same
number of steps west when ``j < i``. The diagonal children -- ``A1``,
``B2``, ``C3``, ``D4``, ``E5`` -- straddle the meridian and are
themselves isosceles triangles with base twice height, so the
construction is closed under refinement and no western cell is ever
created.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, cast

from .constants import (
    ANTEMERIDIAN_LON,
    EXTENSION_ZONES,
    MERIDIAN_QUADRANT,
    RES1_MAX_INDEX,
    RES1_SEPARATOR,
    CellShape,
    ExtensionZone,
    refinement_alphabet,
)
from .exceptions import ResolutionError
from .geodesy import geodetic_to_sinusoidal, sinusoidal_to_geodetic
from .index import is_atomic, is_valid_index, iter_cells, split_components
from .resolutions import cell_size, linear_refinement_ratio

if TYPE_CHECKING:
    from shapely.geometry.base import BaseGeometry

__all__ = [
    "cell_shape",
    "is_valid_cell",
    "is_boundary_cell",
    "is_triangular_cell",
    "is_trapezoidal_cell",
    "is_extension_cell",
    "extension_zone",
    "extension_zone_for_point",
    "extension_bounds",
    "crosses_antemeridian",
    "is_equal_area_cell",
    "absorbs_border",
    "last_lattice_column",
]

_L1: float = cell_size(1)
"""Resolution-1 cell side, in metres."""

_AREA_EPSILON_M2: float = 1e-6
"""Below this a clipped ring counts as empty, in square metres.

One square micrometre: far under the square centimetre of resolution 13,
so it separates a genuinely empty clip from the rounding noise of the
intersection arithmetic without ever swallowing a real cell.
"""

_CLIP_EPSILON_M: float = 1e-6
"""Slack allowed before a vertex counts as outside the border, in metres.

The same micrometre the cell layer floors with. A vertex that the
projection round trip leaves a nanometre past the border must not turn an
untouched parallelogram into a trapezoid.
"""


# --------------------------------------------------------------------------
# The domain border
# --------------------------------------------------------------------------


def _rows_of_zone(zone: str) -> tuple[int, int]:
    """First and last resolution-1 row an extension zone occupies.

    The declared latitude limits of Figure 5 are multiples of half a
    degree and fall inside rows rather than on their edges. A limit
    crossing a row would split that row's cells into a part governed by
    the extension and a part governed by the antemeridian, and the two
    have different borders: the cell would need a stepped edge and could
    no longer be clipped by a single line.

    The zone is therefore realized on whole rows, every row the declared
    band touches. Growing rather than shrinking is the cadastrally safe
    direction, since it can only add ocean and never drop declared land.
    """
    spec = EXTENSION_ZONES[zone]
    first = abs(geodetic_to_sinusoidal(0.0, spec.lat_min)[1])
    last = abs(geodetic_to_sinusoidal(0.0, spec.lat_max)[1])
    low, high = min(first, last), max(first, last)
    return int(math.floor(low / _L1)), int(math.ceil(high / _L1)) - 1


ZONE_ROWS: dict[str, tuple[int, int]] = {
    name: _rows_of_zone(name) for name in EXTENSION_ZONES
}
"""Resolution-1 row span each zone is realized on, inclusive."""


def _zone_of_row(quadrant: str, row: int) -> str | None:
    """Extension zone governing one resolution-1 row of a quadrant, if any.

    Both quadrants of a hemisphere are governed: the eastern one because
    it grows, the western one because it gives way. ``NE`` reaches 190.5
    degrees over Chukotka and ``NW`` stops at 169.5 west in the same rows.
    """
    for name, spec in EXTENSION_ZONES.items():
        if quadrant[0] != spec.quadrant[0]:
            continue
        first, last = ZONE_ROWS[name]
        if first <= row <= last:
            return name
    return None


def _lon_limit(quadrant: str, row: int) -> float:
    """Absolute longitude bounding a quadrant in one row, in degrees.

    180 everywhere except in an extension zone, where the eastern
    quadrant reaches past the antemeridian by exactly as much as the
    western quadrant gives up.
    """
    zone = _zone_of_row(quadrant, row)
    if zone is None:
        return ANTEMERIDIAN_LON
    reach = ANTEMERIDIAN_LON - abs(EXTENSION_ZONES[zone].lon_limit)
    if quadrant[1] == "E":
        return ANTEMERIDIAN_LON + reach
    return ANTEMERIDIAN_LON - reach


def _x_border(quadrant: str, row: int, y_abs: float) -> float:
    """Border abscissa of a quadrant at one ordinate, in metres.

    Positive, measured from the prime meridian. Equation (1) is linear in
    longitude, so it evaluates past 180 degrees without special-casing --
    which is why :func:`~itacart.geodesy.geodetic_to_sinusoidal` does not
    range-check longitude.
    """
    if y_abs >= MERIDIAN_QUADRANT:
        return 0.0
    lat = sinusoidal_to_geodetic(0.0, y_abs)[1]
    return abs(geodetic_to_sinusoidal(_lon_limit(quadrant, row), lat)[0])


# --------------------------------------------------------------------------
# Plane geometry
# --------------------------------------------------------------------------


def _signed_area(ring: list[tuple[float, float]]) -> float:
    """Signed shoelace area, positive counter-clockwise, in square metres.

    Shifted onto the first vertex before summing. Absolute plane
    coordinates run to twenty million metres while a resolution-13 cell
    is a centimetre across, so differencing them last rather than first
    is catastrophic cancellation.
    """
    if len(ring) < 3:
        return 0.0
    origin_x, origin_y = ring[0]
    count = len(ring)
    return (
        math.fsum(
            (ring[i][0] - origin_x) * (ring[(i + 1) % count][1] - origin_y)
            - (ring[(i + 1) % count][0] - origin_x) * (ring[i][1] - origin_y)
            for i in range(count)
        )
        / 2.0
    )


def ring_area(ring: list[tuple[float, float]]) -> float:
    """Unsigned area of a plane ring, in square metres."""
    return abs(_signed_area(ring))


def ring_centroid(ring: list[tuple[float, float]]) -> tuple[float, float]:
    """Area centroid of a plane ring, in metres.

    The shoelace-weighted centroid, not the mean of the vertices. The two
    agree for a triangle and for a parallelogram and part company for a
    trapezoid, where the mean is pulled toward whichever base carries
    more vertices rather than more area.

    Shifted onto the first vertex for the same reason :func:`ring_area`
    is: the moments would otherwise difference twenty-million-metre
    coordinates at the very end.
    """
    origin_x, origin_y = ring[0]
    count = len(ring)
    cross = [
        (ring[i][0] - origin_x) * (ring[(i + 1) % count][1] - origin_y)
        - (ring[(i + 1) % count][0] - origin_x) * (ring[i][1] - origin_y)
        for i in range(count)
    ]
    twice = math.fsum(cross)
    if abs(twice) < _AREA_EPSILON_M2:  # pragma: no cover - guarded by callers
        return (
            math.fsum(x for x, _ in ring) / count,
            math.fsum(y for _, y in ring) / count,
        )
    offset_x = (
        math.fsum(
            (ring[i][0] - origin_x + ring[(i + 1) % count][0] - origin_x) * cross[i]
            for i in range(count)
        )
        / 3.0
        / twice
    )
    offset_y = (
        math.fsum(
            (ring[i][1] - origin_y + ring[(i + 1) % count][1] - origin_y) * cross[i]
            for i in range(count)
        )
        / 3.0
        / twice
    )
    return origin_x + offset_x, origin_y + offset_y


def _dedupe(ring: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Drop vertices repeated within the clip slack of their predecessor."""
    out: list[tuple[float, float]] = []
    for vertex in ring:
        if out and math.hypot(vertex[0] - out[-1][0], vertex[1] - out[-1][1]) <= (
            _CLIP_EPSILON_M
        ):
            continue
        out.append(vertex)
    while (
        len(out) > 1
        and math.hypot(out[0][0] - out[-1][0], out[0][1] - out[-1][1])
        <= _CLIP_EPSILON_M
    ):
        out.pop()
    return out


def _clip_half_plane(
    ring: list[tuple[float, float]],
    start: tuple[float, float],
    end: tuple[float, float],
    keep: tuple[float, float],
) -> list[tuple[float, float]]:
    """Sutherland-Hodgman clip of a ring against one line.

    The line runs through ``start`` and ``end``; the side holding ``keep``
    survives. A convex cell clipped by a single line stays convex and
    simple, which is what makes the vertex count of a trapezoid
    predictable.
    """
    start_x, start_y = start
    end_x, end_y = end

    def side(point: tuple[float, float]) -> float:
        return (end_x - start_x) * (point[1] - start_y) - (end_y - start_y) * (
            point[0] - start_x
        )

    reference = side(keep)
    if reference == 0.0:  # pragma: no cover - the keep point never lies on the line
        return list(ring)
    sign = math.copysign(1.0, reference)
    tolerance = _CLIP_EPSILON_M * (math.hypot(end_x - start_x, end_y - start_y) or 1.0)

    out: list[tuple[float, float]] = []
    for position, current in enumerate(ring):
        previous = ring[position - 1]
        distance_current = sign * side(current)
        distance_previous = sign * side(previous)
        current_in = distance_current >= -tolerance
        previous_in = distance_previous >= -tolerance
        if current_in != previous_in:
            fraction = distance_previous / (distance_previous - distance_current)
            out.append(
                (
                    previous[0] + fraction * (current[0] - previous[0]),
                    previous[1] + fraction * (current[1] - previous[1]),
                )
            )
        if current_in:
            out.append(current)
    return _dedupe(out)


# --------------------------------------------------------------------------
# The triangle fold (Figura 4)
# --------------------------------------------------------------------------


def child_position(code: str, level: int) -> tuple[int, int]:
    """Row and column of a refinement code inside its parent, zero-based.

    Row counts from the parent's equator-side edge and column from its
    meridian-side edge, as Figures 3(c) and 3(d) draw them.
    """
    alphabet = refinement_alphabet(level)
    size = linear_refinement_ratio(level)
    position = alphabet.index(code)
    return position // size, position % size


def child_code(row: int, column: int, level: int) -> str:
    """Refinement code of a child at one row and column.

    Inverse of :func:`child_position`.
    """
    size = linear_refinement_ratio(level)
    return refinement_alphabet(level)[row * size + column]


def meridian_child(row: int, column: int, size: int) -> tuple[int, int]:
    """Sub-row and signed offset of a child inside a prime-meridian triangle.

    The fold of Figure 4(c): the ordinary grid is reflected onto itself
    about its diagonal, so ``(i, j)`` and ``(j, i)`` land in the same
    sub-row on opposite sides of the meridian. The sub-row is
    ``min(i, j)`` and the offset is ``j - i``, counted in cell sides east
    of the meridian and negative west of it.

    Sub-row ``k`` holds ``2 * (size - k) - 1`` children: one triangle
    straddling the meridian and ``size - k - 1`` parallelograms on each
    side. Summed over the sub-rows that is exactly ``size * size``, which
    is what makes the refinement a partition.
    """
    return min(row, column), column - row


def meridian_child_grid(sub_row: int, offset: int) -> tuple[int, int]:
    """Grid row and column of a child. Inverse of :func:`meridian_child`."""
    if offset >= 0:
        return sub_row, sub_row + offset
    return sub_row - offset, sub_row


def meridian_geometry(cell: str) -> tuple[str, float, float, float, float, float]:
    """Descend the prime-meridian column of one atomic index.

    Returns ``(shape, x, y, side, x_sign, y_sign)``. For a triangle ``x``
    is zero -- a meridian triangle is always centred on the line -- and
    ``y`` is its base ordinate. For a parallelogram they are the anchor,
    and ``x_sign`` says which side of the meridian it fell on, which need
    not be the side its quadrant code names.
    """
    quadrant, base_code, *codes = split_components(cell)
    row = int(base_code.partition(RES1_SEPARATOR)[2])
    y_sign = -1.0 if quadrant[0] == "S" else 1.0
    side = _L1
    base = y_sign * row * _L1

    for position, code in enumerate(codes):
        level = 2 + position
        size = linear_refinement_ratio(level)
        sub = side / size
        sub_row, offset = meridian_child(*child_position(code, level), size)
        base += y_sign * sub_row * sub
        side = sub
        if offset == 0:
            continue

        # Off the meridian: an ordinary parallelogram, and every deeper
        # code is the ordinary subdivision from here on.
        x_sign = math.copysign(1.0, offset)
        sheared_u = abs(offset) * sub + abs(base)
        sheared_v = abs(base)
        for deeper, deeper_code in enumerate(codes[position + 1 :]):  # noqa: E203
            deeper_level = level + 1 + deeper
            deeper_size = linear_refinement_ratio(deeper_level)
            side /= deeper_size
            child_row, child_column = child_position(deeper_code, deeper_level)
            sheared_u += child_column * side
            sheared_v += child_row * side
        return (
            "parallelogram",
            x_sign * (sheared_u - sheared_v),
            y_sign * sheared_v,
            side,
            x_sign,
            y_sign,
        )

    return "triangle", 0.0, base, side, 1.0, y_sign


# --------------------------------------------------------------------------
# Nominal and effective rings
# --------------------------------------------------------------------------


def _nominal_ring(
    cell: str,
) -> tuple[CellShape, list[tuple[float, float]], str, int]:
    """Nominal shape, plane ring, quadrant and resolution-1 column of a cell.

    The ring is counter-clockwise and signed, in the convention
    :func:`itacart.cells.cell_to_boundary` publishes. Clipping has not
    been applied.

    The triangle's first vertex is the midpoint of its base, because that
    is the point the paper says the index denotes: "the cell index
    intersects with the prime meridian and functions as the midpoint of
    the base". It is not the lower-left corner that the parallelogram
    reports.
    """
    from .cells import _vertices_on_plane

    components = split_components(cell)
    if len(components) < 2:
        raise ResolutionError(
            f"{cell!r} is a whole quadrant (resolution 0), which has no "
            "boundary geometry"
        )
    quadrant = components[0]
    column = int(components[1].partition(RES1_SEPARATOR)[0])
    if column != 0:
        return "parallelogram", _vertices_on_plane(cell), quadrant, column

    shape, x, y, side, x_sign, y_sign = meridian_geometry(cell)
    if shape == "parallelogram":
        ring = [
            (x, y),
            (x + x_sign * side, y),
            (x, y + y_sign * side),
            (x - x_sign * side, y + y_sign * side),
        ]
        if _signed_area(ring) < 0.0:
            ring.reverse()
            ring.insert(0, ring.pop())
        return "parallelogram", ring, quadrant, column

    ring = [(-side, y), (side, y), (0.0, y + y_sign * side)]
    if _signed_area(ring) < 0.0:
        ring[1], ring[2] = ring[2], ring[1]
    return "triangle", ring, quadrant, column


_RES1_MAX_COLUMN: int = int(RES1_MAX_INDEX.partition(RES1_SEPARATOR)[0])
"""Greatest resolution-1 column the index grammar admits (Table 1)."""


def _side_quadrant(quadrant: str, x: float) -> str:
    """Quadrant whose border governs one abscissa, in one hemisphere.

    A cell off the prime meridian lies wholly on its own quadrant's side
    and this answers that quadrant. A meridian triangle straddles the
    line, and its western half is bounded by the western quadrant of the
    same hemisphere, which is a different border wherever an extension
    zone makes the east reach further than the west.
    """
    return quadrant[0] + ("W" if x < 0.0 else "E")


def _trapezoid_bases(
    ring: list[tuple[float, float]], quadrant: str, row: int
) -> tuple[float, float]:
    """Length of a cell's two bases once its outer side lies on the border.

    Measured from the vertex nearest the prime meridian on each of the
    cell's two horizontal sides out to the border at that same ordinate.
    Negative means the near vertex is itself outside the border, so no
    base of that ordinate lies inside the domain at all.

    With ``b`` the border abscissa, ``x`` the anchor and ``s`` the side,
    these are ``b(y0) - x`` and ``b(y1) - (x - s)``. The lean of the cell
    carries the upper side one whole ``s`` towards the meridian while the
    border retreats by much less, so the upper base is the one that
    extends and the lower is the one that shrinks. Exactly one of them
    exceeds ``s``.
    """
    ordinates = [y for _, y in ring]
    bases: list[float] = []
    for ordinate in (min(ordinates, key=abs), max(ordinates, key=abs)):
        inner = min(abs(x) for x, y in ring if y == ordinate)
        bases.append(_x_border(quadrant, row, abs(ordinate)) - inner)
    return bases[0], bases[1]


def _absorb_border(
    ring: list[tuple[float, float]], quadrant: str, row: int
) -> list[tuple[float, float]] | None:
    """Carry a ring's outer vertices onto the border at their own ordinates.

    The rule the paper states for the antemeridian, applied wherever the
    domain ends: "if a vertex of the parallelogram on the side opposite
    the prime meridian exceeds the boundary, then that side shall be
    constrained within the boundary line, thereby extending the longer
    base of the trapezoid". The outer side is *moved onto* the border, not
    cut back from it, so the strip between the last addressable cell and
    the border is covered rather than dropped, and the result keeps four
    vertices instead of collapsing to three or splitting into five.

    Each ordinate is treated on its own, and each side of the prime
    meridian on its own, so a meridian triangle is carried east and west
    against two different borders. A vertex on the meridian never moves.

    Returns ``None`` when no vertex exceeds its border, which is the
    ordinary interior cell and the signal that the nominal ring stands.
    """
    outers: list[tuple[int, float, float]] = []
    exceeds = False
    for ordinate in {y for _, y in ring}:
        for sign in (1.0, -1.0):
            positions = [
                position
                for position, (x, y) in enumerate(ring)
                if y == ordinate and x != 0.0 and math.copysign(1.0, x) == sign
            ]
            if not positions:
                continue
            outer = max(positions, key=lambda position: abs(ring[position][0]))
            border = _x_border(_side_quadrant(quadrant, sign), row, abs(ordinate))
            outers.append((outer, sign * border, ordinate))
            exceeds = exceeds or abs(ring[outer][0]) > border + _CLIP_EPSILON_M

    if not exceeds:
        return None
    moved = list(ring)
    for position, abscissa, ordinate in outers:
        moved[position] = (abscissa, ordinate)
    return _dedupe(moved)


def last_lattice_column(quadrant: str, row: int, side: float) -> int:
    """Greatest column of a square lattice row that exists inside the border.

    The lattice of side ``side`` is the one every resolution shares: cells
    are nested, so the same two base tests decide existence at any depth.
    A column exists when its equator-side base is not negative and its
    polar-side base is positive, which bounds it by the border at the
    lower ordinate and by the border one side up, plus one for the lean.
    The tighter of the two wins.

    Above roughly 18.5 degrees of latitude the border retreats faster than
    one side per side, so the upper bound is the binding one and the last
    column stands short of where the lower ordinate alone would put it.
    The cell that does survive absorbs everything out to the border, so
    the row is still covered.

    Args:
        quadrant: Quadrant code.
        row: Lattice row index, counted from the equator.
        side: Lattice cell side in metres.

    Returns:
        The greatest existing column index, which is ``-1`` when the row
        holds no cell at all.
    """
    lower_y = row * side
    res1_row = int(math.floor((lower_y + _CLIP_EPSILON_M) / _L1))
    lower = _x_border(quadrant, res1_row, lower_y) / side
    upper = _x_border(quadrant, res1_row, lower_y + side) / side + 1.0
    return min(int(math.floor(lower)), int(math.ceil(upper)) - 1)


def _exists_inside_border(
    ring: list[tuple[float, float]], quadrant: str, row: int, column: int
) -> bool:
    """Whether a cell has any part of itself inside the domain border.

    Off the prime meridian the test is on the two bases of
    :func:`_trapezoid_bases`: the equator-side base must not be negative,
    which is the same as saying the anchor is inside the border, and the
    polar-side base must be positive. The second is not implied by the
    first. The border leans away from the meridian faster than one cell
    side above roughly 18.5 degrees of latitude, so a cell can have its
    anchor inside and its whole upper side outside.

    On the prime meridian the two halves lie against two different
    borders and the apex is on the line, so the test is instead that the
    base spans a positive width. That is what keeps the polar row alive
    after its apex has been cut away.
    """
    if column == 0:
        ordinate = min((y for _, y in ring), key=abs)
        width = sum(
            _x_border(quadrant[0] + side, row, abs(ordinate)) for side in ("E", "W")
        )
        return width > _CLIP_EPSILON_M
    lower, upper = _trapezoid_bases(ring, quadrant, row)
    return lower >= -_CLIP_EPSILON_M and upper > _CLIP_EPSILON_M


def plane_ring(cell: str) -> tuple[CellShape, list[tuple[float, float]]]:
    """Effective plane ring of one atomic cell, and the shape it takes.

    The entry point every geometric operation in the package goes
    through. Builds the nominal ring, caps it at the pole, and carries any
    vertex that exceeds the domain border out onto the border itself. A
    cell that absorbs the border is a ``"trapezoid"``: it no longer
    honours the nominal area of its resolution.

    Existence is decided before the geometry is published, and by the
    bases rather than by leftover area. A cell exists when the base on its
    equator side is not negative -- that is, when its anchor is inside the
    border -- and the base on its polar side is positive. The prime
    meridian column is judged instead on the width of its base, since its
    two halves lie against two different borders and its apex is on the
    line.

    Args:
        cell: Atomic compositional index string.

    Returns:
        ``(shape, ring)``, the ring counter-clockwise on the parallels
        plane in metres. The ring is empty when the cell has no area
        inside the domain, that is, when it does not exist.
    """
    shape, ring, absorbed = _ring_state(cell)
    return shape, ring


def _ring_state(cell: str) -> tuple[CellShape, list[tuple[float, float]], bool]:
    """:func:`plane_ring` plus whether the cell absorbed the border.

    Absorption is what breaks the equal-area guarantee, and it is not the
    same fact as the shape name: the polar cell is still a triangle and
    still absorbs. Both callers that care about area go through here
    rather than testing ``shape != "trapezoid"``.
    """
    shape, ring, quadrant, column = _nominal_ring(cell)
    ordinates = [abs(y) for _, y in ring]
    working = list(ring)
    row = int(math.floor(min(ordinates) / _L1 + _CLIP_EPSILON_M))

    if max(ordinates) > MERIDIAN_QUADRANT:
        y_sign = -1.0 if quadrant[0] == "S" else 1.0
        cap = y_sign * MERIDIAN_QUADRANT
        working = _clip_half_plane(
            working,
            (-_L1, cap),
            (_L1, cap),
            (0.0, y_sign * MERIDIAN_QUADRANT / 2.0),
        )

    if not _exists_inside_border(working, quadrant, row, column):
        return shape, [], False

    absorbed = _absorb_border(working, quadrant, row)
    if absorbed is None:
        return shape, working, False
    if column == 0:
        # The meridian triangle keeps its name even when the border has
        # come nearer than one side and collapsed its apex onto the pole.
        return shape, absorbed, True
    return "trapezoid", absorbed, True


def absorbs_border(cell: str) -> bool | list[bool]:
    """Whether a cell's outer side has been carried onto the domain border.

    The construction-level fact behind :func:`is_equal_area_cell` and
    :func:`itacart.resolutions.effective_cell_area`. True for every
    trapezoid and, additionally, for the polar triangle, whose apex the
    border has collapsed onto the pole.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    values: list[object] = []
    for atom in iter_cells(cell):
        try:
            values.append(_ring_state(atom)[2])
        except (ResolutionError, ValueError, KeyError, IndexError):
            values.append(False)
    return cast("bool | list[bool]", _per_cell(cell, values))


def _safe_ring(cell: str) -> tuple[CellShape, list[tuple[float, float]]]:
    """:func:`plane_ring` answering an empty ring instead of raising."""
    try:
        return plane_ring(cell)
    except (ResolutionError, ValueError, KeyError, IndexError):
        return "parallelogram", []


def to_geodetic(x: float, y: float) -> tuple[float, float]:
    """Invert a plane vertex, answering the pole itself at the quadrant.

    The meridian quadrant is the open end of
    :func:`~itacart.geodesy.inverse_meridian_arc`, and the polar row's
    clipped cells put vertices exactly on it.
    """
    if abs(y) >= MERIDIAN_QUADRANT:
        return 0.0, math.copysign(90.0, y)
    return sinusoidal_to_geodetic(x, y)


def _per_cell(index: str, values: list[object]) -> object:
    """Scalar for an atomic index, aligned list for a composed one."""
    return values[0] if is_atomic(index) else values


def _reaches_past_antemeridian(ring: list[tuple[float, float]]) -> bool:
    """Whether any vertex of a ring lies past the antemeridian."""
    return any(abs(to_geodetic(x, y)[0]) > ANTEMERIDIAN_LON for x, y in ring)


# --------------------------------------------------------------------------
# Shape classification
# --------------------------------------------------------------------------


def cell_shape(cell: str) -> CellShape | list[CellShape]:
    """Geometric class of a cell.

    Drives vertex count in :func:`itacart.cells.cell_to_boundary` and area
    computation in :func:`itacart.resolutions.effective_cell_area`.

    ``"triangle"`` is inherited and ``"trapezoid"`` is not, and the paper
    says both in the same breath. The meridian triangle is refined into
    triangles at every resolution, since the indexing "does not create
    separate western cells"; the trapezoid rule applies "solely to the
    cells within this specified condition and not to all subsequent
    resolutions", so a child of a clipped cell that falls wholly inside
    the border is an ordinary parallelogram again.

    A triangle that is itself clipped -- which happens only in the polar
    row, where the border comes nearer the meridian than one cell width
    -- reports ``"trapezoid"``, because the classification exists to
    answer whether the nominal area still holds.

    Args:
        cell: Compositional index string.

    Returns:
        One of ``"parallelogram"``, ``"triangle"`` or ``"trapezoid"`` for
        a single cell, or a positionally aligned list.
    """
    values = [_safe_ring(atom)[0] for atom in iter_cells(cell)]
    return cast("CellShape | list[CellShape]", _per_cell(cell, list(values)))


def is_triangular_cell(cell: str) -> bool | list[bool]:
    """Whether a cell is a prime-meridian isosceles triangle.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    values = [_safe_ring(atom)[0] == "triangle" for atom in iter_cells(cell)]
    return cast("bool | list[bool]", _per_cell(cell, list(values)))


def is_trapezoidal_cell(cell: str) -> bool | list[bool]:
    """Whether a cell has been clipped at a boundary line.

    Trapezoidal cells break the equal-area guarantee; this is the check
    that should gate any area-sensitive cadastral computation.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    values = [_safe_ring(atom)[0] == "trapezoid" for atom in iter_cells(cell)]
    return cast("bool | list[bool]", _per_cell(cell, list(values)))


def is_equal_area_cell(cell: str) -> bool | list[bool]:
    """Whether a cell honours the nominal area of its resolution.

    Asked of the construction, not of the shape name. A parallelogram and
    a prime-meridian triangle both carry the nominal area of their
    resolution, and a cell that has absorbed the border does not. The
    distinction is not the same as ``shape != "trapezoid"``: the polar
    cell is still named a triangle and still absorbs, and it holds an
    eighth of what its resolution promises. Comparing measured area
    against nominal would answer the same question through a tolerance;
    this answers it exactly.

    The triangle qualifies: base twice height, halved, is ``l * l``, the
    very area of the parallelogram it replaces. That is the "analogous
    properties to a parallelogram regarding base and height" the paper
    invokes to justify the choice, and it is why equal area survives the
    prime meridian. The exception the paper declares is the antemeridian,
    and the trapezoid is the only shape that meets it.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    values = [not absorbs_border(atom) for atom in iter_cells(cell)]
    return cast("bool | list[bool]", _per_cell(cell, list(values)))


# --------------------------------------------------------------------------
# Domain validity (OGC Req 8-10)
# --------------------------------------------------------------------------


def is_valid_cell(cell: str) -> bool | list[bool]:
    """Whether a cell exists in the ITACaRT domain.

    Stricter than :func:`itacart.index.is_valid_index`, which only checks
    syntax. This additionally applies:

    - the western-quadrant ``X = 0`` non-existence rule,
    - the prefix-only status of column ``2004``,
    - the resolution-1 extent bounds per quadrant,
    - antemeridian coverage, so oceanic cells beyond 180 degrees outside
      an extension zone are rejected.

    Column 2004 stands one past Table 1 and the grammar admits it only so
    that the fifth child of a trapezoidal cell of column 2003 can be
    spelled. It is a prefix, never a resolution-1 cell of its own, so this
    predicate refuses it at resolution 1 and accepts it deeper.

    The extent bounds and antemeridian coverage are one test: a cell
    exists when its ring still has area
    after the clip. The index space of resolution 1 is not a rectangle,
    and no bound of the form ``X <= 2003`` describes it. The greatest
    addressable column shrinks with the cosine of the latitude, from 2003
    on the equator to 313 in row 900 and none at all in the polar row.
    ``RES1_MAX_INDEX`` transcribes Table 1 faithfully but names the
    corner of a bounding box, not the last cell of any row.

    A syntactically malformed index is answered ``False`` rather than
    raised on, so the predicate composes.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    if not is_valid_index(cell):
        return False
    values: list[object] = []
    for atom in iter_cells(cell):
        components = split_components(atom)
        if len(components) < 2:
            values.append(True)
            continue
        quadrant = components[0]
        column = int(components[1].partition(RES1_SEPARATOR)[0])
        if column == 0 and quadrant[1] == "W":
            values.append(False)
            continue
        if column > _RES1_MAX_COLUMN and len(components) < 3:
            values.append(False)
            continue
        values.append(ring_area(_safe_ring(atom)[1]) > _AREA_EPSILON_M2)
    return cast("bool | list[bool]", _per_cell(cell, values))


def is_boundary_cell(cell: str) -> bool | list[bool]:
    """Whether a cell touches any grid discontinuity.

    True for prime-meridian triangles, antemeridian cells and
    extension-zone edge cells. A coarse screen; use :func:`cell_shape` or
    :func:`extension_zone` when the specific condition matters.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    values: list[object] = []
    for atom in iter_cells(cell):
        shape, ring = _safe_ring(atom)
        values.append(shape != "parallelogram" or _reaches_past_antemeridian(ring))
    return cast("bool | list[bool]", _per_cell(cell, values))


# --------------------------------------------------------------------------
# Extension zones
# --------------------------------------------------------------------------


def is_extension_cell(cell: str) -> bool | list[bool]:
    """Whether a cell lies inside an eastern-quadrant extension.

    True both for a cell wholly past the antemeridian and for the clipped
    cell on the zone's own meridian border, which is what makes it the
    companion of :func:`is_trapezoidal_cell` rather than a synonym for it.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    values: list[object] = []
    for atom in iter_cells(cell):
        values.append(_reaches_past_antemeridian(_safe_ring(atom)[1]))
    return cast("bool | list[bool]", _per_cell(cell, values))


def extension_zone(cell: str) -> ExtensionZone | None | list[ExtensionZone | None]:
    """Which extension zone contains a cell, if any.

    Callers that need to branch on the specific zone should use this
    rather than :func:`is_extension_cell`, since the two zones sit in
    different quadrants and have different longitude limits.

    Agrees with :func:`extension_zone_for_point` on the cell's centroid by
    construction: zones are realized on whole resolution-1 rows, so no
    cell straddles a latitude limit and none can disagree with its own
    interior.

    Args:
        cell: Compositional index string.

    Returns:
        ``"FIJI"``, ``"CHUKOTKA"`` or ``None`` for a single cell, or a
        positionally aligned list.
    """
    values: list[object] = []
    for atom in iter_cells(cell):
        components = split_components(atom)
        if len(components) < 2:
            values.append(None)
            continue
        row = int(components[1].partition(RES1_SEPARATOR)[2])
        values.append(_zone_of_row(components[0], row))
    return cast(
        "ExtensionZone | None | list[ExtensionZone | None]", _per_cell(cell, values)
    )


def extension_zone_for_point(lon: float, lat: float) -> ExtensionZone | None:
    """Which extension zone contains a geodetic position, if any.

    The point-level counterpart of :func:`extension_zone`, called during
    quantization before any cell exists. What quantization needs from it
    is which meridian bounds the quadrant the position falls in, so the
    zone is read as the region it governs: the latitude band, on the
    hemisphere the zone belongs to, out to the longitude the extension
    reaches.

    Rejected: restricting the zone to the strip between the antemeridian
    and the extension limit. That is the literal reading of the shaded
    box in Figure 5, and it leaves most of Viti Levu outside the zone
    named after Fiji -- which cannot be the intent of a rule written to
    bring those islands into a single quadrant.

    Args:
        lon: Longitude in decimal degrees.
        lat: Latitude in decimal degrees.

    Returns:
        The zone name, or ``None`` outside both zones.
    """
    ordinate = abs(geodetic_to_sinusoidal(0.0, lat)[1])
    row = int(math.floor(ordinate / _L1))
    zone = _zone_of_row(("S" if lat < 0.0 else "N") + "E", row)
    if zone is None:
        return None
    if lon >= 0.0 or lon <= EXTENSION_ZONES[zone].lon_limit:
        return cast(ExtensionZone, zone)
    return None


def extension_bounds(zone: ExtensionZone) -> tuple[float, float, float, float]:
    """Geodetic bounds of an extension zone.

    The declared limits of Figure 5, not the rows they are realized on.
    All six numbers are multiples of half a degree, the precision the
    paper adopted as the coarsest that clears every landmass except
    Antarctica.

    Args:
        zone: ``"FIJI"`` or ``"CHUKOTKA"``.

    Returns:
        ``(min_lon, min_lat, max_lon, max_lat)`` in decimal degrees, with
        ``min_lon`` on the western side of the antemeridian.

    Raises:
        ValueError: If ``zone`` is not a defined zone.
    """
    if zone not in EXTENSION_ZONES:
        raise ValueError(
            f"{zone!r} is not a defined extension zone; expected one of "
            f"{sorted(EXTENSION_ZONES)}"
        )
    spec = EXTENSION_ZONES[zone]
    return (-ANTEMERIDIAN_LON, spec.lat_min, spec.lon_limit, spec.lat_max)


def crosses_antemeridian(geometry: "BaseGeometry") -> bool:
    """Whether a geometry spans the 180th meridian.

    Called by :func:`itacart.geometry.polyfill` before descent, so the
    caller gets :class:`~itacart.exceptions.AntemeridianError` rather than
    a silently truncated cell set.

    Asked of the geometry as given, not of any cell it may cover. A ring
    written with longitudes past 180 degrees -- the natural way to
    describe an extension-zone footprint -- does not cross: it lies on
    one side already. What crosses is a ring whose consecutive vertices
    jump more than half the globe, which is what wrapping into
    ``[-180, 180]`` does to a shape spanning the line.

    Args:
        geometry: A Shapely geometry in EPSG:4326.

    Returns:
        ``True`` if the geometry crosses 180 degrees longitude.
    """
    from shapely.geometry import LineString, Polygon

    if geometry.is_empty:
        return False
    for part in getattr(geometry, "geoms", [geometry]):
        if isinstance(part, Polygon):
            sequences = [part.exterior.coords]
            sequences.extend(ring.coords for ring in part.interiors)
        elif isinstance(part, LineString):
            sequences = [part.coords]
        else:
            continue
        for coords in sequences:
            longitudes = [point[0] for point in coords]
            if any(
                abs(longitudes[position + 1] - longitudes[position]) > ANTEMERIDIAN_LON
                for position in range(len(longitudes) - 1)
            ):
                return True
    return False
