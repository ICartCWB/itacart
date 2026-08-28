"""Tests for :mod:`itacart.boundary`, one section per acceptance criterion.

The module carries OGC DGGS Core requirements 8, 9 and 10 on its own, so
the tests are written as properties over enumerated families wherever the
family is finite. Boundary sets have effectively zero measure: sampling
finds them by luck and reports a number that means nothing.
"""

from __future__ import annotations

import math

import pytest
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import unary_union

from itacart import boundary, cells
from itacart.constants import (
    ANTEMERIDIAN_LON,
    EXTENSION_LON_PRECISION,
    EXTENSION_ZONES,
    MAX_RESOLUTION,
    MERIDIAN_QUADRANT,
    TRIANGLE_BASE_TO_HEIGHT_RATIO,
    refinement_alphabet,
)
from itacart.exceptions import NonExistentCellError
from itacart.geodesy import geodetic_to_sinusoidal
from itacart.resolutions import cell_size, effective_cell_area, nominal_cell_area

L1 = cell_size(1)


def plane_polygon(cell: str) -> Polygon:
    """Cell as a Shapely polygon on the projection plane.

    The plane is where the cell's edges are straight and its area is the
    true ellipsoidal one. Testing containment in longitude and latitude
    instead would measure the chord error of an undensified edge, which
    is F7's subject and not this module's.
    """
    return Polygon(boundary.plane_ring(cell)[1])


def descend(cell: str, code: str) -> str:
    """Append one refinement code to an atomic index."""
    depth = cell.count("(")
    return cell[:-depth] + "(" + code + ")" * (depth + 1)


def children(cell: str, level: int) -> list[str]:
    """Every child of an atomic cell at the next resolution."""
    return [descend(cell, code) for code in refinement_alphabet(level)]


def _last_column(quadrant: str, row: int) -> int:
    """Greatest existing resolution-1 column of a row.

    Reads :func:`boundary.last_lattice_column` rather than walking the
    row. The walk costs two thousand ``is_valid_cell`` calls per row and
    the enumerations below run it over every row of the globe;
    :func:`test_the_column_limit_agrees_with_a_walk_along_the_row` pins
    the two together so the shortcut cannot drift from the predicate.
    """
    return boundary.last_lattice_column(quadrant, row, L1)


def _last_column_by_walking(quadrant: str, row: int) -> int:
    """The same number, obtained only from the public predicate."""
    column = 0
    while boundary.is_valid_cell(f"{quadrant}({column + 1:04d}/{row:04d})"):
        column += 1
    return column


@pytest.mark.parametrize("row", [0, 1, 17, 100, 200, 465, 750, 900, 999, 1000])
@pytest.mark.parametrize("quadrant", ["NE", "NW", "SE", "SW"])
def test_the_column_limit_agrees_with_a_walk_along_the_row(
    quadrant: str, row: int
) -> None:
    """The closed form and the predicate answer the same row.

    ``last_lattice_column`` derives the bound from the two bases;
    ``is_valid_cell`` builds each ring and asks whether it has area.
    Agreement across the quadrants is what lets the enumerations use the
    cheap one.
    """
    assert _last_column(quadrant, row) == _last_column_by_walking(quadrant, row)


# --------------------------------------------------------------------------
# Criterion 1 -- the prime-meridian triangle
# --------------------------------------------------------------------------


@pytest.mark.parametrize("quadrant", ["NE", "SE"])
@pytest.mark.parametrize("row", [0, 1, 137, 500, 999])
def test_criterion_1_the_meridian_cell_is_an_isosceles_triangle(
    quadrant: str, row: int
) -> None:
    """Base twice the height, mirrored about the meridian, three vertices."""
    cell = f"{quadrant}(0000/{row:04d})"
    shape, ring = boundary.plane_ring(cell)
    assert shape == "triangle"
    assert len(ring) == 3

    base = [vertex for vertex in ring if abs(vertex[1]) == row * L1]
    apex = [vertex for vertex in ring if abs(vertex[1]) != row * L1]
    assert len(base) == 2 and len(apex) == 1

    width = abs(base[0][0] - base[1][0])
    height = abs(apex[0][1] - base[0][1])
    assert width == pytest.approx(TRIANGLE_BASE_TO_HEIGHT_RATIO * height)
    assert sorted(vertex[0] for vertex in base) == [-height, height]
    assert apex[0][0] == 0.0


@pytest.mark.parametrize("resolution", range(1, 14))
def test_criterion_1_the_triangle_carries_the_nominal_area(resolution: int) -> None:
    """Base ``2l`` and height ``l``, halved, is ``l * l``.

    This is the whole reason the paper chose a triangle: it has
    "analogous properties to a parallelogram regarding base and height",
    which is to say the same area. The prime meridian therefore costs the
    grid nothing in equal area, and the antemeridian remains the sole
    declared exception.
    """
    cell = "NE(0000/0500)"
    for level in range(2, resolution + 1):
        cell = descend(cell, refinement_alphabet(level)[0])
    area = boundary.ring_area(boundary.plane_ring(cell)[1])
    # The shoelace differences coordinates five million metres from the
    # origin. At resolution 13 the cell is a centimetre across, so the
    # last surviving bits are the float64 floor of that subtraction and
    # not the geometry: 5e6 * 2**-52 over a 1e-2 side is about 1e-7.
    assert area == pytest.approx(nominal_cell_area(resolution), rel=1e-6)
    assert boundary.is_equal_area_cell(cell) is True
    assert effective_cell_area(cell) == nominal_cell_area(resolution)


def test_criterion_1_the_triangle_tiles_its_row_with_no_gap_or_overlap() -> None:
    """The triangle is exactly the two half-columns it replaces.

    Its eastern edge is the western edge of ``NE(0001/Y)`` and its
    western edge is the eastern edge of ``NW(0001/Y)``. Three cells, one
    simple polygon, three times the nominal area.
    """
    parts = [
        plane_polygon(cell)
        for cell in ("NE(0000/0500)", "NE(0001/0500)", "NW(0001/0500)")
    ]
    for left in range(len(parts)):
        for right in range(left + 1, len(parts)):
            assert parts[left].intersection(parts[right]).area == 0.0
    row = unary_union(parts)
    assert row.geom_type == "Polygon"
    assert row.area == pytest.approx(3 * nominal_cell_area(1), rel=1e-12)


# --------------------------------------------------------------------------
# Criterion 2 -- finer resolutions create no western cells
# --------------------------------------------------------------------------


@pytest.mark.parametrize("level", [2, 3, 4, 5])
def test_criterion_2_triangle_children_partition_their_parent(level: int) -> None:
    """One triangle per sub-row, parallelograms beside it, and no gap.

    The triangle is a local anomaly of the meridian line, not a shape
    that propagates through a subtree. Sub-row ``k`` of a refinement by
    ``n`` holds one triangle straddling the line and ``n - k - 1``
    ordinary parallelograms on each side, which sums over the sub-rows to
    exactly ``n²`` -- so the count is what makes the partition possible
    in the first place.
    """
    parent = "NE(0000/0500)"
    for step in range(2, level):
        parent = descend(parent, refinement_alphabet(step)[0])
    kids = children(parent, level)
    polygons = [plane_polygon(cell) for cell in kids]
    size = int(math.isqrt(len(kids)))

    triangles = [c for c in kids if boundary.cell_shape(c) == "triangle"]
    assert len(triangles) == size
    assert all(len(boundary.plane_ring(c)[1]) == 3 for c in triangles)
    assert all(
        len(boundary.plane_ring(c)[1]) == 4
        for c in kids
        if boundary.cell_shape(c) == "parallelogram"
    )

    covered = unary_union(polygons)
    assert covered.symmetric_difference(plane_polygon(parent)).area == pytest.approx(
        0.0, abs=1e-9
    )
    assert sum(polygon.area for polygon in polygons) == pytest.approx(
        plane_polygon(parent).area, rel=1e-12
    )


@pytest.mark.parametrize("quadrant", ["NE", "SE"])
def test_criterion_2_every_triangle_points_at_the_pole(quadrant: str) -> None:
    """Base toward the equator, apex toward the pole, at every resolution.

    The refinement never produces an inverted triangle. Figure 4(b) of
    the paper draws one, and draws no interior parallelograms at all;
    Figure 4(c), which carries all 25 labels, is consistent with this and
    is the panel the package follows.
    """
    cell = f"{quadrant}(0000/0500)"
    sign = -1.0 if quadrant[0] == "S" else 1.0
    for level in range(2, 10):
        assert boundary.cell_shape(cell) == "triangle"
        ring = boundary.plane_ring(cell)[1]
        base = [v for v in ring if len([w for w in ring if w[1] == v[1]]) == 2]
        apex = [v for v in ring if v not in base][0]
        assert len(base) == 2
        assert (apex[1] - base[0][1]) * sign > 0.0, "apex must lie poleward"
        assert base[0][0] == -base[1][0], "base must be centred on the meridian"
        assert apex[0] == 0.0
        # the child that stays on the line is the diagonal one
        cell = descend(cell, refinement_alphabet(level)[0])


def test_criterion_2_the_diagonal_children_straddle_the_meridian() -> None:
    """The fold of Figure 4(c) closes: ``A1``..``E5`` are meridian cells.

    The children whose grid row and column agree are centred on the line,
    so refining a triangle always yields another triangle sitting on the
    meridian and never a western index. The off-diagonal children lie
    wholly on one side.
    """
    parent = "NE(0000/0500(1))"
    straddling = {
        code
        for code in refinement_alphabet(3)
        if plane_polygon(descend(parent, code)).bounds[0]
        < 0.0
        < plane_polygon(descend(parent, code)).bounds[2]
    }
    assert straddling == {"A1", "B2", "C3", "D4", "E5"}


def test_criterion_2_the_fold_is_a_bijection_in_both_alphabets() -> None:
    """Every grid cell reaches one position and every position one cell.

    The fold reads Figure 4(c): sub-row ``min(i, j)``, offset ``j - i``
    east of the meridian. Sub-row ``k`` must come out holding one cell at
    offset zero and ``size - k - 1`` on each side.
    """
    for size in (2, 5):
        seen = {
            boundary.meridian_child(row, column, size)
            for row in range(size)
            for column in range(size)
        }
        assert len(seen) == size * size
        for row in range(size):
            for column in range(size):
                folded = boundary.meridian_child(row, column, size)
                assert boundary.meridian_child_grid(*folded) == (row, column)
        per_row: dict[int, set[int]] = {}
        for sub_row, offset in seen:
            per_row.setdefault(sub_row, set()).add(offset)
        for sub_row, offsets in per_row.items():
            reach = size - sub_row - 1
            assert offsets == set(range(-reach, reach + 1))


# --------------------------------------------------------------------------
# Criterion 3 -- is_valid_cell
# --------------------------------------------------------------------------


@pytest.mark.parametrize("quadrant", ["NW", "SW"])
@pytest.mark.parametrize("row", [0, 1, 400, 999])
def test_criterion_3_the_western_meridian_column_does_not_exist(
    quadrant: str, row: int
) -> None:
    """ "Cells in resolution 1 with the X index equal to 0 in the western
    quadrants will be non-existent"."""
    assert boundary.is_valid_cell(f"{quadrant}(0000/{row:04d})") is False
    assert boundary.is_valid_cell(f"{quadrant}(0001/{row:04d})") is True


def test_criterion_3_the_column_limit_shrinks_with_the_cosine() -> None:
    """The resolution-1 index space is not the rectangle Table 1 quotes.

    ``RES1_MAX_INDEX`` transcribes the table correctly and names the
    corner of a bounding box. The region inside it is sheared: the last
    addressable column falls away with ``cos(phi)``, and past the last
    one the inverse shear still answers a number -- it is the longitude
    that leaves the planet.
    """
    limits = {}
    for row in (0, 300, 900, 1000):
        column = 0
        while boundary.is_valid_cell(f"NE({column + 1:04d}/{row:04d})"):
            column += 1
        limits[row] = column
    assert limits[0] == 2003
    assert limits[300] == 1784
    assert limits[900] == 311
    assert limits[1000] == 0
    assert limits[0] > limits[300] > limits[900] > limits[1000]

    # Two bounds, and above 18.5 degrees the upper one is the binding
    # one: the anchor of column 312 in row 900 is still inside the
    # border, but the cell's whole upper side is outside it, so no part
    # of that side lies in the domain and the cell cannot be closed.
    lower, upper = boundary._trapezoid_bases(
        boundary._nominal_ring("NE(0312/0900)")[1], "NE", 900
    )
    assert lower > 0.0
    assert upper < 0.0

    corner_lon, _ = cells.cell_to_anchor("NE(2003/1000)")
    assert abs(corner_lon) > ANTEMERIDIAN_LON
    assert boundary.is_valid_cell("NE(2003/1000)") is False


def test_criterion_3_a_malformed_index_is_answered_not_raised_on() -> None:
    """The predicate composes: syntax failure is ``False``, never an exception."""
    assert boundary.is_valid_cell("NE(0000/0000") is False
    assert boundary.is_valid_cell("XX(0001/0001)") is False
    assert boundary.is_valid_cell("NE(0001/0001(9))") is False


def test_criterion_3_a_bare_quadrant_is_valid() -> None:
    """Resolution 0 is a cell of the grid, and all four of them exist."""
    assert boundary.is_valid_cell("NE") is True
    assert boundary.is_valid_cell("SW") is True


# --------------------------------------------------------------------------
# Criterion 4 -- extension_zone_for_point
# --------------------------------------------------------------------------


def test_criterion_4_named_positions_classify_as_the_paper_names_them(
    suva: tuple[float, float],
    wrangel: tuple[float, float],
    central_pacific: tuple[float, float],
) -> None:
    assert boundary.extension_zone_for_point(*suva) == "FIJI"
    assert boundary.extension_zone_for_point(*wrangel) == "CHUKOTKA"
    assert boundary.extension_zone_for_point(*central_pacific) is None


def test_criterion_4_the_band_alone_is_not_enough() -> None:
    """Inside the latitude band but past the limit is not inside the zone.

    The Chukotka extension reaches 169.5 degrees west. A position further
    east than that, at the same latitude, is ordinary ``NW`` territory
    and the function must say so.
    """
    assert boundary.extension_zone_for_point(-172.0, 68.0) == "CHUKOTKA"
    assert boundary.extension_zone_for_point(-150.0, 68.0) is None
    assert boundary.extension_zone_for_point(-177.0, -18.0) is None
    assert boundary.extension_zone_for_point(177.0, -18.0) == "FIJI"


def test_criterion_4_the_zones_do_not_reach_the_wrong_hemisphere() -> None:
    """Fiji governs a southern band and Chukotka a northern one."""
    assert boundary.extension_zone_for_point(179.0, 18.0) is None
    assert boundary.extension_zone_for_point(179.0, -68.0) is None


# --------------------------------------------------------------------------
# Criterion 5 -- the limits of Figure 5
# --------------------------------------------------------------------------


def test_criterion_5_the_published_limits_are_reproduced() -> None:
    """Fiji to 178 W between 15.5 S and 21.5 S; Chukotka to 169.5 W between
    64 N and 72 N."""
    assert boundary.extension_bounds("FIJI") == (-180.0, -21.5, -178.0, -15.5)
    assert boundary.extension_bounds("CHUKOTKA") == (-180.0, 64.0, -169.5, 72.0)


@pytest.mark.parametrize("zone", ["FIJI", "CHUKOTKA"])
def test_criterion_5_every_limit_is_a_multiple_of_the_adopted_precision(
    zone: str,
) -> None:
    """ "We adopted a precision of 0.5 degrees" -- and it shows in all six.

    Not a coincidence and not a tolerance: half a degree is the coarsest
    step that clears every landmass except Antarctica, so the numbers are
    quantized to it by construction.
    """
    minimum_lon, minimum_lat, maximum_lon, maximum_lat = boundary.extension_bounds(zone)
    for value in (minimum_lat, maximum_lon, maximum_lat):
        assert math.isclose(
            value / EXTENSION_LON_PRECISION,
            round(value / EXTENSION_LON_PRECISION),
            abs_tol=1e-12,
        )
    assert minimum_lon == -ANTEMERIDIAN_LON


def test_criterion_5_an_undefined_zone_is_refused() -> None:
    with pytest.raises(ValueError, match="not a defined extension zone"):
        boundary.extension_bounds("ATLANTIS")  # type: ignore[arg-type]


# --------------------------------------------------------------------------
# Criterion 6 -- effective area against nominal
# --------------------------------------------------------------------------


def test_criterion_6_a_trapezoid_is_not_an_equal_area_cell() -> None:
    """Clipped cells lose the guarantee, and say so before they are measured."""
    cell = "NE(2003/0000)"
    assert boundary.cell_shape(cell) == "trapezoid"
    assert boundary.is_trapezoidal_cell(cell) is True
    assert boundary.is_equal_area_cell(cell) is False
    assert effective_cell_area(cell) != nominal_cell_area(1)


def test_criterion_6_the_triangle_keeps_the_guarantee() -> None:
    """The exception the paper declares is the antemeridian, not the meridian.

    "Preserving equal areas for all cells, with the exception of the
    antemeridian." The prime-meridian triangle is on the other boundary
    and keeps its area exactly.
    """
    for cell in ("NE(0000/0000)", "SE(0000/0731)", "NE(0000/0500(1(A1)))"):
        assert boundary.is_triangular_cell(cell) is True
        assert boundary.is_equal_area_cell(cell) is True
        assert effective_cell_area(cell) == nominal_cell_area(
            len(cells.split_path(cell)) if False else _resolution_of(cell)
        )


def _resolution_of(cell: str) -> int:
    from itacart.resolutions import get_resolution

    return get_resolution(cell)


def test_criterion_6_the_interior_is_untouched() -> None:
    """An ordinary cell far from any boundary carries the nominal area."""
    for cell in ("NE(1400/0374)", "SW(0476/0260)", "SE(0900/0200(3))"):
        assert boundary.cell_shape(cell) == "parallelogram"
        assert boundary.is_equal_area_cell(cell) is True
        assert effective_cell_area(cell) == nominal_cell_area(_resolution_of(cell))


# --------------------------------------------------------------------------
# Criterion 7 -- every position resolves to exactly one cell
# --------------------------------------------------------------------------

COVERAGE_POINTS: dict[str, tuple[float, float]] = {
    "north east quadrant": (100.0, 40.0),
    "north west quadrant": (-100.0, 40.0),
    "south east quadrant": (100.0, -40.0),
    "south west quadrant": (-100.0, -40.0),
    "prime meridian north": (0.0, 42.0),
    "prime meridian south": (0.0, -42.0),
    "just west of the meridian": (-0.02, 51.4779),
    "equator east": (100.0, 0.0),
    "equator west": (-100.0, 0.0),
    "equator on the meridian": (0.0, 0.0),
    "fiji zone, east of the line": (178.4419, -18.1416),
    "fiji zone, west of the line": (-178.5, -18.1416),
    "chukotka zone, east of the line": (179.0, 68.0),
    "chukotka zone, west of the line": (-179.4, 71.2333),
    "antemeridian outside any zone": (180.0, 42.0),
    "antemeridian outside any zone, south": (-180.0, -42.0),
    "polar row north": (12.0, 89.99),
    "polar row south": (-12.0, -89.99),
    "corner: origin": (0.0, 0.0),
    "corner: antemeridian equator": (180.0, 0.0),
    "corner: north pole": (0.0, 89.999),
    "corner: south pole": (0.0, -89.999),
}


@pytest.mark.parametrize("resolution", [1, 2, 3, 7, 13])
@pytest.mark.parametrize("name", sorted(COVERAGE_POINTS))
def test_criterion_7_every_named_position_resolves_to_exactly_one_cell(
    name: str, resolution: int
) -> None:
    """Enumerated, not sampled, and the enumeration names what it covers.

    The four quadrants, both axes, both extension zones on both sides of
    the line, the antemeridian outside any zone, the polar row and the
    corners of the index space. A position claimed by no cell and a
    position claimed by two are both failures, and the containment check
    is what distinguishes them.
    """
    lon, lat = COVERAGE_POINTS[name]
    cell = cells.geo_to_cell(lon, lat, resolution)
    assert boundary.is_valid_cell(cell) is True

    zone = boundary.extension_zone_for_point(lon, lat)
    plane_lon = ANTEMERIDIAN_LON if lon == -ANTEMERIDIAN_LON else lon
    if zone is not None and lon < 0.0 and lon <= EXTENSION_ZONES[zone].lon_limit:
        plane_lon = lon + 2.0 * ANTEMERIDIAN_LON
    # The border edge of a clipped cell is the chord of a meridian, not
    # the meridian itself, so a position lying exactly on the line can sit
    # up to the chord's sagitta outside the polygon. Measured maximum over
    # every resolution-1 row: 6.2 m, against a 10 km side.
    assert (
        plane_polygon(cell)
        .buffer(CHORD_SAGITTA_M)
        .contains(Point(*geodetic_to_sinusoidal(plane_lon, lat)))
    )


CHORD_SAGITTA_M = 10.0
"""Bound on the gap between a clipped cell's straight edge and its meridian."""


@pytest.mark.parametrize("resolution", [1, 2, 3])
def test_criterion_7_neighbouring_columns_do_not_overlap(resolution: int) -> None:
    """Uniqueness, stated as a property over a whole row rather than a point."""
    row = "0400"
    polygons = [
        plane_polygon(f"NE({column:04d}/{row})")
        for column in range(0, 6)
        if boundary.is_valid_cell(f"NE({column:04d}/{row})")
    ]
    for left in range(len(polygons)):
        for right in range(left + 1, len(polygons)):
            assert polygons[left].intersection(polygons[right]).area == pytest.approx(
                0.0, abs=1e-9
            )
    assert unary_union(polygons).geom_type == "Polygon"


# --------------------------------------------------------------------------
# Criterion 8 -- vertex counts
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cell",
    ["NE(0000/0500)", "SE(0000/0500)", "NW(0001/0500)", "SW(0001/0500)"],
)
def test_criterion_8_the_ring_is_counter_clockwise_in_every_quadrant(
    cell: str,
) -> None:
    """Mirroring reverses orientation; the sequence is reversed back."""
    ring = boundary.plane_ring(cell)[1]
    assert boundary._signed_area(ring) > 0.0


@pytest.mark.parametrize("quadrant", ["NE", "NW", "SE", "SW"])
def test_criterion_8_a_trapezoid_ring_is_counter_clockwise(quadrant: str) -> None:
    """Absorption must not flip the winding, in any quadrant.

    The sibling test above checks parallelograms and triangles, whose
    vertices are written in order by construction. A trapezoid's are not:
    two of them are displaced along their own horizontals to land on the
    border, and in the western and southern quadrants the ring has also
    been mirrored, which reverses orientation before the sequence is
    reversed back. A ring that came out clockwise would still have the
    right area under the shoelace but the wrong sign, and Shapely would
    hand a caller a polygon wound the wrong way.

    Enumerated over every row rather than spot-checked, because the
    displacement depends on how far the border has retreated and that
    varies with latitude.
    """
    for row in range(0, 1000):
        cell = f"{quadrant}({_last_column(quadrant, row):04d}/{row:04d})"
        if boundary.cell_shape(cell) != "trapezoid":
            continue
        ring = boundary.plane_ring(cell)[1]
        assert boundary._signed_area(ring) > 0.0, cell
        polygon = Polygon(ring)
        assert polygon.is_valid and polygon.is_simple, cell
        assert polygon.exterior.is_ccw, cell


@pytest.mark.parametrize("quadrant", ["NE", "NW", "SE", "SW"])
def test_criterion_8_the_polar_triangle_ring_is_counter_clockwise(
    quadrant: str,
) -> None:
    """The polar cell too, whose apex the border has collapsed to a point.

    Three vertices instead of four, produced by a different path -- the
    pole clip runs before absorption -- so its winding is worth asserting
    on its own rather than inferring from the trapezoid case.
    """
    cell = f"{quadrant[0]}E(0000/1000)"
    ring = boundary.plane_ring(cell)[1]
    assert len(ring) == 3
    assert boundary._signed_area(ring) > 0.0
    assert Polygon(ring).exterior.is_ccw


@pytest.mark.parametrize("quadrant", ["NE", "SE"])
def test_criterion_8_a_triangle_has_three_vertices(quadrant: str) -> None:
    cell = f"{quadrant}(0000/0500)"
    assert len(cells.cell_to_boundary(cell)) == 3
    assert len(cells.cell_to_boundary(cell, close=True)) == 4


def test_criterion_8_a_trapezoid_has_four_vertices() -> None:
    """Four, always. The border is absorbed, not cut.

    Figure 6's construction moves the outer side onto the border at both
    ordinates of the cell. A quadrilateral whose two outer vertices are
    displaced along their own horizontals is still a quadrilateral, so
    the count cannot change. Cutting the cell with the border line
    instead would leave three, four or five vertices depending on where
    the line entered, and that is the model this replaced.

    The two bases differ by exactly ``s - db``, the cell's own lean minus
    the border's retreat over the same height, so which of them ends up
    longer depends on latitude. Here, near the equator, it is the polar
    one.
    """
    cell = "NE(2003/0000)"
    assert boundary.cell_shape(cell) == "trapezoid"
    assert len(cells.cell_to_boundary(cell)) == 4

    plane = boundary.plane_ring(cell)[1]
    horizontals = [
        (left, right)
        for left, right in zip(plane, plane[1:] + plane[:1])
        if abs(left[1] - right[1]) < 1e-9
    ]
    assert len(horizontals) == 2
    lengths = sorted(abs(left[0] - right[0]) for left, right in horizontals)
    assert lengths[0] < L1 < lengths[1]


def _bases_of(cell: str) -> tuple[float, float]:
    """Equator-side and polar-side base of a cell, in nominal sides."""
    ring = boundary.plane_ring(cell)[1]
    ordinates = [y for _, y in ring]

    def width(ordinate: float) -> float:
        row = [x for x, y in ring if y == ordinate]
        return (max(row) - min(row)) / L1

    return width(min(ordinates, key=abs)), width(max(ordinates, key=abs))


def test_criterion_8_the_two_bases_differ_by_the_lean_minus_the_retreat() -> None:
    """``upper - lower == s - db``, to the last bit, in every row.

    The cell's outer side is carried onto the border at both ordinates, so
    the only things that separate the two bases are the cell's own lean --
    exactly one side per row -- and how far the border withdraws over the
    same height. Nothing else can enter, and the identity is what makes
    the base lengths predictable without measuring the projection twice.
    """
    worst = 0.0
    for row in range(0, 1000):
        cell = f"NE({_last_column('NE', row):04d}/{row:04d})"
        lower, upper = _bases_of(cell)
        ordinates = [abs(y) for _, y in boundary.plane_ring(cell)[1]]
        retreat = (
            boundary._x_border("NE", row, min(ordinates))
            - boundary._x_border("NE", row, max(ordinates))
        ) / L1
        worst = max(worst, abs((upper - lower) - (1.0 - retreat)))
    assert worst < 1e-12


def test_criterion_8_which_base_is_longer_changes_sides_at_the_slope_crossing() -> None:
    """The working statement of the phase was too strong, and this is what holds.

    "Exactly one base extends past the nominal side; the other is
    necessarily shorter" is what ``NE(2003/0000)`` does. It is not the
    rule: in 866 of the 4000 border cells of the globe **both** bases fall
    short. What holds everywhere is that they are never both longer, and
    that the longer of the two is the polar-side base below the latitude
    where the border starts leaning faster than the cell does, and the
    equator-side base above it.

    The crossing is row 205, at 18.5332 degrees -- the same ``k = 1``
    latitude that used to govern the vertex count, arriving here for the
    same reason but with a consequence that is now merely descriptive.
    """
    polar_longer = []
    both_short = 0
    for row in range(0, 1000):
        lower, upper = _bases_of(f"NE({_last_column('NE', row):04d}/{row:04d})")
        assert not (lower > 1.0 and upper > 1.0)
        both_short += int(lower < 1.0 and upper < 1.0)
        if upper > lower:
            polar_longer.append(row)
    assert polar_longer == list(range(0, 205))
    assert both_short == 216

    crossing = cells.cell_to_anchor("NE(0001/0205)")[1]
    assert crossing == pytest.approx(18.5332, abs=1e-4)


def test_criterion_8_every_trapezoid_on_the_globe_has_four_vertices() -> None:
    """Enumerated over every resolution-1 row, not sampled.

    The family of boundary cells has practically zero measure against the
    interior, so sampling finds them by luck and reports a number that
    means nothing. There are exactly 1000 of them in the north-east
    quadrant -- one per row, the polar row excepted, since there the
    meridian triangle is the only cell and it is not a trapezoid.
    """
    counts: dict[int, int] = {}
    per_row: list[int] = []
    for row in range(0, 1001):
        last = _last_column("NE", row)
        trapezoids = [
            f"NE({column:04d}/{row:04d})"
            for column in range(max(last - 6, 1), last + 1)
            if boundary.cell_shape(f"NE({column:04d}/{row:04d})") == "trapezoid"
        ]
        per_row.append(len(trapezoids))
        for cell in trapezoids:
            vertices = len(boundary.plane_ring(cell)[1])
            counts[vertices] = counts.get(vertices, 0) + 1
    assert counts == {4: 1000}
    assert set(per_row[:1000]) == {1}
    assert per_row[1000] == 0


# --------------------------------------------------------------------------
# Criterion 9 -- cell_to_polygon in all three shapes
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cell",
    [
        "NE(1400/0374)",
        "SW(0476/0260)",
        "NE(0000/0500)",
        "SE(0000/0500(3))",
        "NE(2003/0000)",
        "NE(0000/1000)",
    ],
)
def test_criterion_9_the_polygon_is_valid_and_closed_in_all_three_shapes(
    cell: str,
) -> None:
    polygon = cells.cell_to_polygon(cell)
    assert polygon.is_valid
    assert not polygon.is_empty
    assert polygon.exterior.is_ring


@pytest.mark.parametrize(
    "cell",
    [
        "NE(1400/0374)",
        "SW(0476/0260)",
        "NE(0000/0500)",
        "SE(0000/0500(3))",
        "NE(2003/0000)",
        "SE(0171/0100)",
    ],
)
def test_criterion_9_the_measured_area_matches_the_reported_one(cell: str) -> None:
    """Within 0.01 percent, measured on the plane the projection preserves."""
    measured = plane_polygon(cell).area
    reported = effective_cell_area(cell)
    assert measured == pytest.approx(reported, rel=1e-4)


def test_criterion_9_a_cell_with_no_ground_under_it_is_refused() -> None:
    """Geometry for a non-existent cell is an error, not an empty polygon."""
    assert boundary.is_valid_cell("NE(2003/1000)") is False
    with pytest.raises(NonExistentCellError, match="no area inside"):
        cells.cell_to_polygon("NE(2003/1000)")


# --------------------------------------------------------------------------
# Criterion 10 -- anchor and centroid in all three shapes
# --------------------------------------------------------------------------


def test_criterion_10_the_triangle_anchor_is_the_midpoint_of_its_base() -> None:
    """ "The cell index intersects with the prime meridian and functions as
    the midpoint of the base of an isosceles triangle."

    So it is not the lower-left vertex, and it is not a vertex at all:
    the anchor of a triangular cell lies on an edge, halfway along it.
    """
    cell = "NE(0000/0500)"
    anchor = cells.cell_to_anchor(cell)
    assert anchor[0] == 0.0
    plane = cells.cell_to_sinusoidal(cell)
    assert plane == (0.0, 500 * L1)
    assert plane not in boundary.plane_ring(cell)[1]

    ring = boundary.plane_ring(cell)[1]
    base = sorted(vertex for vertex in ring if vertex[1] == 500 * L1)
    assert plane[0] == pytest.approx((base[0][0] + base[1][0]) / 2.0)


def test_criterion_10_the_parallelogram_anchor_is_still_its_lower_left() -> None:
    """The rule the index encodes, unchanged by this phase."""
    assert cells.cell_to_sinusoidal("NE(0003/0007)") == (3 * L1, 7 * L1)


@pytest.mark.parametrize(
    "cell", ["NE(1400/0374)", "NE(0000/0500)", "SE(0000/0500(3))", "NE(2003/0000)"]
)
def test_criterion_10_the_centroid_is_inside_its_own_cell(cell: str) -> None:
    """True of all three shapes, which the mean of the vertices would not be."""
    centroid = cells.cell_to_centroid(cell)
    zone = boundary.extension_zone_for_point(*centroid)
    lon = centroid[0]
    if zone is not None and lon < 0.0 and lon <= EXTENSION_ZONES[zone].lon_limit:
        lon += 2.0 * ANTEMERIDIAN_LON
    assert (
        plane_polygon(cell)
        .buffer(1e-6)
        .contains(Point(*geodetic_to_sinusoidal(lon, centroid[1])))
    )


def test_criterion_10_the_trapezoid_centroid_is_area_weighted() -> None:
    """The shoelace centroid, not the mean of the vertices.

    They agree for a triangle and for a parallelogram and part company
    for a trapezoid, where the mean is pulled toward whichever base
    carries more vertices rather than more area.
    """
    ring = boundary.plane_ring("NE(2003/0000)")[1]
    weighted = boundary.ring_centroid(ring)
    naive = (
        math.fsum(x for x, _ in ring) / len(ring),
        math.fsum(y for _, y in ring) / len(ring),
    )
    assert weighted != naive
    assert Polygon(ring).buffer(1e-6).contains(Point(*weighted))


# --------------------------------------------------------------------------
# Criterion 12 -- construction in longitude
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("zone", "quadrant", "row"), [("CHUKOTKA", "NE", 750), ("FIJI", "SE", 200)]
)
def test_criterion_12_the_eastern_quadrant_accepts_columns_past_the_normal_limit(
    zone: str, quadrant: str, row: int
) -> None:
    """Inside a band the east reaches further, and stops at the zone's limit.

    The extra columns are not a tolerance: they are the arc of the
    parallel between 180 degrees and the zone's meridian, divided by the
    cell side. The count is measured against that arc.
    """
    outside_row = 500 if quadrant == "NE" else 400
    inside = _last_column(quadrant, row)
    outside = _last_column(quadrant, outside_row)
    assert inside > 0 and outside > 0

    latitude = cells.cell_to_anchor(f"{quadrant}(0001/{row:04d})")[1]
    reach = ANTEMERIDIAN_LON - abs(EXTENSION_ZONES[zone].lon_limit)
    plain = abs(geodetic_to_sinusoidal(ANTEMERIDIAN_LON, latitude)[0])
    extended = abs(geodetic_to_sinusoidal(ANTEMERIDIAN_LON + reach, latitude)[0])

    # The claim is about the ground the row covers, not about the number
    # the last index carries. Under absorption those are two different
    # quantities: the last column stands short of the arc and the cell it
    # names stretches out to meet it. Its outer edge is the measurement.
    last_cell = f"{quadrant}({inside:04d}/{row:04d})"
    outer = max(abs(x) for x, _ in boundary.plane_ring(last_cell)[1])
    assert outer == pytest.approx(extended, abs=1.0)
    assert outer > plain
    assert (extended - plain) / L1 == pytest.approx((outer - plain) / L1, abs=1e-6)


@pytest.mark.parametrize(
    ("zone", "quadrant", "row"), [("CHUKOTKA", "NE", 750), ("FIJI", "SE", 200)]
)
def test_criterion_12_the_column_past_the_zone_limit_is_refused(
    zone: str, quadrant: str, row: int
) -> None:
    last = _last_column(quadrant, row)
    assert boundary.is_valid_cell(f"{quadrant}({last:04d}/{row:04d})") is True
    assert boundary.is_valid_cell(f"{quadrant}({last + 1:04d}/{row:04d})") is False
    assert boundary.extension_zone(f"{quadrant}({last:04d}/{row:04d})") == zone


@pytest.mark.parametrize(
    ("zone", "east", "west", "row"),
    [("CHUKOTKA", "NE", "NW", 750), ("FIJI", "SE", "SW", 200)],
)
def test_criterion_12_the_western_quadrant_gives_up_what_the_eastern_gains(
    zone: str, east: str, west: str, row: int
) -> None:
    """The extension moves a border; it does not duplicate ground.

    Inside a band the western quadrant stops at the zone's meridian
    instead of the antemeridian, so the two quadrants still meet on one
    line and no position is claimed twice.
    """
    latitude = abs(cells.cell_to_anchor(f"{east}(0001/{row:04d})")[1])
    limit = abs(EXTENSION_ZONES[zone].lon_limit)
    last_west = _last_column(west, row)
    expected = abs(geodetic_to_sinusoidal(limit, latitude)[0])
    outer = max(
        abs(x) for x, _ in boundary.plane_ring(f"{west}({last_west:04d}/{row:04d})")[1]
    )
    assert outer == pytest.approx(expected, abs=1.0)
    assert boundary.extension_zone(f"{west}({last_west:04d}/{row:04d})") == zone


# --------------------------------------------------------------------------
# Criterion 13 -- construction in latitude, and the discontinuity F6 inherits
# --------------------------------------------------------------------------


@pytest.mark.parametrize(("zone", "quadrant"), [("CHUKOTKA", "NE"), ("FIJI", "SE")])
def test_criterion_13_the_column_count_steps_at_the_band_edge(
    zone: str, quadrant: str
) -> None:
    """The same column is accepted in one row and refused in the next.

    This is the discontinuity F6 inherits: neighbouring rows hold
    different numbers of columns, so a neighbour rule that assumes a
    constant row width is wrong here and only here.
    """
    first, last = boundary.ZONE_ROWS[zone]
    inside_low = _last_column(quadrant, first)
    outside_low = _last_column(quadrant, first - 1)
    inside_high = _last_column(quadrant, last)
    outside_high = _last_column(quadrant, last + 1)

    assert inside_low > outside_low
    assert inside_high > outside_high
    probe = outside_low + 1
    assert boundary.is_valid_cell(f"{quadrant}({probe:04d}/{first:04d})") is True
    assert boundary.is_valid_cell(f"{quadrant}({probe:04d}/{first - 1:04d})") is False


@pytest.mark.parametrize("zone", ["FIJI", "CHUKOTKA"])
def test_criterion_13_the_band_is_realized_on_whole_rows(zone: str) -> None:
    """No cell is split by a latitude limit, so none needs a stepped edge.

    The declared limits are multiples of half a degree and fall inside
    rows. Realizing the zone on the rows the band touches makes its
    latitude boundary a row boundary, and it can only add ocean: the
    realized band contains the declared one.
    """
    first, last = boundary.ZONE_ROWS[zone]
    _, declared_min, _, declared_max = boundary.extension_bounds(zone)
    low = abs(geodetic_to_sinusoidal(0.0, declared_min)[1])
    high = abs(geodetic_to_sinusoidal(0.0, declared_max)[1])
    low, high = min(low, high), max(low, high)
    assert first * L1 <= low
    assert (last + 1) * L1 >= high


@pytest.mark.parametrize("zone", ["FIJI", "CHUKOTKA"])
def test_criterion_13_the_declared_band_is_covered_by_the_realized_rows(
    zone: str,
) -> None:
    """Every declared latitude lands in a row the zone owns."""
    _, declared_min, _, declared_max = boundary.extension_bounds(zone)
    quadrant = EXTENSION_ZONES[zone].quadrant
    for latitude in (declared_min, (declared_min + declared_max) / 2, declared_max):
        ordinate = abs(geodetic_to_sinusoidal(0.0, latitude)[1])
        row = int(ordinate // L1)
        first, last = boundary.ZONE_ROWS[zone]
        assert first <= row <= last
        assert boundary.extension_zone(f"{quadrant}(0001/{row:04d})") == zone


# --------------------------------------------------------------------------
# Criterion 14 -- only the border cell is a trapezoid
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("zone", "quadrant", "row"), [("CHUKOTKA", "NE", 750), ("FIJI", "SE", 200)]
)
def test_criterion_14_inside_the_zone_is_nominal_and_only_the_border_is_clipped(
    zone: str, quadrant: str, row: int
) -> None:
    """Both cells are extension cells; only one of them is a trapezoid."""
    last = _last_column(quadrant, row)
    border = f"{quadrant}({last:04d}/{row:04d})"
    interior = f"{quadrant}({last - 4:04d}/{row:04d})"

    assert boundary.is_extension_cell(border) is True
    assert boundary.is_extension_cell(interior) is True
    assert boundary.is_trapezoidal_cell(border) is True
    assert boundary.is_trapezoidal_cell(interior) is False
    assert effective_cell_area(interior) == nominal_cell_area(1)
    assert boundary.extension_zone(border) == zone


@pytest.mark.parametrize(
    ("quadrant", "row"), [("NE", 750), ("SE", 200), ("NE", 400), ("SW", 300)]
)
def test_criterion_14_a_cell_agrees_with_its_own_centroid(
    quadrant: str, row: int
) -> None:
    """``extension_zone`` on a cell equals ``extension_zone_for_point`` on
    its centroid, and the row realization is what guarantees it."""
    for column in (1, 5, _last_column(quadrant, row)):
        cell = f"{quadrant}({column:04d}/{row:04d})"
        if boundary.is_valid_cell(cell) is not True:
            continue
        centroid = cells.cell_to_centroid(cell)
        assert boundary.extension_zone(cell) == boundary.extension_zone_for_point(
            *centroid
        )


@pytest.mark.parametrize("row", [0, 100, 200, 465, 750, 900, 999])
def test_criterion_14_exactly_one_cell_per_row_is_a_trapezoid(row: int) -> None:
    """ "Only the cell that meets the border", and now it holds everywhere.

    Under clipping it did not. A cell's outer side leans one column per
    row while the border, being a meridian, leans
    ``pi * a * sin(phi) / rho(phi)`` columns per row, which passes 1 at
    about 18.5 degrees. Above that the border line crossed several
    columns within one row's height and cut every one of them, so two or
    three cells per row came back clipped.

    Absorption removes the whole effect. A cell whose upper side lies
    wholly outside the border does not exist at all, and the last one
    that does exist carries its outer side out to the border by itself.
    The count is one, from the equator to the polar row.
    """
    last = _last_column("NE", row)
    trapezoids = [
        column
        for column in range(max(last - 9, 1), last + 1)
        if boundary.cell_shape(f"NE({column:04d}/{row:04d})") == "trapezoid"
    ]
    assert trapezoids == [last]


def test_criterion_14_is_boundary_cell_is_the_coarse_screen() -> None:
    """True for every discontinuity, and false in the plain interior."""
    assert boundary.is_boundary_cell("NE(0000/0500)") is True
    assert boundary.is_boundary_cell("NE(2003/0000)") is True
    assert boundary.is_boundary_cell("NE(1400/0374)") is False
    row = 750
    assert boundary.is_boundary_cell(f"NE({_last_column('NE', row):04d}/0750)") is True


def test_absorbs_border_answers_a_quadrant_rather_than_raising() -> None:
    """It composes like every other predicate here.

    Resolution 0 has no geometry, so the underlying construction raises.
    The predicate answers ``False`` instead, because a caller screening a
    mixed list for area-safe cells should not have to guard each entry.
    """
    assert boundary.absorbs_border("NE") is False
    assert boundary.absorbs_border("NE(2003/0000)") is True
    assert boundary.absorbs_border("NE(1400/0374)") is False
    assert boundary.absorbs_border("NE(0000/1000)") is True
    assert boundary.absorbs_border("NE(1400/0374(1,2))") == [False, False]


# --------------------------------------------------------------------------
# Trapezoid refinement, and the prefix that is not the parent
# --------------------------------------------------------------------------


def _refinement_of(parent: str) -> list[str]:
    """Every existing cell one resolution below a trapezoid.

    A trapezoid reaches east past its own nominal column, so its children
    are not all spelled under its own resolution-1 prefix. The column one
    to the east has to be searched too.
    """
    quadrant, base, *rest = parent.split("(")
    column, _, row = base.partition("/")
    row = row.rstrip(")")
    candidates = [
        f"{quadrant}({int(column) + step:04d}/{row}({code}))"
        for step in (0, 1)
        for code in "1234"
    ]
    return [cell for cell in candidates if boundary.is_valid_cell(cell)]


def test_the_equator_trapezoid_refines_into_five_cells() -> None:
    """Three parallelograms, two trapezoids, and one alien prefix.

    ``NE(2003/0000)`` extends east to 2003.7508 sides, which reaches into
    ground the column 2004 would nominally hold. Its refinement therefore
    has a child whose resolution-1 component is 2004, not 2003 -- the one
    place in ITACaRT where ``prefix(child) != parent(child)``. The count
    is five, not four, and that is what makes the anomaly unavoidable
    rather than a spelling choice.
    """
    children = _refinement_of("NE(2003/0000)")
    assert children == [
        "NE(2003/0000(1))",
        "NE(2003/0000(2))",
        "NE(2003/0000(3))",
        "NE(2003/0000(4))",
        "NE(2004/0000(3))",
    ]
    shapes = [boundary.cell_shape(cell) for cell in children]
    assert shapes.count("parallelogram") == 3
    assert shapes.count("trapezoid") == 2
    assert boundary.cell_shape("NE(2003/0000(2))") == "trapezoid"
    assert boundary.cell_shape("NE(2004/0000(3))") == "trapezoid"

    # The alien prefix is a prefix and nothing more.
    assert boundary.is_valid_cell("NE(2004/0000)") is False
    assert boundary.is_valid_cell("NE(2004/0000(3))") is True


def test_the_five_children_do_not_overlap_one_another() -> None:
    """A partition needs disjointness even when the parent is irregular."""
    polygons = [plane_polygon(cell) for cell in _refinement_of("NE(2003/0000)")]
    for left in range(len(polygons)):
        for right in range(left + 1, len(polygons)):
            assert polygons[left].intersection(polygons[right]).area == pytest.approx(
                0.0, abs=1e-6
            )


def test_a_position_in_the_alien_child_quantises_back_to_it() -> None:
    """Both directions. Geometry alone would not prove the index reachable."""
    inside = plane_polygon("NE(2004/0000(3))").representative_point()
    assert cells.sinusoidal_to_cell(inside.x, inside.y, 2) == "NE(2004/0000(3))"


def test_a_trapezoid_does_not_refine_into_a_fixed_number_of_children() -> None:
    """Five is this cell's answer, not the rule.

    The children of a trapezoid are whichever cells of the finer lattice
    survive inside it, and how many that is depends on how much of the
    parent's own column the border left standing. Enumerated across
    latitude rather than argued.
    """
    counts = {}
    for row in (0, 100, 200, 465, 750, 900):
        parent = f"NE({_last_column('NE', row):04d}/{row:04d})"
        counts[row] = len(_refinement_of(parent))
    assert counts == {0: 5, 100: 2, 200: 2, 465: 5, 750: 6, 900: 4}


def test_the_children_of_a_trapezoid_exceed_it_by_the_chord_sagitta() -> None:
    """The parent is not exactly the union of its children, and cannot be.

    The domain border is a meridian, hence a curve on the plane, and
    every cell approximates it by the chord across its own height. A
    child spans less height, so its chord hugs the curve more closely and
    reaches slightly further out. The union of the children therefore
    covers a sliver the parent does not, bounded by the sagitta of the
    parent's chord.

    Closing the gap would mean giving a boundary cell one edge per
    descendant, which contradicts the four-vertex rule the paper states.
    It is recorded rather than removed, and ``F5`` and ``F7`` inherit it.
    """
    for row in (0, 100, 465, 900):
        parent = f"NE({_last_column('NE', row):04d}/{row:04d})"
        outer = unary_union([plane_polygon(c) for c in _refinement_of(parent)])
        inner = plane_polygon(parent)
        assert outer.area > inner.area
        assert outer.contains(inner.buffer(-1e-6))
        assert (outer.area - inner.area) / inner.area < 2e-3


# --------------------------------------------------------------------------
# The syntactic ceiling of column 2004
# --------------------------------------------------------------------------


def _prefix_of_fine_cell(row: int, fine_row: int, column: int, side: float) -> int:
    """Resolution-1 column a fine-lattice cell's anchor falls in."""
    sheared = (column + fine_row) * side
    return int(math.floor((sheared + 1e-6) / L1)) - row


def _rows_reaching_the_alien_prefix(
    resolution: int, rows: int = 40, every_sub_row: bool = False
) -> list[int]:
    """Rows in which some cell of one resolution is spelled under 2004.

    Only the topmost sub-row of each row is consulted by default. That is
    the one nearest the next row up, so it is where the cell's lean has
    carried it furthest east and where the alien prefix appears first if
    it appears at all;
    :func:`test_the_topmost_sub_row_is_the_one_that_reaches` checks the
    shortcut against full enumeration where full enumeration is
    affordable. At resolution 13 a row holds a million sub-rows, so the
    shortcut is what makes the deepest case testable at all.
    """
    side = cell_size(resolution)
    per_row = int(round(L1 / side))
    subs = range(per_row) if every_sub_row else (per_row - 1,)
    reaching = []
    for row in range(0, rows):
        for sub in subs:
            fine_row = row * per_row + sub
            last = boundary.last_lattice_column("NE", fine_row, side)
            if _prefix_of_fine_cell(row, fine_row, last, side) > RES1_TABLE_MAX:
                reaching.append(row)
                break
    return reaching


@pytest.mark.parametrize("resolution", [2, 3, 4])
def test_the_topmost_sub_row_is_the_one_that_reaches(resolution: int) -> None:
    """The shortcut the deeper measurements rest on, checked exhaustively."""
    assert _rows_reaching_the_alien_prefix(
        resolution, every_sub_row=True
    ) == _rows_reaching_the_alien_prefix(resolution)


def test_the_alien_prefix_reaches_further_as_the_resolution_deepens() -> None:
    """Where 2004 is actually needed, measured per resolution.

    A finer lattice hugs the border more closely, so a cell whose anchor
    falls past 2004's western edge appears in rows where the coarser
    lattice had none. The reach grows and then stops: rows 0 to 9 at
    resolution 2, 0 to 15 at 3 and 4, and 0 to 16 from resolution 5 down
    to the finest, where it settles.
    """
    assert _rows_reaching_the_alien_prefix(2) == list(range(0, 10))
    assert _rows_reaching_the_alien_prefix(3) == list(range(0, 16))
    assert _rows_reaching_the_alien_prefix(4) == list(range(0, 16))
    for resolution in (5, 7, 9, 11, 13):
        assert _rows_reaching_the_alien_prefix(resolution) == list(range(0, 17))


def test_ending_the_row_at_2003_does_not_by_itself_summon_the_prefix() -> None:
    """The stated predicate is necessary and not sufficient.

    "2004/Y is admissible if and only if row Y's last column is 2003"
    holds in one direction only. Eighteen rows end at 2003, and only
    seventeen of them ever spell a cell under 2004: row 17 ends at 2003
    and still leaves no room, at any resolution the system defines.

    Recorded rather than smoothed over, because a guard written from the
    stated predicate would admit an index that names nothing.
    """
    ends_at_2003 = [
        row
        for row in range(0, 40)
        if boundary.last_lattice_column("NE", row, L1) == RES1_TABLE_MAX
    ]
    assert ends_at_2003 == list(range(0, 18))

    deepest = _rows_reaching_the_alien_prefix(MAX_RESOLUTION)
    assert deepest == list(range(0, 17))
    assert set(deepest) < set(ends_at_2003)
    assert 17 in ends_at_2003 and 17 not in deepest

    assert boundary.is_valid_cell("NE(2004/0000(3))") is True
    assert boundary.is_valid_cell("NE(2004/0017(3))") is False


def test_no_row_anywhere_reaches_column_2005() -> None:
    """Enumerated over every row of the globe, extension zones included.

    The zones are where the eastern quadrant grows past 180 degrees, so
    they are the only place a ceiling argument could fail for a reason
    other than ``cos(phi)``. They do not come close: the widest parallel
    is the equator, outside any zone, at 2003.7508 sides.
    """
    widest = 0.0
    for quadrant in ("NE", "NW", "SE", "SW"):
        for row in range(0, 1001):
            for ordinate in (row * L1, (row + 1) * L1):
                widest = max(widest, boundary._x_border(quadrant, row, ordinate) / L1)
    assert int(widest) == RES1_TABLE_MAX
    assert widest < RES1_TABLE_MAX + 1
    assert not boundary.is_valid_cell("NE(2005/0000(3))")


RES1_TABLE_MAX = 2003
"""The last resolution-1 column Table 1 names."""


# --------------------------------------------------------------------------
# crosses_antemeridian
# --------------------------------------------------------------------------


def test_crosses_antemeridian_reads_the_geometry_as_given() -> None:
    """A wrapped ring crosses; a ring written past 180 does not.

    An extension-zone footprint is naturally written with longitudes past
    the line, and it lies on one side. What crosses is the same shape
    wrapped into ``[-180, 180]``, whose vertices then jump half a globe.
    """
    wrapped = Polygon([(179.0, 0.0), (-179.0, 0.0), (-179.0, 1.0), (179.0, 1.0)])
    unwrapped = Polygon([(179.0, 0.0), (181.0, 0.0), (181.0, 1.0), (179.0, 1.0)])
    ordinary = Polygon([(10.0, 0.0), (11.0, 0.0), (11.0, 1.0), (10.0, 1.0)])
    assert boundary.crosses_antemeridian(wrapped) is True
    assert boundary.crosses_antemeridian(unwrapped) is False
    assert boundary.crosses_antemeridian(ordinary) is False


def test_crosses_antemeridian_handles_lines_holes_and_empties() -> None:
    holed = Polygon(
        [(170.0, -5.0), (-170.0, -5.0), (-170.0, 5.0), (170.0, 5.0)],
        [[(178.0, -1.0), (179.0, -1.0), (179.0, 1.0), (178.0, 1.0)]],
    )
    assert boundary.crosses_antemeridian(holed) is True
    assert boundary.crosses_antemeridian(LineString([(179.0, 0.0), (-179.0, 0.0)]))
    assert (
        boundary.crosses_antemeridian(LineString([(10.0, 0.0), (11.0, 0.0)])) is False
    )
    assert boundary.crosses_antemeridian(Point(0.0, 0.0)) is False
    assert boundary.crosses_antemeridian(Polygon()) is False


# --------------------------------------------------------------------------
# Vectorised semantics
# --------------------------------------------------------------------------


def test_every_predicate_is_positionally_aligned_on_a_composed_index() -> None:
    """A composed index answers a list, an atomic one a bare value."""
    composed = "NE(1400/0374(1,2,3,4))"
    for function in (
        boundary.cell_shape,
        boundary.is_valid_cell,
        boundary.is_boundary_cell,
        boundary.is_triangular_cell,
        boundary.is_trapezoidal_cell,
        boundary.is_equal_area_cell,
        boundary.is_extension_cell,
        boundary.extension_zone,
    ):
        answer = function(composed)
        assert isinstance(answer, list)
        assert len(answer) == 4
        assert not isinstance(function("NE(1400/0374)"), list)


def test_a_bare_quadrant_has_no_zone_and_no_geometry() -> None:
    assert boundary.extension_zone("NE") is None
    assert boundary.cell_shape("NE") == "parallelogram"
    assert boundary.is_boundary_cell("NE") is False


def test_the_polar_row_is_clipped_rather_than_deleted() -> None:
    """The meridian quadrant holds 1000.2 cells, so the row exists in part.

    Its triangle survives, cut both by the pole and by a border that has
    come nearer the meridian than one cell width, so it is no longer an
    equal-area cell.
    """
    cell = "NE(0000/1000)"
    assert boundary.is_valid_cell(cell) is True
    assert boundary.is_equal_area_cell(cell) is False
    assert effective_cell_area(cell) < nominal_cell_area(1)
    assert max(abs(y) for _, y in boundary.plane_ring(cell)[1]) <= MERIDIAN_QUADRANT
