"""Crossing the prime meridian and the equator, pinned rather than assumed.

Three things live here that nothing else in the suite measured.

The first is a brute-force scan. Every other neighbour test asks the
package what a cell's neighbours are and checks the answer against a
rule; this one builds the rings of every candidate in a window and asks
shapely which ones touch. That difference is the point: the module used
to decide which quadrants to search from a sentence in a docstring, and
a scan is the only instrument that can falsify a sentence.

The second is a regression floor for the equator. Filling a region that
straddles latitude zero works today and nothing named it, so a change
aimed at the meridian could break it in silence. The longitudes and the
resolutions are written out rather than drawn, because the families that
break at boundaries hold one cell per row over as many as two thousand
columns and sampling finds them by luck.

The third asks a question the suite had never asked: does the rest of
the package accept what the package itself produces? Every earlier
instrument fed a module inputs that module knew how to make.
"""

from __future__ import annotations

import pytest
from shapely.geometry import LineString, Point, Polygon, box

import itacart
from itacart import topology

# --------------------------------------------------------------------------
# The corner where both axes meet
# --------------------------------------------------------------------------

#: The two spellings of the corner triangle the oracle admits. The western
#: spellings, ``NW(0000/0000)`` and ``SW(0000/0000)``, name the same two
#: cells and ``is_valid_cell`` denies them.
CORNER_TRIANGLES = ("NE(0000/0000)", "SE(0000/0000)")

#: Cells next to the corner that reach one axis but not both. They are the
#: control: whatever the corner needs, these must not change.
NEAR_CORNER = ("NE(0001/0000)", "NW(0001/0000)", "SE(0001/0000)", "SW(0001/0000)")


def _ring(cell: str) -> Polygon:
    return Polygon(itacart.cell_to_boundary(cell))


def _scan(cell: str, reach: int = 4) -> set[str]:
    """Cells in a window whose ring actually meets this one, by measurement.

    Independent of the package's own neighbour machinery on purpose: it
    spells candidates, builds their rings and intersects. Nothing here
    consults ``get_neighbor``, ``grid_disk`` or the contact cache, so a
    fault shared by all three cannot hide from it.
    """
    here = _ring(cell)
    found = set()
    for quadrant in ("NE", "NW", "SE", "SW"):
        for column in range(reach):
            for row in range(reach):
                other = f"{quadrant}({column:04d}/{row:04d})"
                if other == cell or not itacart.is_valid_cell(other):
                    continue
                if here.intersects(_ring(other)):
                    found.add(other)
    return found


@pytest.mark.parametrize("cell", CORNER_TRIANGLES)
def test_the_corner_triangle_touches_the_opposite_quadrant(cell: str) -> None:
    """The claim that no cell reaches the diagonally opposite quadrant.

    It was written in the candidate generator's docstring and was false.
    ``SE(0000/0000)`` and ``NW(0001/0000)`` share the point at longitude
    -0.08983152841195215, latitude 0, and both rings carry it with
    identical bits, so this is a vertex the tessellation has rather than
    a sliver floating point invented.
    """
    quadrant = cell[:2]
    opposite = topology._ACROSS_MERIDIAN[topology._ACROSS_EQUATOR[quadrant]]
    touching = _scan(cell)
    assert any(other.startswith(opposite) for other in touching), (
        f"{cell} touches nothing in {opposite}, so the scan cannot see the "
        "contact this test exists to protect"
    )


@pytest.mark.parametrize("cell", CORNER_TRIANGLES + NEAR_CORNER)
def test_the_measured_contacts_are_exactly_what_the_scan_finds(cell: str) -> None:
    """The package's neighbour set, against a set built by measurement.

    Equality in both directions. A missing contact loses a neighbour in
    silence, which is how the opposite-quadrant fault survived; a surplus
    one would put a cell in a disk that does not touch it.
    """
    edges, vertices = topology._contacts(cell)
    assert set(edges) | set(vertices) == _scan(cell)


def test_the_corner_names_its_neighbour_across_both_axes() -> None:
    """A direction, not just a contact.

    Being in the contact set is not enough: ``grid_disk`` reads the set
    directly and would look complete while no direction could name the
    cell. The step from the southern corner triangle to ``NW(0001/0000)``
    crosses the meridian and the equator at once, which is the step the
    package had no way to take.
    """
    assert itacart.get_neighbor("SE(0000/0000)", "NW") == "NW(0001/0000)"
    assert itacart.get_neighbor("NE(0000/0000)", "SW") == "SW(0001/0000)"


@pytest.mark.parametrize("cell", NEAR_CORNER)
def test_the_disk_at_the_corner_is_the_same_size_in_every_quadrant(
    cell: str,
) -> None:
    """Symmetry across the four quadrants, at the radius that broke it.

    The western quadrants used to return eight cells at ``k=1`` where the
    eastern ones returned nine, and the missing one was the contact
    across both axes. The asymmetry appeared only at ``k=1``: from
    ``k=2`` the same cell was reached by another path, so a test written
    at ``k=2`` or beyond would have passed throughout.
    """
    disk = set(itacart.grid_disk(cell, 1, dedupe=True, flatten=True))
    assert len(disk) == 9


@pytest.mark.parametrize("k", range(1, 5))
def test_the_disk_is_the_closure_of_the_measured_contacts(k: int) -> None:
    """What ``grid_disk`` returns, against walking the contact sets by hand.

    This is the property ``grid_disk`` actually has. It is worth pinning
    because the obvious expectation -- that a disk of radius ``k`` holds
    ``(2k + 1) ** 2`` cells -- is a statement about a square lattice and
    is false here: near the corner the meridian column carries one
    triangle where a square lattice would carry two half cells, so the
    neighbourhood is genuinely smaller. Checking cardinality against
    ``(2k + 1) ** 2`` measures the wrong object and reports a defect that
    is not there.
    """
    cell = "NE(0001/0000)"
    seen = {cell}
    frontier = {cell}
    for _ in range(k):
        nxt: set[str] = set()
        for member in frontier:
            edges, vertices = topology._contacts(member)
            nxt |= set(edges) | set(vertices)
        frontier = nxt - seen
        seen |= nxt
    assert set(itacart.grid_disk(cell, k, dedupe=True, flatten=True)) == seen


def test_the_interior_disk_is_untouched_by_the_corner_fix() -> None:
    """The control. Away from both axes the square-lattice count holds."""
    for k in range(1, 5):
        disk = itacart.grid_disk("NE(0500/0300)", k, dedupe=True, flatten=True)
        assert len(set(disk)) == (2 * k + 1) ** 2


# --------------------------------------------------------------------------
# The equator: a floor, not a fix
# --------------------------------------------------------------------------

#: Longitudes far from both the meridian and the seam, written out so that
#: the case list cannot drift with a random seed.
EQUATOR_LONGITUDES = (-170.0, -90.0, -30.0, 30.0, 90.0, 170.0)


@pytest.mark.parametrize("longitude", EQUATOR_LONGITUDES)
@pytest.mark.parametrize("resolution", (1, 2))
def test_a_box_on_the_equator_fills_symmetrically(
    longitude: float, resolution: int
) -> None:
    """Filling across latitude zero, and the same count on each side.

    Working behaviour with no test on it is the first thing to break in
    silence when the meridian descent is rewritten, and the equator had
    none: ``grep -c equator tests/unit/test_geometry.py`` returned zero.
    """
    half = 0.5
    region = box(longitude - half, -half, longitude + half, half)
    filled = itacart.polyfill(region, resolution)
    north = sum(
        1 for cell in itacart.iter_cells(filled) if cell.startswith(("NE", "NW"))
    )
    south = sum(
        1 for cell in itacart.iter_cells(filled) if cell.startswith(("SE", "SW"))
    )
    assert north > 0 and south > 0, "the region straddles the equator"
    assert north == south


@pytest.mark.parametrize("longitude", EQUATOR_LONGITUDES)
@pytest.mark.parametrize("resolution", (1, 2))
def test_counting_agrees_with_filling_across_the_equator(
    longitude: float, resolution: int
) -> None:
    """``count_internal_cells`` against ``count_cells(polyfill(...))``.

    Two functions reaching the same number by different routes is the
    property; either one alone only says it is self-consistent.
    """
    half = 0.5
    region = box(longitude - half, -half, longitude + half, half)
    assert itacart.count_internal_cells(region, resolution) == itacart.count_cells(
        itacart.polyfill(region, resolution)
    )


def test_a_box_away_from_both_axes_is_the_control() -> None:
    """The same two checks where nothing is crossed."""
    region = box(29.5, 9.5, 30.5, 10.5)
    filled = itacart.polyfill(region, 2)
    assert all(cell.startswith("NE") for cell in itacart.iter_cells(filled))
    assert itacart.count_internal_cells(region, 2) == itacart.count_cells(filled)


# --------------------------------------------------------------------------
# Totality: does the package accept what the package produces?
# --------------------------------------------------------------------------


def _produced_indices() -> list[str]:
    """Indices built by public functions rather than written by hand.

    The distinction is the whole instrument. An index a test author types
    is an index the author already believed was representable; an index
    ``compose`` returns is one the package commits to.
    """
    produced = [
        itacart.compose(["NE(0001/0001)", "NE(0002/0001)"]),
        itacart.compose(["NE(0001/0001)", "NW(0003/0004)"]),
        itacart.compose(["NE(0001/0000)", "SE(0001/0000)"]),
        itacart.normalize(itacart.compose(["NE(0002/0001)", "NE(0001/0001)"])),
        itacart.polyfill(box(29.5, 9.5, 30.5, 10.5), 2),
        itacart.polyfill(box(29.5, -0.5, 30.5, 0.5), 1),
        itacart.recompose_to_prefix_form("NE(0001/0001)"),
        "".join(itacart.compose(list(itacart.get_children("NE(0001/0001)"))[0])),
    ]
    return produced


@pytest.mark.parametrize("index", _produced_indices())
def test_every_index_the_package_produces_survives_a_round_trip(index: str) -> None:
    """What one function emits, the others take.

    Written against the shape of ``B-8.4``: a codec refused the canonical
    output of ``compose`` for nearly a year of phases because every test
    it had fed it indices it had built itself.
    """
    assert itacart.count_cells(index) >= 1
    assert itacart.normalize(index)
    assert itacart.decompose(index)
    assert itacart.encode_tree(index)


# --------------------------------------------------------------------------
# The antimeridian, declared rather than discovered
# --------------------------------------------------------------------------

#: Cells whose ring reaches neither axis, one of them next to the seam. The
#: widened candidate window must not reach any of them.
AWAY_FROM_THE_ORIGIN = ("NE(0000/0110)", "NE(0500/0300)", "NE(1999/0000)")


@pytest.mark.parametrize("cell", AWAY_FROM_THE_ORIGIN)
def test_the_widened_window_stops_at_the_cells_that_reach_both_axes(
    cell: str,
) -> None:
    """The fourth quadrant is searched only where a ring reaches both axes.

    Scoped deliberately. The seam is not this phase's subject and its
    behaviour is not this phase's to change, so what is pinned is that
    the change cannot reach it: a cell at the last lattice column takes
    the same three-quadrant window it always took.
    """
    west, south, east, north = topology._ring(cell).bounds
    assert not (west <= 0.0 <= east and south <= 0.0 <= north)


# --------------------------------------------------------------------------
# The extension zones
# --------------------------------------------------------------------------

#: Positions the two zones govern, and one just outside each, written out
#: because the edges are where the two functions used to part company.
EXTENSION_PROBES = (
    ("FIJI interior", -179.5, -20.0),
    ("FIJI beside the dateline", -179.98, -20.0),
    ("FIJI in the extra row band", -179.5, -21.51),
    ("CHUKOTKA interior", -179.0, 66.0),
    ("CHUKOTKA beside the dateline", -179.98, 66.0),
    ("CHUKOTKA in the extra row band", -179.0, 63.95),
    ("CHUKOTKA below the band", -179.0, 63.5),
    ("east of 180, no zone needed", 179.5, -20.0),
)


@pytest.mark.parametrize(("label", "longitude", "latitude"), EXTENSION_PROBES)
def test_the_fill_resolves_an_extension_zone_the_way_the_point_does(
    label: str, longitude: float, latitude: float
) -> None:
    """``polyfill`` against ``geo_to_cell`` inside and around both zones.

    A zone lets a row's domain run past 180 degrees so that land beyond
    the line is reached from the near side. ``geo_to_cell`` honours that
    by adding 360 degrees to a position inside a zone; the fill projected
    at the literal longitude, so the same ground landed in the western
    quadrant at a column past the end of its row, and every cell it named
    there failed ``is_valid_cell``. Eighty-eight of eighty-eight, for a
    region half a degree across over Fiji.

    The probes include a latitude inside the row band but outside the
    declared degrees. The zone is realized on whole rows, so the two
    differ by up to a row, and a window drawn on the declared numbers
    would leave that strip resolving one way for the point and another
    for the fill.
    """
    target = itacart.geo_to_cell(longitude, latitude, 1)
    step = 1e-7
    region = box(longitude - step, latitude - step, longitude + step, latitude + step)
    filled = itacart.decompose(itacart.polyfill(region, 1, containment="intersects"))
    assert filled == [target]
    assert itacart.is_valid_cell(filled[0])


def test_a_region_inside_a_zone_is_filled_with_cells_that_exist() -> None:
    """The areal form, which is what the fill is for.

    A point probe cannot show that the cut is placed correctly, only that
    the shift happens. This fills half a degree of Fiji and asks the
    package whether every cell it named is a cell.
    """
    region = box(-179.8, -20.5, -179.2, -19.5)
    for resolution in (1, 2):
        cells = itacart.decompose(
            itacart.polyfill(region, resolution, containment="intersects")
        )
        assert cells
        assert all(itacart.is_valid_cell(cell) for cell in cells)
        assert all(cell.startswith("SE") for cell in cells)
        assert itacart.count_internal_cells(region, resolution) == itacart.count_cells(
            itacart.polyfill(region, resolution)
        )


@pytest.mark.parametrize("longitude", (-178.0001, -177.9999))
def test_a_position_at_a_zone_limit_is_refused_for_now(longitude: float) -> None:
    """The limit is the domain edge, and the cell there absorbs the border.

    Lifting a position at ``lon_limit`` puts it at the longitude bounding
    the eastern quadrant in those rows, so both sides of the limit sit in
    the last lattice column of their row. That column absorbs the strip
    between itself and the border and is not the sheared square the
    descent tests, so the fill refuses it.

    **This is a limitation of the fill, not of the grid.** The cell
    exists: ``geo_to_cell`` names it, ``boundary`` builds its ring, and
    the hierarchy addresses its children. The refusal stands only until
    the fill can descend a trapezoid, and the test moves to an equality
    against ``geo_to_cell`` when it can.
    """
    step = 1e-7
    region = box(longitude - step, -20.0 - step, longitude + step, -20.0 + step)
    assert itacart.is_valid_cell(itacart.geo_to_cell(longitude, -20.0, 1))
    with pytest.raises(itacart.NonExistentCellError, match="last lattice column"):
        itacart.polyfill(region, 1, containment="intersects")


def test_a_zone_limit_is_the_domain_edge_of_its_rows() -> None:
    """Why a region straddling the limit is refused rather than split.

    Lifting a position at ``lon_limit`` puts it at ``lon_limit + 360``,
    which is exactly the longitude bounding the eastern quadrant in those
    rows -- the extension limit and the domain edge are one line. So the
    lifted part of a straddling region lands in the last lattice column,
    which is a refused family for reasons that have nothing to do with
    the zone, and the refusal names that family.
    """
    from itacart.boundary import ZONE_ROWS, _lon_limit

    for name, spec in itacart.EXTENSION_ZONES.items():
        rows = ZONE_ROWS[name]
        edge = _lon_limit(spec.quadrant, rows[0])
        assert edge == spec.lon_limit + 360.0, name

    with pytest.raises(itacart.NonExistentCellError, match="last lattice column"):
        itacart.polyfill(box(-178.4, -20.4, -177.6, -19.6), 1)


def test_the_zone_window_follows_the_meridian_and_not_its_chord() -> None:
    """The cut edge is a curve on the plane, and is cut as one.

    A meridian projects to a curve, so the straight edge of a box in
    degrees would miss it by 26.6 km over Fiji's band -- more than two
    resolution-1 cells, measured. The window carries its own vertices
    along that edge; this reads the worst departure that survives.
    """
    from itacart.geodesy import geodetic_to_sinusoidal
    from itacart.geometry import _extension_windows

    worst = 0.0
    for window in _extension_windows():
        edge = max(x for x, _ in window.exterior.coords)
        along = sorted({y for x, y in window.exterior.coords if x == edge})
        assert len(along) > 1000, "the limit edge carries its own vertices"
        for first, second in zip(along, along[1:]):
            start = geodetic_to_sinusoidal(edge, first)
            stop = geodetic_to_sinusoidal(edge, second)
            for step in (0.25, 0.5, 0.75):
                here = geodetic_to_sinusoidal(edge, first + (second - first) * step)
                ratio = (here[1] - start[1]) / (stop[1] - start[1])
                chord = start[0] + ratio * (stop[0] - start[0])
                worst = max(worst, abs(here[0] - chord))
    assert worst < 0.01, f"the window edge departs from the meridian by {worst:.5f} m"


def test_land_east_of_the_dateline_is_addressed_without_crossing_it() -> None:
    """Why refusing the antimeridian is a design property and not a gap.

    An extension zone lets a row's domain run past 180 degrees so that
    the land beyond the line is reached from the near side. Without this
    written down, a later phase reads the refusal as a defect and
    "fixes" it.
    """
    fiji = itacart.geo_to_cell(179.9, -16.6, 4)
    chukotka = itacart.geo_to_cell(-179.9, 66.0, 4)
    assert itacart.is_valid_cell(fiji)
    assert itacart.is_valid_cell(chukotka)
    assert itacart.extension_zone(fiji) == "FIJI"
    assert itacart.extension_zone(chukotka) == "CHUKOTKA"


# --------------------------------------------------------------------------
# Column zero, which only the east may name
# --------------------------------------------------------------------------

#: Each western root and the eastern one that has to end up carrying its
#: column 0. Written as a pair so the southern half cannot be forgotten:
#: every fault this file records was first found in one hemisphere only.
MERIDIAN_TWINS = (("NW", "NE"), ("SW", "SE"))

#: Every row column 0 has. Measured rather than transcribed: row 1000 is
#: the polar cap and exists, row 1001 does not.
MERIDIAN_ROWS = range(1001)


@pytest.mark.parametrize(("west", "east"), MERIDIAN_TWINS)
def test_compose_folds_the_meridian_column_to_the_east(west: str, east: str) -> None:
    """The fold, in the operation the defect was reported against."""
    assert itacart.compose([f"{west}(0000/0110)"]) == f"{east}(0000/0110)"
    assert (
        itacart.compose([f"{east}(0000/0110)", f"{west}(0000/0110)"])
        == f"{east}(0000/0110)"
    )


@pytest.mark.parametrize(("west", "east"), MERIDIAN_TWINS)
def test_normalize_folds_the_meridian_column_to_the_east(west: str, east: str) -> None:
    """The same fold in the other folding operation, padding included.

    ``0/110`` and ``0000/0110`` are one cell, so the fold has to read the
    column as a number; a string comparison against ``"0000"`` would miss
    the short spelling and leave it in the west.
    """
    assert itacart.normalize(f"{west}(0000/0110)") == f"{east}(0000/0110)"
    assert itacart.normalize(f"{west}(0/110)") == f"{east}(0000/0110)"


@pytest.mark.parametrize(("west", "east"), MERIDIAN_TWINS)
def test_every_row_of_the_meridian_column_folds_to_a_cell_that_exists(
    west: str, east: str
) -> None:
    """Enumerated, not sampled, because the exceptions are one row deep.

    The polar cap sits in row 1000 alone and the corner triangle in row 0
    alone; a test drawn from the middle of the column describes neither.
    """
    for row in MERIDIAN_ROWS:
        folded = itacart.compose([f"{west}(0000/{row:04d})"])
        assert folded == f"{east}(0000/{row:04d})"
        assert itacart.is_valid_cell(folded)


@pytest.mark.parametrize(("west", "east"), MERIDIAN_TWINS)
def test_the_emptied_western_root_is_dropped_and_not_left_bare(
    west: str, east: str
) -> None:
    """The trap in moving a node between roots.

    A root with no children denotes its whole quadrant, so a western root
    that has just given away its only base cell must go, not stay. Leaving
    it would turn one triangle into a hemisphere and nothing downstream
    would object: the result parses, validates and encodes.
    """
    folded = itacart.compose([f"{west}(0000/0110)"])
    assert west not in folded
    assert itacart.count_cells(folded) == 1


@pytest.mark.parametrize(("west", "east"), MERIDIAN_TWINS)
def test_a_whole_eastern_quadrant_still_absorbs_the_folded_cell(
    west: str, east: str
) -> None:
    """Containment is decided in one place, after the fold, not twice.

    The eastern root already holds the triangle, so the folded cell adds
    nothing. This passes because the fold hands the merge a fresh root and
    lets the childless-occurrence rule settle it; a fold that appended to
    an existing bare root would shrink the region to a single cell here.
    """
    assert itacart.compose([east, f"{west}(0000/0110)"]) == east
    assert itacart.normalize(f"{east},{west}(0000/0110)") == east


@pytest.mark.parametrize(("west", "east"), MERIDIAN_TWINS)
def test_the_fold_carries_the_whole_refinement_subtree(west: str, east: str) -> None:
    """Children of the triangle are spelled from the east like their parent.

    They reach across the line -- that is what the triangle is -- but the
    resolution-1 component of that place is ``0000/Y`` on the eastern side
    whatever side a child falls on, so the subtree travels intact.
    """
    assert itacart.compose([f"{west}(0000/0110(1))"]) == f"{east}(0000/0110(1))"
    assert (
        itacart.compose([f"{east}(0000/0110(1))", f"{west}(0000/0110(2))"])
        == f"{east}(0000/0110(1,2))"
    )


@pytest.mark.parametrize(("west", "east"), MERIDIAN_TWINS)
def test_a_western_root_keeps_every_column_that_is_not_zero(
    west: str, east: str
) -> None:
    """The control. Only column 0 moves, and the rest of the root survives."""
    folded = itacart.compose([f"{west}(0000/0110,0001/0110)"])
    assert set(itacart.decompose(folded)) == {
        f"{east}(0000/0110)",
        f"{west}(0001/0110)",
    }
    assert itacart.compose([f"{west}(0001/0110)"]) == f"{west}(0001/0110)"


@pytest.mark.parametrize(("west", "east"), MERIDIAN_TWINS)
def test_the_fold_does_not_reach_the_layer_that_judges_existence(
    west: str, east: str
) -> None:
    """Where the fold deliberately stops, pinned so it is not extended.

    ``is_valid_cell`` reads the quadrant off what ``iter_cells`` yields,
    so a fold in the parse layer would answer True for a cell the codec
    refuses -- the predicate would stop measuring existence and start
    measuring its own rewrite. Folding belongs to the operations that
    emit an index, not to the ones that report what they were given.
    """
    spelling = f"{west}(0000/0110)"
    assert itacart.is_valid_cell(spelling) is False
    assert itacart.decompose(spelling) == [spelling]
    assert list(itacart.iter_cells(spelling)) == [spelling]


def test_composing_what_decompose_produced_is_still_an_identity() -> None:
    """The fold must not move anything that was already addressable."""
    for index in (
        "NE(0000/0110)",
        "SE(0000/0000)",
        "NE(0000/0110(1,2))",
        "NW(0001/0110),NE(0000/0110)",
        "NE(0500/0300),SW(0002/0004)",
    ):
        assert itacart.compose(itacart.decompose(index)) == index


# --------------------------------------------------------------------------
# Who owns a position exactly on an axis
# --------------------------------------------------------------------------

#: Longitudes far from the meridian and the seam, reused from the block
#: above so the two sections cannot drift apart.
AXIS_LONGITUDES = (-170.0, -30.0, 30.0, 170.0)


@pytest.mark.parametrize("longitude", AXIS_LONGITUDES)
def test_a_point_on_the_equator_is_filled_once_and_to_the_north(
    longitude: float,
) -> None:
    """``polyfill`` and ``geo_to_cell`` must name the same owner.

    They disagreed. Shapely can only clip with a closed box, so all four
    quadrant windows carried their own boundary and a point on the axis
    landed in two of them; ``polyfill`` returned both the northern and the
    southern cell where ``geo_to_cell`` returns one. This is the smallest
    form the fault takes: one point, one cell, two answers.
    """
    point = Point(longitude, 0.0)
    filled = itacart.polyfill(point, 1, containment="intersects")
    assert itacart.decompose(filled) == [itacart.geo_to_cell(longitude, 0.0, 1)]


@pytest.mark.parametrize("longitude", AXIS_LONGITUDES)
def test_a_line_lying_on_the_equator_is_filled_once(longitude: float) -> None:
    """The same rule for a feature with extent rather than a single point."""
    line = LineString([(longitude - 0.5, 0.0), (longitude + 0.5, 0.0)])
    filled = itacart.polyfill(line, 1, containment="intersects")
    quadrants = {cell[:2] for cell in itacart.iter_cells(filled)}
    assert quadrants == {"NE" if longitude > 0 else "NW"}


@pytest.mark.parametrize("longitude", AXIS_LONGITUDES)
def test_a_line_crossing_the_equator_still_reaches_both_sides(
    longitude: float,
) -> None:
    """The control. Disowning the axis must not disown the crossing.

    A rule that dropped the far side whenever a feature met the axis would
    pass the two tests above and be wrong here, so the length is checked
    as well as the quadrants: the southern part has to survive intact.
    """
    line = LineString([(longitude, -0.5), (longitude, 0.5)])
    filled = itacart.polyfill(line, 1, containment="intersects")
    quadrants = {cell[:2] for cell in itacart.iter_cells(filled)}
    assert len(quadrants) == 2
    north = sum(1 for c in itacart.iter_cells(filled) if c.startswith(("NE", "NW")))
    south = sum(1 for c in itacart.iter_cells(filled) if c.startswith(("SE", "SW")))
    assert north == south


@pytest.mark.parametrize("longitude", AXIS_LONGITUDES)
def test_a_point_just_off_the_equator_keeps_its_own_side(longitude: float) -> None:
    """The other control. The convention applies to the axis, not near it."""
    for latitude, hemisphere in ((1e-9, "N"), (-1e-9, "S")):
        filled = itacart.polyfill(
            Point(longitude, latitude), 1, containment="intersects"
        )
        assert all(c[0] == hemisphere for c in itacart.iter_cells(filled))


# --------------------------------------------------------------------------
# A polygon smaller than the cell that holds it
# --------------------------------------------------------------------------

#: A square about two centimetres across, far inside any cell at any
#: resolution the package addresses.
INSIDE_ONE_CELL = 1e-7

#: Positions in all four quadrants and near both axes, written out because
#: the interesting failures are the ones next to a discontinuity.
SITES = (
    (30.0, 10.0),
    (-30.0, 10.0),
    (30.0, -10.0),
    (-30.0, -10.0),
    (0.05, 10.0),
    (-0.05, 10.0),
    (30.0, 0.05),
    (30.0, -0.05),
)


@pytest.mark.parametrize("resolution", (1, 3, 5))
def test_a_polygon_inside_one_cell_names_that_cell(resolution: int) -> None:
    """``polyfill`` against ``geo_to_cell`` on the smallest region there is.

    A polygon wholly inside one cell has exactly one right answer, and
    ``geo_to_cell`` already knows it, so this compares two routes to the
    same fact rather than checking one against a rule written by hand.
    """
    for longitude, latitude in SITES:
        region = box(
            longitude - INSIDE_ONE_CELL,
            latitude - INSIDE_ONE_CELL,
            longitude + INSIDE_ONE_CELL,
            latitude + INSIDE_ONE_CELL,
        )
        filled = itacart.polyfill(region, resolution, containment="intersects")
        assert itacart.decompose(filled) == [
            itacart.geo_to_cell(longitude, latitude, resolution)
        ]


@pytest.mark.parametrize("containment", ("center", "contains"))
def test_a_polygon_inside_one_cell_covers_nothing_under_the_other_modes(
    containment: str,
) -> None:
    """The measured consequence of the two modes, pinned rather than fixed.

    Under ``center`` a cell is kept when its centre falls inside, and a
    region smaller than a cell holds no centre; under ``contains`` it holds
    no whole cell either. Both are the documented meaning, and both mean
    that ``polyfill`` is not a generalisation of ``geo_to_cell`` outside
    ``intersects``. Recorded here so that a later phase changing it knows
    it is changing something, and that ``center`` is the default.
    """
    region = box(
        30.0 - INSIDE_ONE_CELL,
        10.0 - INSIDE_ONE_CELL,
        30.0 + INSIDE_ONE_CELL,
        10.0 + INSIDE_ONE_CELL,
    )
    with pytest.raises(itacart.GeometryError, match="covers no cell"):
        itacart.polyfill(region, 3, containment=containment)


# --------------------------------------------------------------------------
# Filling the prime-meridian column
# --------------------------------------------------------------------------

#: Latitudes reaching both hemispheres and both refinement steps, written
#: out because the seam's exceptions sit in single rows.
SEAM_LATITUDES = (0.2, 9.5, 44.5, -0.2, -9.5, -44.5)


@pytest.mark.parametrize("latitude", SEAM_LATITUDES)
def test_a_region_astride_the_meridian_fills_the_column(latitude: float) -> None:
    """Column zero is filled, and every cell it emits exists.

    The fill used to refuse here. What makes the refusal removable is
    that a triangle straddling the line has one spelling, the eastern
    one, so the two halves of the region do not each name a cell.
    """
    region = box(-0.3, latitude, 0.3, latitude + 0.4)
    cells = itacart.decompose(itacart.polyfill(region, 1, containment="intersects"))
    seam = [cell for cell in cells if cell[3:7] == "0000"]
    assert seam
    assert all(itacart.is_valid_cell(cell) for cell in cells)
    assert all(cell.startswith(("NE", "SE")) for cell in seam)
    assert len(set(cells)) == len(cells)


@pytest.mark.parametrize("latitude", SEAM_LATITUDES)
def test_the_meridian_fill_agrees_with_geo_to_cell(latitude: float) -> None:
    """A point on the line, against the function that already owns it.

    ``geo_to_cell`` awards a meridian position to one cell. Until the
    column could be filled this was untestable through ``polyfill``, so
    the agreement is checked now that both routes exist.
    """
    for resolution in (1, 2, 3):
        filled = itacart.polyfill(
            Point(0.0, latitude), resolution, containment="intersects"
        )
        assert itacart.decompose(filled) == [
            itacart.geo_to_cell(0.0, latitude, resolution)
        ]


@pytest.mark.parametrize("latitude", SEAM_LATITUDES)
def test_the_containment_modes_still_nest_across_the_meridian(
    latitude: float,
) -> None:
    """``contains`` inside ``center`` inside ``intersects``, on triangles.

    The chain holds on a square because the centre is a point of the
    square. It has to be shown again here rather than assumed: the
    triangle is a different figure, its centroid sits exactly on the
    line the quadrant windows cut along, and the three modes are decided
    against the geometry before that cut.
    """
    region = box(-0.3, latitude, 0.3, latitude + 0.4)
    sets = {
        mode: set(itacart.decompose(itacart.polyfill(region, 1, containment=mode)))
        for mode in ("intersects", "center", "contains")
    }
    assert sets["contains"] <= sets["center"] <= sets["intersects"]
    assert any(cell[3:7] == "0000" for cell in sets["contains"])


@pytest.mark.parametrize("latitude", SEAM_LATITUDES)
def test_counting_agrees_with_filling_across_the_meridian(latitude: float) -> None:
    """``_count_meridian_node`` against ``_fill_meridian_node``.

    Two descents of the same triangles that must reach the same number,
    one naming cells and one not. Written for the meridian because the
    equator already had it and the two paths diverge only here.
    """
    region = box(-0.3, latitude, 0.3, latitude + 0.4)
    for resolution in (1, 2):
        assert itacart.count_internal_cells(region, resolution) == itacart.count_cells(
            itacart.polyfill(region, resolution)
        )


@pytest.mark.parametrize("latitude", SEAM_LATITUDES)
def test_every_cell_the_meridian_fill_emits_meets_the_region(
    latitude: float,
) -> None:
    """The cells are measured against the region, not against the descent.

    ``intersects`` must emit nothing the region misses and ``contains``
    nothing the region fails to hold, checked by rebuilding each cell's
    own ring rather than by trusting the coordinates the walk carried.
    """
    from shapely.geometry import Polygon as _Polygon

    from itacart.geodesy import geodetic_to_sinusoidal

    region = box(-0.3, latitude, 0.3, latitude + 0.4)
    plane = _Polygon([geodetic_to_sinusoidal(x, y) for x, y in region.exterior.coords])

    def ring(cell: str) -> _Polygon:
        return _Polygon(
            [geodetic_to_sinusoidal(x, y) for x, y in itacart.cell_to_boundary(cell)]
        )

    touching = itacart.decompose(itacart.polyfill(region, 1, containment="intersects"))
    assert all(ring(cell).intersects(plane) for cell in touching)
    inside = itacart.decompose(itacart.polyfill(region, 1, containment="contains"))
    assert all(plane.contains(ring(cell)) for cell in inside)


def test_a_region_wholly_west_of_the_line_still_reaches_the_column() -> None:
    """The half of the seam cell that only the west can see.

    A region that ends at the meridian covers the western half of the
    triangle and none of the eastern one, so no eastern piece survives
    the clip. The column still has to be filled, and from the east,
    which is why the walk is driven by the unsplit plane rather than by
    the pieces the split produced.
    """
    region = box(-0.4, 9.5, 0.0, 10.0)
    cells = itacart.decompose(itacart.polyfill(region, 1, containment="intersects"))
    seam = [cell for cell in cells if cell[3:7] == "0000"]
    assert seam
    assert all(cell.startswith("NE") for cell in seam)
    assert all(itacart.is_valid_cell(cell) for cell in cells)


def test_counting_a_region_west_of_the_line_still_counts_the_column() -> None:
    """The counting twin, where only one side of the fold has a lattice.

    A region ending at the meridian leaves no eastern piece, so the
    eastern parallelogram children of the seam have nowhere to be
    measured and contribute nothing. They must be skipped rather than
    looked up, and the count still has to match the fill.
    """
    region = box(-0.4, 9.5, 0.0, 10.0)
    for resolution in (1, 2):
        assert itacart.count_internal_cells(region, resolution) == itacart.count_cells(
            itacart.polyfill(region, resolution)
        )


def test_the_origin_is_filled_once_and_from_the_north() -> None:
    """The one position both discontinuities claim.

    The seam walk is driven by the unsplit plane, so it does not pass
    through the clip that enforces the half-open rule and has to apply
    that rule itself. Without it the point at the origin is filled twice,
    once from each hemisphere, while ``geo_to_cell`` names one cell.

    The southern control is a nanodegree away, which is the smallest
    statement that the rule applies to the axis and not near it.
    """
    assert itacart.decompose(
        itacart.polyfill(Point(0.0, 0.0), 1, containment="intersects")
    ) == [itacart.geo_to_cell(0.0, 0.0, 1)]
    assert itacart.decompose(
        itacart.polyfill(Point(0.0, -1e-9), 1, containment="intersects")
    ) == [itacart.geo_to_cell(0.0, -1e-9, 1)]


def test_the_seam_children_come_out_in_index_order() -> None:
    """One walk, not two fills reconciled afterwards.

    The fold sends ``(i, j)`` and ``(j, i)`` to opposite sides of the
    line, so a descent that filled each side separately would have to
    merge and re-sort. Walking the grid in alphabet order emits the
    children already ordered, and this is what says so.
    """
    region = box(-0.3, 9.5, 0.3, 9.9)
    filled = itacart.polyfill(region, 2, containment="intersects")
    for root in filled.split(","):
        if not root.startswith(("NE", "SE")):
            continue
        assert "0000/" in root
        assert root.index("0000/") < len(root)
    cells = [c for c in itacart.decompose(filled) if c[3:7] == "0000"]
    assert cells == sorted(cells)


def test_the_western_spelling_of_the_triangle_does_not_reach_the_codec() -> None:
    """What the folding operations emit, the codec takes.

    Was an ``xfail``. The western spelling is still not a cell, which is
    the second assertion; what changed is that ``compose`` no longer
    carries it through to output, so the pair round-trips.
    """
    east, west = "NE(0000/0110)", "NW(0000/0110)"
    assert itacart.is_valid_cell(east)
    assert not itacart.is_valid_cell(west)
    itacart.encode_tree(itacart.compose([east, west]))
