"""The rebuilt lattice distance, pinned against a measured oracle.

Nothing here asserts a number a test author believed. The oracle is
breadth-first search over ``topology._touching``, which is the same step
relation ``grid_disk`` and ``are_neighbor_cells`` already answer to, so a
disagreement is a real disagreement between two parts of the package
rather than between the package and an opinion about it.

Three families of defect are named below because each one survived a
phase by being invisible to the instruments then in use: refinement codes
of a triangle read as a square grid, the eastern shear applied west of the
meridian, and the quadrant prefix taken for the side of the line.
"""

from __future__ import annotations

from collections import deque

import pytest

import itacart
from itacart import topology
from itacart.exceptions import DomainError

# --------------------------------------------------------------------------
# The oracle
# --------------------------------------------------------------------------

#: The domain every count below enumerates over, written out because a count
#: without its domain cannot be checked. Columns one to five and rows zero to
#: five of each quadrant, at the resolution named by the test.
#:
#: What this domain excludes is as much a part of it: the last addressable
#: column of a row, and the rows nearest the pole. Both are measured and both
#: are wrong, pinned below rather than left to a handoff.
WINDOW = 5


def _spell(quadrant: str, column: int, row: int) -> str:
    return f"{quadrant}({column:04d}/{row:04d})"


def _neighbourhood(quadrants: tuple[str, ...], reach: int) -> list[str]:
    cells = [
        _spell(quadrant, column, row)
        for quadrant in quadrants
        for column in range(reach + 1)
        for row in range(reach + 1)
    ]
    return [cell for cell in cells if itacart.is_valid_cell(cell)]


def _oracle(origin: str, metric: str, reach: int) -> dict[str, int]:
    """Shortest path lengths by search over the measured contact sets."""
    seen = {origin: 0}
    queue = deque([origin])
    while queue:
        cell = queue.popleft()
        if seen[cell] >= reach:
            continue
        for other in topology._touching(cell, metric):
            if other not in seen:
                seen[other] = seen[cell] + 1
                queue.append(other)
    return seen


# --------------------------------------------------------------------------
# What a distance is
# --------------------------------------------------------------------------

SPREAD = (
    "NE(0500/0300)",
    "NE(0000/0110)",
    "NE(0000/0110(3))",
    "NE(0000/0110(1(E5)))",
    "NW(0001/0000)",
    "SE(0000/0000)",
    "SW(0002/0004)",
    "NE(0500/0300(1(C3(2))))",
)


@pytest.mark.parametrize("cell", SPREAD)
@pytest.mark.parametrize("metric", ("chebyshev", "manhattan"))
def test_a_cell_is_no_distance_from_itself(cell: str, metric: str) -> None:
    assert itacart.grid_distance(cell, cell, metric) == 0


@pytest.mark.parametrize("metric", ("chebyshev", "manhattan"))
def test_the_measure_reads_the_same_in_both_directions(metric: str) -> None:
    """Symmetry, including across the line where the basis mirrors.

    The crossing is resolved by minimising over the row at which the line
    is crossed, and the two legs enter that sum in opposite order. An
    asymmetry would mean the minimisation depended on which end was named
    first, which is how an earlier crossing reported different numbers for
    a pair depending on where the origin lay.
    """
    cells = _neighbourhood(("NE", "NW", "SE", "SW"), 3)
    for origin in cells:
        for destination in cells:
            assert itacart.grid_distance(
                origin, destination, metric
            ) == itacart.grid_distance(destination, origin, metric)


def test_an_edge_neighbour_is_one_step_under_either_metric() -> None:
    """``are_neighbor_cells`` is edge adjacency, and an edge is one step.

    Both metrics, because an edge is a step in both. The vertex case is the
    opposite implication and is pinned separately below.
    """
    cells = _neighbourhood(("NE", "NW", "SE", "SW"), 3)
    checked = 0
    for origin in cells:
        for destination in cells:
            if not itacart.are_neighbor_cells(origin, destination):
                continue
            checked += 1
            assert itacart.grid_distance(origin, destination, "chebyshev") == 1
            assert itacart.grid_distance(origin, destination, "manhattan") == 1
    assert checked > 0, "the window holds no edge neighbours to check"


def test_a_shared_vertex_is_a_step_only_under_chebyshev() -> None:
    """Why the apex crossing exists in one metric and not the other.

    The two cells that meet at a meridian triangle's apex share a single
    point. Charging one Manhattan step for it would report a path that
    cannot be walked, so the crossing there has to go through the triangle
    and the distance is two.
    """
    east, west = "NE(0001/0300)", "NW(0001/0300)"
    assert not itacart.are_neighbor_cells(east, west)
    assert itacart.grid_distance(east, west, "chebyshev") == 1
    assert itacart.grid_distance(east, west, "manhattan") == 2


# --------------------------------------------------------------------------
# Against the oracle
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "origin",
    (
        "NE(0003/0003)",
        "NE(0000/0000)",
        "NE(0001/0000)",
        "SE(0000/0000)",
        "NW(0001/0000)",
        "SW(0002/0000)",
        "NE(0000/0003)",
    ),
)
@pytest.mark.parametrize("metric", ("chebyshev", "manhattan"))
def test_the_measure_agrees_with_search_at_resolution_one(
    origin: str, metric: str
) -> None:
    reach = 4 if metric == "chebyshev" else 3
    for cell, steps in _oracle(origin, metric, reach).items():
        if itacart.get_resolution(cell) != 1:
            continue
        assert itacart.grid_distance(origin, cell, metric) == steps


@pytest.mark.parametrize(
    "origin",
    (
        "NE(0000/0110(1))",
        "NE(0000/0110(3))",
        "NE(0001/0000(1))",
        "SE(0000/0000(4))",
        "NW(0001/0000(2))",
    ),
)
@pytest.mark.parametrize("metric", ("chebyshev", "manhattan"))
def test_the_measure_agrees_with_search_at_resolution_two(
    origin: str, metric: str
) -> None:
    for cell, steps in _oracle(origin, metric, 3).items():
        if itacart.get_resolution(cell) != 2:
            continue
        assert itacart.grid_distance(origin, cell, metric) == steps


@pytest.mark.parametrize(
    "origin",
    (
        "NE(0000/0110(1(A1)))",
        "NE(0000/0110(1(E5)))",
        "NW(0001/0000(2(C3)))",
    ),
)
@pytest.mark.parametrize("metric", ("chebyshev", "manhattan"))
def test_the_measure_agrees_with_search_at_resolution_three(
    origin: str, metric: str
) -> None:
    for cell, steps in _oracle(origin, metric, 2).items():
        if itacart.get_resolution(cell) != 3:
            continue
        assert itacart.grid_distance(origin, cell, metric) == steps


@pytest.mark.parametrize("metric", ("chebyshev", "manhattan"))
def test_the_measure_never_reports_a_path_shorter_than_one_that_exists(
    metric: str,
) -> None:
    """The direction of any residual error, pinned as a property.

    Excess is a path that exists and is not the shortest; a deficit is a
    path that cannot be walked. The rebuilt measure has neither in the
    window below, and this test exists so that a later change cannot trade
    one for the other in silence.
    """
    for origin in ("NE(0001/0000)", "NE(0000/0110(3))", "SE(0000/0000)"):
        for cell, steps in _oracle(origin, metric, 3).items():
            if itacart.get_resolution(cell) != itacart.get_resolution(origin):
                continue
            assert itacart.grid_distance(origin, cell, metric) >= steps


# --------------------------------------------------------------------------
# The three named defects
# --------------------------------------------------------------------------


def test_two_children_of_one_triangle_that_touch_are_one_step_apart() -> None:
    """A triangle's codes do not lay out as a square refinement grid.

    Read as a five-by-five grid these two land four columns and five rows
    apart and the measure returned nine steps under Manhattan. They share a
    boundary. The rows of a triangle's twenty-five children hold nine,
    seven, five, three and one.
    """
    east, west = "NE(0000/0110(1(E5)))", "NE(0000/0110(3(E1)))"
    assert itacart.grid_distance(east, west, "chebyshev") == 1
    assert itacart.are_neighbor_cells(east, west)


def test_a_child_west_of_the_meridian_takes_the_western_shear() -> None:
    """Orientation is geometric and the prefix does not carry it.

    ``NE(0000/0110(3))`` has an eastern prefix and lies wholly west of the
    line. Its own children run westward as their digit rises, and the
    poleward step moves them east, which is the mirror of the eastern rule.
    """
    first, last = "NE(0000/0110(3(A1)))", "NE(0000/0110(3(A5)))"
    assert itacart.grid_distance(first, last, "chebyshev") == 4
    assert topology._lattice_descent(first)[0] > topology._lattice_descent(last)[0]


def test_a_shared_prefix_does_not_mean_a_shared_side_of_the_line() -> None:
    """The crossing is decided by the descended column, not the quadrant.

    Both cells spell ``NE`` and they face each other across the meridian.
    Deciding from the prefix that no line separates them collapsed the two
    sides onto one point and returned zero steps between distinct cells.
    """
    east, west = "NE(0000/0110(2))", "NE(0000/0110(3))"
    assert topology._lattice_descent(east)[0] > 0
    assert topology._lattice_descent(west)[0] < 0
    assert itacart.grid_distance(east, west, "chebyshev") == 1


# --------------------------------------------------------------------------
# What the path may not touch
# --------------------------------------------------------------------------

#: Every door out of arithmetic that ``topology`` still holds open. Ring
#: construction, the shapely predicates the contact sets are built from, and
#: the projection the replaced hinge sampled. The module imports several of
#: these for its other functions, so their presence proves nothing about this
#: path; only firing them does.
GEOMETRY_DOORS = (
    ("itacart.topology", "cell_to_boundary"),
    ("itacart.topology", "_ring"),
    ("itacart.topology", "_contacts"),
    ("itacart.topology", "_touching"),
    ("itacart.topology", "Polygon"),
    ("itacart.topology", "LineString"),
    ("itacart.topology", "snap"),
    ("itacart.topology", "translate"),
    ("itacart.cells", "cell_to_sinusoidal"),
    ("itacart.cells", "sinusoidal_to_cell"),
)


def _arm(doors: tuple[tuple[str, str], ...], patch: pytest.MonkeyPatch) -> None:
    """Replace every named door with something that fails if it is called."""
    import importlib

    def detonator(module_name: str, label: str) -> object:
        def fired(*args: object, **kwargs: object) -> object:
            raise AssertionError(f"the distance reached {module_name}.{label}")

        return fired

    for module_name, label in doors:
        module = importlib.import_module(module_name)
        patch.setattr(module, label, detonator(module_name, label))


@pytest.mark.parametrize("metric", ("chebyshev", "manhattan"))
@pytest.mark.parametrize(
    "origin,destination",
    (
        ("NE(0500/0300)", "NE(0505/0304)"),
        ("NE(0000/0110(3))", "NE(0000/0110(2))"),
        ("NE(0003/0000)", "SW(0001/0002)"),
        ("NE(0000/0110(1(E5)))", "NW(0002/0110(1(A1)))"),
        ("SE(0000/0000(4))", "SW(0002/0004(1))"),
    ),
)
def test_the_distance_reaches_no_geometry_and_no_contact_set(
    origin: str,
    destination: str,
    metric: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The property is about the path, so the evidence has to be too.

    That two names are absent from the module says nothing: ``topology``
    still imports shapely and builds rings for its neighbour machinery, and
    a distance that quietly reached one of them would leave the imports
    looking exactly the same. Arming the doors and running the measure is
    the only instrument that can tell the two apart.

    The pairs span the four cases the measure distinguishes -- within one
    side, across the meridian under a shared prefix, across both axes, and
    across the meridian at depth -- so no branch of the path escapes by
    being untaken.
    """
    _arm(GEOMETRY_DOORS, monkeypatch)
    assert itacart.grid_distance(origin, destination, metric) >= 0


def test_the_armed_doors_are_doors(monkeypatch: pytest.MonkeyPatch) -> None:
    """The control. An armed door that no longer exists arms nothing.

    Without this, deleting a name from ``topology`` would silently shrink
    the guarantee above while leaving it green.
    """
    import importlib

    for module_name, label in GEOMETRY_DOORS:
        assert hasattr(
            importlib.import_module(module_name), label
        ), f"{module_name}.{label} is armed but absent"
    _arm(GEOMETRY_DOORS, monkeypatch)
    with pytest.raises(AssertionError, match="the distance reached"):
        itacart.are_neighbor_cells("NE(0500/0300)", "NE(0501/0300)")


# --------------------------------------------------------------------------
# The minimisation
# --------------------------------------------------------------------------


def test_the_search_finds_what_enumeration_finds() -> None:
    """Ternary search against every value, on a convex cost.

    Each crossing cost is a sum of maxima of absolute affine functions of
    the crossing row, so it is convex and its breakpoints are integers.
    This checks the consequence rather than the derivation.
    """
    for shift in range(-4, 5):

        def cost(row: int, shift: int = shift) -> int:
            return abs(row - 3) + max(abs(row + shift), abs(row - 7))

        assert topology._minimise(cost, -20, 20) == min(
            cost(row) for row in range(-20, 21)
        )


def test_the_crossing_row_is_searched_in_each_hemisphere_separately() -> None:
    """The seam bends at the equator, so each branch is bounded on its own.

    A pair far to the north leaves the southern interval with nothing in
    range, and it is clamped to the nearest row that exists rather than
    left empty.
    """
    north, south = topology._rows_to_search((900, 900), (905, 903), 1)
    assert north[0] <= north[1]
    assert south[0] <= south[1]
    assert south[1] <= -1


@pytest.mark.parametrize(
    "resolution,rows", ((1, 1000), (2, 2000), (3, 10000), (4, 20000))
)
def test_the_row_limit_follows_the_refinement(resolution: int, rows: int) -> None:
    assert topology._seam_row_limit(resolution) == rows


# --------------------------------------------------------------------------
# The border, declared rather than omitted
# --------------------------------------------------------------------------

#: Trapezoids of the last addressable column, one per row measured. The
#: measure refuses them: the row above narrows, so the poleward neighbour of
#: a last-column cell lies several columns in and the uniform shear stops
#: describing it. Land beyond the line is addressed through an extension zone
#: rather than by stepping across it.
REFUSED_ROWS = (0, 300, 600, 900)


@pytest.mark.parametrize("row", REFUSED_ROWS)
@pytest.mark.parametrize("metric", ("chebyshev", "manhattan"))
def test_a_trapezoid_is_refused_at_the_entrance(row: int, metric: str) -> None:
    trapezoid = f"NE({topology._last_column('NE', row):04d}/{row:04d})"
    assert itacart.cell_shape(trapezoid) == "trapezoid"
    with pytest.raises(DomainError, match="trapezoid"):
        itacart.grid_distance(trapezoid, "NE(0500/0300)", metric)
    with pytest.raises(DomainError, match="trapezoid"):
        itacart.grid_distance("NE(0500/0300)", trapezoid, metric)


def test_the_refusal_is_a_limit_of_the_measure_and_not_of_the_adjacency() -> None:
    """The disk still answers where the distance declines to.

    Written the same way as the refusal it replaces: a measure that stops
    where the lattice stops is not a tessellation that stops there. The
    contact sets carry the trapezoid and ``grid_disk`` reads them.
    """
    trapezoid = f"NE({topology._last_column('NE', 600):04d}/0600)"
    disk = itacart.grid_disk(trapezoid, 1, dedupe=True, flatten=True)
    assert len(set(disk)) > 1
    assert any(itacart.are_neighbor_cells(trapezoid, cell) for cell in disk)


#: Pairs the trapezoid refusal does not see -- ordinary parallelograms two
#: columns inside the last one -- that are nonetheless closer around the far
#: side of the lattice. Measured under the oracle that skips trapezoids, so
#: the short route is not one that walks the border cells: at the edges of
#: the FIJI and CHUKOTKA extension zones the seam contact needs no trapezoid
#: at all. Refused rather than answered, because reaching the land beyond the
#: line is what an extension zone is for.
ANTIMERIDIAN_REFUSED = (
    ("NE(1173/0600)", "NW(1173/0600)", "chebyshev"),
    ("NE(1173/0600)", "NW(1172/0600)", "chebyshev"),
    ("SE(1931/0171)", "SW(1931/0170)", "chebyshev"),
    ("SE(1864/0237)", "SW(1863/0238)", "manhattan"),
)


@pytest.mark.parametrize("origin,destination,metric", ANTIMERIDIAN_REFUSED)
def test_a_pair_closer_around_the_far_side_is_refused(
    origin: str, destination: str, metric: str
) -> None:
    """The test is a comparison, not a proximity.

    Both cells are ordinary parallelograms and neither is near the prime
    meridian, so no earlier guard sees them. What decides is that the route
    out through each cell's own last column costs less than the reading
    through the prime meridian, which for these pairs runs into the
    thousands against an oracle distance of two or three.
    """
    assert itacart.cell_shape(origin) == "parallelogram"
    assert itacart.cell_shape(destination) == "parallelogram"
    with pytest.raises(DomainError, match="antimeridian seam"):
        itacart.grid_distance(origin, destination, metric)


def test_the_antimeridian_refusal_spares_pairs_near_the_prime_meridian() -> None:
    """The control. A refusal by proximity would have taken these too."""
    assert itacart.grid_distance("NE(0500/0300)", "NW(0500/0300)") == 500
    assert itacart.grid_distance("NE(0001/0000)", "NW(0001/0001)", "manhattan") == 2


# --------------------------------------------------------------------------
# A decision whose boundary has to fail out loud
# --------------------------------------------------------------------------

#: Trapezoids whose refinement grid is not square, with the child count each
#: one actually has. Enumerated across four rows rather than sampled: row 300
#: has three children and none spelled elsewhere, so a test written from one
#: row would have described the wrong thing.
TRAPEZOID_ROWS = ((0, 5, 1), (300, 3, 0), (600, 4, 1), (900, 4, 1))


@pytest.mark.parametrize("row,children,foreign", TRAPEZOID_ROWS)
def test_a_trapezoid_carries_the_lattice_but_not_the_refinement_grid(
    row: int, children: int, foreign: int
) -> None:
    """Where "a trapezoid is an ordinary parallelogram" stops being true.

    It is true of position: a trapezoid's children sit exactly where the
    parallelogram recursion puts them. It is false of the refinement grid,
    and this test fails if anyone extends the decision that far. A trapezoid
    keeps only the codes its clipped area still reaches, so the child count
    is not the square of the side, and on three of these four rows one child
    is spelled under the *next* resolution-1 column.
    """
    trapezoid = f"NE({topology._last_column('NE', row):04d}/{row:04d})"
    assert itacart.cell_shape(trapezoid) == "trapezoid"

    found = list(itacart.get_children(trapezoid))[0]
    assert len(found) == children

    outside = [
        child
        for child in found
        if itacart.index.split_components(child)[:-1]
        != itacart.index.split_components(trapezoid)
    ]
    assert len(outside) == foreign
    assert (len(found) != topology._grid_side(2) ** 2) or foreign, (
        "a trapezoid whose children were a full square grid all spelled "
        "under it would make the decision safe to extend, and none is"
    )


# --------------------------------------------------------------------------
# The cell that closes a hemisphere
# --------------------------------------------------------------------------


@pytest.mark.parametrize("cap", ("NE(0000/1000)", "SE(0000/1000)"))
def test_the_polar_cap_has_one_way_out(cap: str) -> None:
    """Why the cap is reduced rather than measured.

    Three cells touch it and two are trapezoids of the last column, which
    this measure does not step on. One contact remains, so the distance
    from the cap is one more than the distance from that cell -- a cut of
    the graph at a single vertex, not a route chosen among several.
    """
    contacts = topology._touching(cap, "chebyshev")
    free = [c for c in contacts if itacart.cell_shape(c) != "trapezoid"]
    assert len(contacts) == 3
    assert free == [topology._polar_equatorward_triangle(cap)]
    assert itacart.cell_shape(free[0]) == "triangle"


@pytest.mark.parametrize(
    "cap,other,steps",
    (
        ("NE(0000/1000)", "NE(0003/0998)", 3),
        ("NE(0000/1000)", "NW(0003/0998)", 3),
        ("SE(0000/1000)", "SE(0003/0998)", 3),
        ("SE(0000/1000)", "SW(0003/0998)", 3),
        ("NE(0000/1000)", "NE(0000/0999)", 1),
    ),
)
def test_the_reduction_reads_the_cap_as_search_does(
    cap: str, other: str, steps: int
) -> None:
    """Both directions, against the oracle that skips trapezoids.

    Read without the reduction these came out one short, which is the
    direction that names a path nobody can walk.
    """
    assert itacart.grid_distance(cap, other, "chebyshev") == steps
    assert itacart.grid_distance(other, cap, "chebyshev") == steps


@pytest.mark.parametrize("cap", ("NE(0000/1000)", "SE(0000/1000)"))
def test_no_manhattan_path_leaves_the_polar_cap(cap: str) -> None:
    """The one contact is a shared vertex, and a vertex is not a step.

    Returning a number here would name a walk that does not exist. The cap
    is still zero steps from itself, which is not a path.
    """
    assert not [
        c
        for c in topology._touching(cap, "manhattan")
        if itacart.cell_shape(c) != "trapezoid"
    ]
    with pytest.raises(DomainError, match="closes its hemisphere"):
        itacart.grid_distance(cap, "NE(0001/0998)", "manhattan")
    assert itacart.grid_distance(cap, cap, "manhattan") == 0


#: What the reduction does not reach, pinned as numbers. Refining the cap
#: destroys the property the reduction rests on: at resolution 2 it has three
#: contacts free of trapezoids and at resolution 3 it has eight, so there is
#: no single vertex to cut and the generic path measures it instead.
DEEP_CAP_RESIDUE = (
    ("NE(0000/1000(1))", "NE(0002/0999(1))", 3, 2),
    ("NE(0000/1000(1))", "NW(0002/0999(1))", 3, 2),
    ("SE(0000/1000(1))", "SE(0002/0999(1))", 3, 2),
    ("NE(0000/1000(1(A1)))", "NE(0000/1000(1(A4)))", 2, 3),
    ("NE(0000/1000(1(A1)))", "NE(0000/1000(1(D1)))", 2, 3),
)


@pytest.mark.parametrize("origin,destination,true,measured", DEEP_CAP_RESIDUE)
def test_a_refined_cap_is_still_read_by_the_generic_path(
    origin: str, destination: str, true: int, measured: int
) -> None:
    """Pinned so that closing it turns this red, in both directions.

    Note the sign changes with depth: resolution 2 reads one short and
    resolution 3 one long. Two defects, not one, and a fix aimed at either
    alone will move this test.
    """
    assert itacart.grid_distance(origin, destination, "chebyshev") == measured
    assert true != measured


@pytest.mark.parametrize("resolution", (2, 3))
def test_a_refined_cap_is_not_treated_as_a_cap(resolution: int) -> None:
    """The restriction to resolution 1 is the measurement, not a shortcut."""
    cap = "NE(0000/1000)"
    for _ in range(resolution - 1):
        cap = [
            c
            for c in list(itacart.get_children(cap))[0]
            if itacart.index.split_components(c)[:-1]
            == itacart.index.split_components(cap)
        ][0]
    assert itacart.get_resolution(cap) == resolution
    assert not topology._is_polar_cap(cap)
    assert topology._is_polar_cap("NE(0000/1000)")
