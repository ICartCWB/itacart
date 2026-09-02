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

# --------------------------------------------------------------------------
# The oracle
# --------------------------------------------------------------------------

#: Small enough to enumerate exhaustively, wide enough to hold the meridian,
#: the equator and the corner where they meet.
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
