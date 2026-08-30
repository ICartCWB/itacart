"""Tests for :mod:`itacart.topology`.

The lattice step is pinned against shared edges, computed with shapely from
the cell rings, and never against anchors or against symmetry. Both weaker
instruments accept the wrong rule: the cell diagonally above is symmetric too,
and its anchor differs from the true neighbour's only in the component the
error already got right.
"""

from __future__ import annotations

import pytest
from shapely.affinity import translate
from shapely.geometry import Polygon
from shapely.ops import snap

import itacart
from itacart import topology
from itacart.boundary import ZONE_ROWS, last_lattice_column
from itacart.exceptions import DomainError, InvalidIndexError, ResolutionError
from itacart.hierarchy import _parent_cell, _render
from itacart.index import split_components

BASE_SIDE_M = 10_000.0
SNAP_TOLERANCE = 1e-9

# The seam is written as +-180 but stored a unit in the last place either side
# of it, so an exact intersection test reports a gap of about six nanometres on
# roughly a third of the rows. Snapping first makes the comparison a statement
# about the tessellation rather than about binary64.


def _polygon(cell: str, longitude_offset: float = 0.0) -> Polygon:
    ring = Polygon(itacart.cell_to_boundary(cell))
    return translate(ring, xoff=longitude_offset) if longitude_offset else ring


def shared_edge_length(first: str, second: str) -> float:
    """Length of the boundary two cells hold in common, zero when none.

    Tries the mirrored quadrant shifted by a full turn as well, so a pair that
    meets across the antemeridian is measured rather than missed.
    """
    left = _polygon(first)
    for offset in (0.0, 360.0, -360.0):
        right = _polygon(second, offset)
        if left.distance(right) > 1.0:
            continue
        overlap = left.intersection(snap(right, left, SNAP_TOLERANCE))
        if overlap.length > 1e-8:
            return float(overlap.length)
    return 0.0


def share_an_edge(first: str, second: str) -> bool:
    return shared_edge_length(first, second) > 0.0


def share_only_a_vertex(first: str, second: str) -> bool:
    left, right = _polygon(first), _polygon(second)
    overlap = left.intersection(snap(right, left, SNAP_TOLERANCE))
    return (not overlap.is_empty) and overlap.length == 0.0


def last_column(quadrant: str, row: int) -> int:
    return last_lattice_column(quadrant, row, BASE_SIDE_M)


INTERIOR = "NE(0500/0300)"


# --------------------------------------------------------------------------
# Criterion 1: the lattice step
# --------------------------------------------------------------------------


def test_north_step_moves_one_column_toward_the_meridian() -> None:
    """The falsifier for the paper's item (a) and for its first amendment.

    The paper decrements ``Y`` for the cell above. The amendment corrected the
    sense of ``Y`` and dropped the term entirely. Both spell a cell that meets
    the origin at one vertex; only ``(X - 1, Y + 1)`` shares an edge.
    """
    assert topology.get_neighbor(INTERIOR, "N") == "NE(0499/0301)"
    assert share_an_edge(INTERIOR, "NE(0499/0301)")
    assert not share_an_edge(INTERIOR, "NE(0500/0301)")
    assert share_only_a_vertex(INTERIOR, "NE(0500/0301)")


def test_south_step_moves_one_column_away_from_the_meridian() -> None:
    assert topology.get_neighbor(INTERIOR, "S") == "NE(0501/0299)"
    assert share_an_edge(INTERIOR, "NE(0501/0299)")
    assert not share_an_edge(INTERIOR, "NE(0500/0299)")


@pytest.mark.parametrize(
    ("direction", "expected"),
    [
        ("N", "NE(0499/0301)"),
        ("S", "NE(0501/0299)"),
        ("E", "NE(0501/0300)"),
        ("W", "NE(0499/0300)"),
    ],
)
def test_cardinal_neighbours_share_an_edge(direction: str, expected: str) -> None:
    assert topology.get_neighbor(INTERIOR, direction) == expected
    assert share_an_edge(INTERIOR, expected)


@pytest.mark.parametrize(
    ("direction", "expected"),
    [
        ("NE", "NE(0500/0301)"),
        ("NW", "NE(0498/0301)"),
        ("SE", "NE(0502/0299)"),
        ("SW", "NE(0500/0299)"),
    ],
)
def test_diagonal_neighbours_share_only_a_vertex(direction: str, expected: str) -> None:
    """Criterion 6 rests on this: a diagonal is not an edge neighbour."""
    assert topology.get_neighbor(INTERIOR, direction) == expected
    assert share_only_a_vertex(INTERIOR, expected)


@pytest.mark.parametrize("quadrant", ["NE", "NW", "SE", "SW"])
def test_step_arithmetic_is_identical_in_every_quadrant(quadrant: str) -> None:
    """Mirroring is in the geography, not in the index arithmetic."""
    origin = f"{quadrant}(0500/0300)"
    for direction, (shift_x, shift_y) in topology.LATTICE_STEP.items():
        expected = f"{quadrant}({500 + shift_x:04d}/{300 + shift_y:04d})"
        assert topology.get_neighbor(origin, direction) == expected


@pytest.mark.slow
def test_edge_neighbours_account_for_the_whole_perimeter() -> None:
    """An ordinary row, enumerated rather than sampled, with no exceptions.

    If the named neighbours cover the perimeter exactly, no further cell can
    be edge-adjacent. That is a stronger statement than finding four, and it
    is what makes the neighbour set a fact about the row instead of a fact
    about the search window.

    Every cell closes, the triangle at the meridian and the trapezoid at the
    seam included. They used to be left open because the lexical step could
    name only one cell per side and their sides carry more than one.
    """
    row = 300
    unclosed: list[tuple[str, float]] = []
    for column in range(last_column("NE", row) + 1):
        cell = f"NE({column:04d}/{row:04d})"
        perimeter = _polygon(cell).exterior.length
        covered = sum(
            shared_edge_length(cell, neighbour)
            for neighbour in topology._edge_neighbors(cell)
        )
        if abs(perimeter - covered) > 1e-9:
            unclosed.append((cell, perimeter - covered))
    assert unclosed == []


# --------------------------------------------------------------------------
# Criterion 9: derivation is lexical, existence is a separate question
# --------------------------------------------------------------------------


def test_derivation_never_consults_geometry(monkeypatch: pytest.MonkeyPatch) -> None:
    """Proved by construction, not by inspection.

    Every function the module could use to look at a coordinate is replaced by
    one that raises. The derivation still answers, so it cannot have been
    calling them.
    """

    def refuse(*args: object, **kwargs: object) -> object:
        raise AssertionError("the derivation consulted geometry")

    monkeypatch.setattr(itacart.cells, "cell_to_anchor", refuse)
    monkeypatch.setattr(itacart.cells, "cell_to_boundary", refuse)
    monkeypatch.setattr(itacart.cells, "cell_to_centroid", refuse)
    monkeypatch.setattr(itacart.geodesy, "geodetic_to_sinusoidal", refuse)
    monkeypatch.setattr(itacart.geodesy, "inverse_geodesic", refuse)

    target, _, _ = topology._lexical_target(INTERIOR, "N")
    assert target == "NE(0499/0301)"


def test_derivation_survives_without_the_existence_predicates(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The other half of the same property.

    ``is_valid_cell`` and ``last_lattice_column`` answer whether a neighbour
    exists. They must not be answering *where* it is.
    """

    def refuse(*args: object, **kwargs: object) -> object:
        raise AssertionError("the derivation consulted an existence predicate")

    monkeypatch.setattr(topology, "is_valid_cell", refuse)
    monkeypatch.setattr(topology, "last_lattice_column", refuse)

    for direction, (shift_x, shift_y) in topology.LATTICE_STEP.items():
        target, target_x, target_y = topology._lexical_target(INTERIOR, direction)
        assert target == f"NE({500 + shift_x:04d}/{300 + shift_y:04d})"
        assert (target_x, target_y) == (500 + shift_x, 300 + shift_y)


# --------------------------------------------------------------------------
# Criterion 5: the boundary deflections
# --------------------------------------------------------------------------


def test_crossing_the_prime_meridian_from_the_eastern_triangle() -> None:
    assert topology.get_neighbor("NE(0000/0300)", "W") == "NW(0001/0300)"
    assert share_an_edge("NE(0000/0300)", "NW(0001/0300)")


def test_crossing_the_prime_meridian_into_the_eastern_triangle() -> None:
    """Western quadrants have no column zero, so the step lands in the east."""
    assert topology.get_neighbor("NW(0001/0300)", "W") == "NE(0000/0300)"
    assert share_an_edge("NW(0001/0300)", "NE(0000/0300)")


def test_poleward_from_a_western_column_one_lands_on_the_meridian_triangle() -> None:
    """The same deflection, reached by a vertical step rather than a horizontal one."""
    assert topology.get_neighbor("NW(0001/0300)", "N") == "NE(0000/0301)"
    assert share_an_edge("NW(0001/0300)", "NE(0000/0301)")


def test_the_meridian_triangle_has_no_poleward_neighbour() -> None:
    """Its apex meets the triangle above at one point, so there is no edge."""
    assert topology.get_neighbor("NE(0000/0300)", "N") is None
    assert not share_an_edge("NE(0000/0300)", "NE(0000/0301)")


def test_the_meridian_triangle_keeps_its_equatorward_step() -> None:
    assert topology.get_neighbor("NE(0000/0300)", "S") == "NE(0001/0299)"
    assert share_an_edge("NE(0000/0300)", "NE(0001/0299)")


@pytest.mark.parametrize(
    ("origin", "expected"),
    [
        ("NE(0500/0000)", "SE(0500/0000)"),
        ("SE(0500/0000)", "NE(0500/0000)"),
        ("NW(0500/0000)", "SW(0500/0000)"),
        ("SW(0500/0000)", "NW(0500/0000)"),
    ],
)
def test_crossing_the_equator_keeps_the_column(origin: str, expected: str) -> None:
    """The hemispheres meet column to column; the shear does not carry across."""
    assert topology.get_neighbor(origin, "S") == expected
    assert share_an_edge(origin, expected)


@pytest.mark.parametrize("quadrant", ["NE", "NW", "SE", "SW"])
def test_the_polar_cell_is_reached_from_either_quadrant(quadrant: str) -> None:
    """One polar cell per hemisphere, in the eastern quadrant, spanning 360."""
    origin = f"{quadrant}(0001/0999)"
    expected = f"{quadrant[0]}E(0000/1000)"
    assert topology.get_neighbor(origin, "N") == expected
    assert share_an_edge(origin, expected)


def test_the_antemeridian_seam_outside_an_extension_zone() -> None:
    row = 300
    east = f"NE({last_column('NE', row):04d}/{row:04d})"
    west = f"NW({last_column('NW', row):04d}/{row:04d})"
    assert topology.get_neighbor(east, "E") == west
    assert topology.get_neighbor(west, "E") == east
    assert share_an_edge(east, west)


def test_the_antemeridian_seam_inside_an_extension_zone() -> None:
    """The two last columns differ here, so the seam cannot be a longitude test."""
    row = ZONE_ROWS["CHUKOTKA"][0]
    east = f"NE({last_column('NE', row):04d}/{row:04d})"
    west = f"NW({last_column('NW', row):04d}/{row:04d})"
    assert (east, west) == ("NE(0933/0709)", "NW(0830/0709)")
    assert topology.get_neighbor(east, "E") == west
    assert topology.get_neighbor(west, "E") == east
    assert share_an_edge(east, west)


# --------------------------------------------------------------------------
# The geometric exception set
# --------------------------------------------------------------------------


def test_zone_limit_rows_are_the_four_rows_of_each_zone_boundary() -> None:
    expected = frozenset({170, 171, 237, 238, 708, 709, 799, 800})
    assert topology.ZONE_LIMIT_ROWS == expected


def test_poleward_overflow_is_resolved_rather_than_refused() -> None:
    """The row above is shorter, so the lexical target names no cell.

    The neighbour is real and sits two columns west, not one, and the offset
    does not follow from the two row limits. The module measures it instead of
    guessing, and instead of refusing as it did while the fallback was
    missing.
    """
    cell = f"NE({last_column('NE', 300):04d}/0300)"
    neighbour = topology.get_neighbor(cell, "N")
    assert neighbour == "NE(1782/0301)"
    assert share_an_edge(cell, neighbour)


def test_a_zone_latitude_limit_sends_the_step_into_the_other_quadrant() -> None:
    """Where a zone starts, one quadrant is extended and its mirror truncated.

    The cell poleward of this one is therefore not in its own quadrant at all,
    which no lexical rule can name and which neither of the two outcomes the
    opening package predicted would have produced.
    """
    neighbour = topology.get_neighbor("NW(0883/0708)", "N")
    assert neighbour == "NE(0884/0709)"
    assert share_an_edge("NW(0883/0708)", neighbour)


# --------------------------------------------------------------------------
# Vectorised semantics and input validation
# --------------------------------------------------------------------------


def test_a_composite_index_returns_one_answer_per_cell() -> None:
    composite = "NE(0500/0300,0501/0300)"
    assert topology.get_neighbor(composite, "E") == [
        "NE(0501/0300)",
        "NE(0502/0300)",
    ]


def test_a_single_cell_returns_a_scalar() -> None:
    assert isinstance(topology.get_neighbor(INTERIOR, "E"), str)


def test_an_unknown_direction_is_refused() -> None:
    with pytest.raises(ValueError, match="not a lattice direction"):
        topology.get_neighbor(INTERIOR, "up")  # type: ignore[arg-type]


def test_a_refined_cell_resolves_through_its_parent() -> None:
    """Resolution 2 and deeper go through the refinement grid, not the pair."""
    assert topology.get_neighbor("NE(0500/0300(1))", "N") == "NE(0500/0300(3))"


def test_a_quadrant_has_no_lattice_pair() -> None:
    with pytest.raises(ResolutionError, match="not a resolution-1 cell"):
        topology._res1_parts("NE")


def test_a_pair_without_a_separator_is_refused() -> None:
    """Rejected upstream, by the index grammar, before topology sees it."""
    with pytest.raises(InvalidIndexError, match="must be 'X/Y'"):
        topology._res1_parts("NE(0500)")


def test_deflect_normalises_a_diagonal_across_two_boundaries() -> None:
    """The case the one-dimensional dispatch could not express.

    This step leaves through the prime meridian and through the row above at
    once. Deciding on one overflow at a time answered ``None``; normalising
    the target through both frames answers the cell.
    """
    assert topology.deflect("NE(0000/0300)", "NW") == "NW(0002/0301)"


def test_the_polar_cell_has_no_poleward_neighbour() -> None:
    assert topology.deflect("NE(0000/1000)", "N") is None


def test_deflect_is_the_identity_on_a_step_that_crosses_nothing() -> None:
    """Normalisation of a target already in range leaves it alone.

    Worth pinning: it is what makes ``deflect`` safe to call unconditionally,
    which is how the module now uses it.
    """
    assert topology.deflect(INTERIOR, "N") == "NE(0499/0301)"


# --------------------------------------------------------------------------
# Criteria 2, 3 and 4: descent through the refinement alphabets
# --------------------------------------------------------------------------

QUINARY_PARENT = "NE(0500/0300(1))"


def child(parent: str, code: str) -> str:
    return _render(split_components(parent) + [code])


@pytest.mark.parametrize(
    ("code", "direction", "expected"),
    [
        ("A1", "E", "NE(0500/0300(1(A2)))"),
        ("A1", "N", "NE(0500/0300(1(B1)))"),
        ("C3", "W", "NE(0500/0300(1(C2)))"),
        ("C3", "S", "NE(0500/0300(1(B3)))"),
    ],
)
def test_a_step_inside_a_quinary_parent_keeps_the_prefix(
    code: str, direction: str, expected: str
) -> None:
    """Letter is the row, digit is the column, and the grid is square.

    The parent's own frame shears with the lattice, so the poleward step keeps
    the column here. Expressed in absolute resolution-1 columns the same step
    moves one west; the two statements are the same geometry in two frames.
    """
    assert topology.get_neighbor(child(QUINARY_PARENT, code), direction) == expected
    assert share_an_edge(child(QUINARY_PARENT, code), expected)


@pytest.mark.parametrize(
    ("code", "direction", "expected"),
    [
        ("1", "E", "NE(0500/0300(2))"),
        ("1", "N", "NE(0500/0300(3))"),
        ("4", "W", "NE(0500/0300(3))"),
        ("4", "S", "NE(0500/0300(2))"),
    ],
)
def test_a_step_inside_a_quaternary_parent(
    code: str, direction: str, expected: str
) -> None:
    """The quad is 1 south-west, 2 south-east, 3 north-west, 4 north-east.

    Vertical adjacency shifts by two and horizontal by one, as the paper says,
    but the horizontal shift needs the row guard: the eastern step off code 2
    leaves the parent rather than landing on code 3.
    """
    origin = f"NE(0500/0300({code}))"
    assert topology.get_neighbor(origin, direction) == expected
    assert share_an_edge(origin, expected)


def test_the_eastern_step_off_child_two_leaves_the_parent() -> None:
    """The falsification the amendment to criterion 2 asked for."""
    neighbour = topology.get_neighbor("NE(0500/0300(2))", "E")
    assert neighbour == "NE(0501/0300(1))"
    assert neighbour != "NE(0500/0300(3))"
    assert share_an_edge("NE(0500/0300(2))", neighbour)


@pytest.mark.parametrize(
    ("code", "direction", "expected"),
    [
        ("A1", "W", "NE(0499/0300(2(A5)))"),
        ("A1", "S", "NE(0501/0299(3(E1)))"),
        ("E1", "N", "NE(0500/0300(3(A1)))"),
        ("E5", "E", "NE(0500/0300(2(E1)))"),
    ],
)
def test_leaving_the_parent_wraps_to_the_far_side(
    code: str, direction: str, expected: str
) -> None:
    """Criterion 4, and the reason criterion 3's wrap-around is not a wrap.

    Joining column five to column one of the *same* parent would return a cell
    four columns away. The wrap is to the neighbouring parent, and reaching it
    can take more than one level of ascent: the equatorward step off ``A1``
    walks up twice before it lands.
    """
    origin = child(QUINARY_PARENT, code)
    assert topology.get_neighbor(origin, direction) == expected
    assert share_an_edge(origin, expected)


@pytest.mark.slow
def test_every_step_of_a_whole_refinement_agrees_with_the_geometry() -> None:
    """One resolution-1 cell refined to resolution 3: 100 cells, 8 directions.

    Enumerated, not sampled. Cardinal steps must share an edge and diagonal
    steps exactly one vertex; either failing would mean the index arithmetic
    and the tessellation have parted company.
    """
    cells = list(itacart.get_children("NE(0700/0420)", target_res=3, flatten=True))
    assert len(cells) == 100
    edges = vertices = 0
    for cell in cells:
        for direction in ("N", "S", "E", "W"):
            neighbour = topology.get_neighbor(cell, direction)
            assert isinstance(neighbour, str)
            assert share_an_edge(cell, neighbour)
            edges += 1
        for direction in ("NE", "NW", "SE", "SW"):
            neighbour = topology.get_neighbor(cell, direction)
            assert isinstance(neighbour, str)
            assert share_only_a_vertex(cell, neighbour)
            vertices += 1
    assert (edges, vertices) == (400, 400)


def test_opposite_steps_return_to_the_origin_in_the_interior() -> None:
    origin = child(QUINARY_PARENT, "C3")
    for outward, back in (("N", "S"), ("E", "W"), ("NE", "SW"), ("NW", "SE")):
        away = topology.get_neighbor(origin, outward)
        assert isinstance(away, str)
        assert topology.get_neighbor(away, back) == origin


def test_a_triangular_parent_falls_through_to_measurement() -> None:
    """Its twenty-five codes pack into rows of nine, seven, five, three, one.

    That packing is a bijection, not a square grid, so column-and-row
    arithmetic would name the wrong child rather than fail loudly. The module
    detects that the index can no longer represent the step and measures.
    """
    neighbour = topology.get_neighbor("NE(0000/0300(1))", "E")
    assert isinstance(neighbour, str)
    assert share_an_edge("NE(0000/0300(1))", neighbour)


def test_an_absorbing_parent_falls_through_to_measurement() -> None:
    """A trapezoid keeps only the codes its clipped area still reaches."""
    row = ZONE_ROWS["CHUKOTKA"][0]
    trapezoid = f"NE({last_column('NE', row):04d}/{row:04d})"
    origin = child(trapezoid, "1")
    neighbour = topology.get_neighbor(origin, "E")
    assert isinstance(neighbour, str)
    assert share_an_edge(origin, neighbour)


def test_a_quadrant_has_no_neighbour() -> None:
    with pytest.raises(ResolutionError, match="quadrant has no lattice neighbour"):
        topology.get_neighbor("NE", "N")


def test_a_refined_diagonal_across_the_meridian_resolves() -> None:
    """The case that exposed the unverified neighbour set.

    This cell touches ``NW(0001/0300(3(E1)))`` at its north-west corner,
    across the prime meridian and one level down. While the index route was
    verified in ``get_neighbor`` and not in ``_contacts``, the two disagreed:
    the step was rejected and the fallback searched a set built from the
    unverified route, so the contact went unreported and looked like a
    property of the meridian triangle. It was a property of the code.
    """
    origin = "NE(0001/0300(3(E1)))"
    touching = "NW(0001/0300(3(E1)))"
    assert share_only_a_vertex(origin, touching)
    assert topology.get_neighbor(origin, "NW") == touching
    assert touching in topology._vertex_neighbors(origin)
    assert touching in topology.grid_disk(origin, 1)

    perimeter = _polygon(origin).exterior.length
    covered = sum(
        shared_edge_length(origin, other) for other in topology._edge_neighbors(origin)
    )
    assert abs(perimeter - covered) < 1e-9


def test_a_child_beside_a_triangle_is_measured_not_wrapped() -> None:
    """Wrapping the code into a triangle would name the wrong child silently.

    The triangle holds the same twenty-five codes in a different arrangement,
    so ``A5`` exists there — it is simply not the cell to the west. The
    descent detects that and hands over rather than answering wrongly.
    """
    origin = "NE(0001/0300(1(A1)))"
    neighbour = topology.get_neighbor(origin, "W")
    assert isinstance(neighbour, str)
    assert share_an_edge(origin, neighbour)


# --------------------------------------------------------------------------
# Criteria 6 and 7: edge adjacency and disk cardinality
# --------------------------------------------------------------------------


def test_edge_adjacency_excludes_vertex_contact() -> None:
    """Criterion 6. The diagonal touches, but touching is not adjacency."""
    assert topology.are_neighbor_cells(INTERIOR, "NE(0499/0301)")
    assert not topology.are_neighbor_cells(INTERIOR, "NE(0500/0301)")
    assert share_only_a_vertex(INTERIOR, "NE(0500/0301)")


def test_adjacency_is_answered_from_both_sides() -> None:
    """A trapezoid names one of its three neighbours on a side; they name it."""
    row = 300
    trapezoid = f"NE({last_column('NE', row):04d}/{row:04d})"
    below = topology.get_neighbor(trapezoid, "S")
    assert isinstance(below, str)
    assert topology.are_neighbor_cells(below, trapezoid)
    assert topology.are_neighbor_cells(trapezoid, below)


def test_adjacency_needs_one_resolution() -> None:
    with pytest.raises(ResolutionError, match="one resolution"):
        topology.are_neighbor_cells(INTERIOR, "NE(0500/0300(1))")


def test_unit_disks_hold_nine_and_five_in_the_interior() -> None:
    """Criterion 7, stated for the interior as the criterion itself is."""
    assert len(topology.grid_disk(INTERIOR, 1)) == 9
    assert len(topology.grid_disk(INTERIOR, 1, "manhattan")) == 5
    assert len(topology.grid_ring(INTERIOR, 1)) == 8
    assert len(topology.grid_ring(INTERIOR, 1, "manhattan")) == 4


def test_the_same_cardinalities_hold_at_a_refined_resolution() -> None:
    refined = child("NE(0700/0420(1))", "C3")
    assert len(topology.grid_disk(refined, 1)) == 9
    assert len(topology.grid_disk(refined, 1, "manhattan")) == 5


def test_a_disk_of_radius_zero_is_the_origin_alone() -> None:
    assert topology.grid_disk(INTERIOR, 0) == [INTERIOR]
    assert topology.grid_ring(INTERIOR, 0) == [INTERIOR]


def test_disk_and_ring_agree_shell_by_shell() -> None:
    disk = set(topology.grid_disk(INTERIOR, 3))
    shells: set[str] = set()
    for k in range(4):
        shells |= set(topology.grid_ring(INTERIOR, k))
    assert disk == shells


def test_a_negative_radius_is_refused() -> None:
    with pytest.raises(ValueError, match="must not be negative"):
        topology.grid_disk(INTERIOR, -1)


def test_an_unknown_metric_is_refused() -> None:
    with pytest.raises(ValueError, match="not a lattice metric"):
        topology.grid_disk(INTERIOR, 1, "euclidean")  # type: ignore[arg-type]


def test_several_origins_give_one_list_each() -> None:
    composite = "NE(0500/0300,0800/0400)"
    result = topology.grid_disk(composite, 1)
    assert isinstance(result, list) and len(result) == 2
    assert all(len(group) == 9 for group in result)


def test_flatten_merges_the_lists() -> None:
    composite = "NE(0500/0300,0501/0300)"
    flat = topology.grid_disk(composite, 1, flatten=True)
    assert isinstance(flat[0], str)
    assert flat == sorted(set(flat))


def test_dedupe_drops_the_overlap() -> None:
    """Two adjacent origins share cells; dedupe keeps only what is unique."""
    composite = "NE(0500/0300,0501/0300)"
    kept = topology.grid_disk(composite, 1, dedupe=True)
    plain = topology.grid_disk(composite, 1)
    assert isinstance(kept, list) and isinstance(plain, list)
    assert all(len(k) < len(p) for k, p in zip(kept, plain))
    assert set(kept[0]).isdisjoint(kept[1])


# --------------------------------------------------------------------------
# Criterion 8: symmetry
# --------------------------------------------------------------------------


@pytest.mark.slow
def test_membership_of_a_unit_disk_is_symmetric() -> None:
    """Encouraged to falsify, and it does not fall in this form.

    Adjacency is symmetric by construction and `grid_disk` composes steps into
    a *set*, so set membership cannot disagree with itself. What is not
    single-valued is `get_neighbor` on the exception set, where one side of a
    trapezoid carries up to three cells; that is a statement about directions,
    not about disks.
    """
    cells = list(itacart.get_children("NE(0700/0420)", target_res=3, flatten=True))
    asymmetric = []
    for cell in cells:
        for other in topology.grid_disk(cell, 1):
            if cell not in topology.grid_disk(other, 1):
                asymmetric.append((cell, other))
    assert asymmetric == []


# --------------------------------------------------------------------------
# grid_distance
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("destination", "chebyshev", "manhattan"),
    [
        ("NE(0500/0300)", 0, 0),
        ("NE(0499/0301)", 1, 1),
        ("NE(0501/0300)", 1, 1),
        ("NE(0500/0301)", 1, 2),
        ("NE(0505/0303)", 8, 11),
    ],
)
def test_grid_distance_counts_steps(
    destination: str, chebyshev: int, manhattan: int
) -> None:
    """Counted against a step count, never against `d(a, b) == d(b, a)`.

    A symmetric test cannot see an error that is wrong by the same amount from
    both sides, and the shear produces exactly that kind of error.
    """
    assert topology.grid_distance(INTERIOR, destination) == chebyshev
    assert topology.grid_distance(INTERIOR, destination, "manhattan") == manhattan


@pytest.mark.slow
def test_grid_distance_agrees_with_the_disk_it_would_take_to_reach() -> None:
    """The destination is in the disk of that radius and not in the smaller one."""
    for shift_x in range(-4, 5):
        for shift_y in range(-4, 5):
            destination = f"NE({500 + shift_x:04d}/{300 + shift_y:04d})"
            for metric in ("chebyshev", "manhattan"):
                distance = topology.grid_distance(INTERIOR, destination, metric)
                assert destination in topology.grid_disk(INTERIOR, distance, metric)
                if distance:
                    assert destination not in topology.grid_disk(
                        INTERIOR, distance - 1, metric
                    )


def test_grid_distance_is_refused_across_quadrants() -> None:
    """The declared restriction, and the reason for it.

    The two cells are adjacent and the disk finds them, so this is a limit of
    the measure, not of the adjacency.
    """
    with pytest.raises(DomainError, match="different quadrants"):
        topology.grid_distance("NE(0000/0300)", "NW(0001/0300)")
    assert topology.are_neighbor_cells("NE(0000/0300)", "NW(0001/0300)")
    assert "NW(0001/0300)" in topology.grid_disk("NE(0000/0300)", 1, "manhattan")


def test_grid_distance_needs_one_resolution() -> None:
    with pytest.raises(ResolutionError, match="cells of one resolution"):
        topology.grid_distance(INTERIOR, "NE(0500/0300(1))")


@pytest.mark.parametrize(
    ("origin", "destination", "chebyshev", "manhattan"),
    [
        ("NE(0500/0300(1))", "NE(0500/0300(3))", 1, 1),
        ("NE(0500/0300(3))", "NE(0499/0301(1))", 1, 1),
    ],
)
def test_grid_distance_counts_steps_below_resolution_one(
    origin: str, destination: str, chebyshev: int, manhattan: int
) -> None:
    """The deep lattice, including a step that leaves its resolution-1 cell.

    The shear sits at the top level only, so a poleward step keeps the column
    inside a parent and moves a whole resolution-1 column when it leaves. Both
    must cost one, and the second case is the one a naive deep coordinate gets
    wrong: it would charge the width of a resolution-1 cell.
    """
    assert topology.grid_distance(origin, destination) == chebyshev
    assert topology.grid_distance(origin, destination, "manhattan") == manhattan


@pytest.mark.slow
def test_deep_grid_distance_agrees_with_the_disk_over_a_whole_refinement() -> None:
    cells = list(itacart.get_children("NE(0700/0420)", target_res=3, flatten=True))
    origin = cells[len(cells) // 2]
    for cell in cells:
        for metric in ("chebyshev", "manhattan"):
            distance = topology.grid_distance(origin, cell, metric)
            assert cell in topology.grid_disk(origin, distance, metric)


def test_a_quadrant_has_no_position_on_the_lattice() -> None:
    with pytest.raises(ResolutionError, match="no position on the lattice"):
        topology.grid_distance("NE", "NW")


def test_grid_distance_refuses_an_unknown_metric() -> None:
    with pytest.raises(ValueError, match="not a lattice metric"):
        topology.grid_distance(
            INTERIOR, INTERIOR, "euclidean"  # type: ignore[arg-type]
        )


# --------------------------------------------------------------------------
# Directed edges
# --------------------------------------------------------------------------


def test_a_directed_edge_round_trips() -> None:
    edge = topology.cells_to_directed_edge(INTERIOR, "NE(0499/0301)")
    assert edge == "NE(0500/0300)>NE(0499/0301)"
    assert topology.directed_edge_to_cells(edge) == (INTERIOR, "NE(0499/0301)")


def test_a_directed_edge_needs_edge_adjacency() -> None:
    with pytest.raises(DomainError, match="do not share an edge"):
        topology.cells_to_directed_edge(INTERIOR, "NE(0500/0301)")


def test_directed_edges_pair_by_position() -> None:
    edges = topology.cells_to_directed_edge(
        "NE(0500/0300,0501/0300)", "NE(0499/0300,0502/0300)"
    )
    assert edges == [
        "NE(0500/0300)>NE(0499/0300)",
        "NE(0501/0300)>NE(0502/0300)",
    ]
    assert topology.directed_edge_to_cells(",".join(edges)) == [
        ("NE(0500/0300)", "NE(0499/0300)"),
        ("NE(0501/0300)", "NE(0502/0300)"),
    ]


def test_mismatched_counts_are_refused() -> None:
    with pytest.raises(ValueError, match="same number of cells"):
        topology.cells_to_directed_edge("NE(0500/0300,0501/0300)", "NE(0499/0300)")


def test_a_malformed_edge_is_refused() -> None:
    with pytest.raises(InvalidIndexError, match="not a directed edge"):
        topology.directed_edge_to_cells("NE(0500/0300)")


def test_a_cell_has_one_edge_per_cardinal_neighbour() -> None:
    edges = topology.cell_to_edges(INTERIOR)
    assert len(edges) == 4
    for edge in edges:
        tail, head = topology.directed_edge_to_cells(edge)
        assert tail == INTERIOR
        assert share_an_edge(tail, head)


def test_the_meridian_triangle_has_four_edges_over_three_sides() -> None:
    """The paper says three, and is wrong on both halves.

    The triangle has three sides but four edge neighbours: its base is divided
    between the two quadrants. All four are returned; the poleward direction
    has none, and the equatorward direction has two, of which
    :func:`get_neighbor` returns the one sharing the longer boundary.
    """
    edges = topology.cell_to_edges("NE(0000/0300)")
    assert len(edges) == 4
    heads = {topology.directed_edge_to_cells(edge)[1] for edge in edges}
    assert heads == {
        "NE(0001/0299)",
        "NE(0001/0300)",
        "NW(0001/0299)",
        "NW(0001/0300)",
    }
    assert topology.get_neighbor("NE(0000/0300)", "N") is None


def test_edges_of_several_cells_are_positionally_aligned() -> None:
    grouped = topology.cell_to_edges("NE(0500/0300,0800/0400)")
    assert isinstance(grouped, list) and len(grouped) == 2
    assert all(len(group) == 4 for group in grouped)


# --------------------------------------------------------------------------
# The coupling to hierarchy, pinned rather than promoted
# --------------------------------------------------------------------------


def test_parent_cell_steps_west_when_the_prefix_names_no_cell() -> None:
    """The private helper this module depends on, and the case it exists for.

    ``topology`` uses ``hierarchy._parent_cell`` rather than ``get_parent``
    because a border trapezoid's children can be spelled under the column
    immediately east, and that column holds no cell of its own. The lexical
    prefix is then a string that names nothing, and composing neighbourhood
    through it would ascend into the void.

    Promoting the helper would add a public name and pull ``__init__.py`` into
    this phase. Pinning it costs one test and fails loudly if a later phase
    changes the signature or the westward step.
    """
    assert _parent_cell("NE(0500/0300(1))") == "NE(0500/0300)"

    row = ZONE_ROWS["CHUKOTKA"][0]
    trapezoid = f"NE({last_column('NE', row):04d}/{row:04d})"
    spilled = "NE(0934/0709(1))"
    assert spilled in itacart.get_children(trapezoid, flatten=True)
    assert itacart.get_parent(spilled) == "NE(0934/0709)"
    assert not itacart.is_valid_cell("NE(0934/0709)")
    assert _parent_cell(spilled) == trapezoid


# --------------------------------------------------------------------------
# The private helpers, at their contracts rather than through the public API
# --------------------------------------------------------------------------


def test_a_quadrant_counts_as_exceptional() -> None:
    """It has no lattice pair, so nothing can be derived from it.

    The public entry points refuse a quadrant on resolution before ever asking
    this, but the predicate has to answer for itself: it is called from
    ``_contacts``, which is reached by other routes.
    """
    assert topology._is_exceptional("NE")


def test_an_exceptional_ancestor_is_found_through_the_chain() -> None:
    """A parallelogram can descend from a triangle.

    Child ``2`` of a meridian triangle is an ordinary parallelogram, so the
    cell itself passes every shape test. Its neighbourhood is still not a
    grid, because the triangle above it in the chain does not tile as one, and
    only walking the prefixes finds that.
    """
    cell = "NE(0000/0300(2(C3)))"
    assert itacart.cell_shape(cell) == "parallelogram"
    assert not topology._is_exceptional(cell)
    assert topology._has_exceptional_ancestor(cell)

    neighbour = topology.get_neighbor(cell, "E")
    assert isinstance(neighbour, str)
    assert share_an_edge(cell, neighbour)


@pytest.mark.parametrize("quadrant", ["NE", "NW", "SE", "SW"])
def test_normalising_past_the_last_row_gives_the_polar_cell(quadrant: str) -> None:
    """Reached by ``deflect`` alone.

    ``get_neighbor`` never gets here: the cells of row 999 absorb the border,
    so they are in the exception set and are measured instead. The
    normalisation still has to be right, because it is the public function.
    """
    expected = f"{quadrant[0]}E(0000/1000)"
    assert topology.deflect(f"{quadrant}(0001/0999)", "N") == expected
    assert topology._eastern(quadrant) == f"{quadrant[0]}E"


def test_normalising_an_overshoot_the_mirror_row_cannot_absorb() -> None:
    """The seam mirror runs out of columns before the overshoot does.

    Rows near the pole hold one or two cells, so a target several columns past
    the last has no counterpart on the other side of the seam. Answering
    ``None`` is right; wrapping again would land back in the row it came from.
    """
    assert topology._normalize_target("NE", 3, 999, 0) == ("NW", 0, 999)
    assert topology._normalize_target("NE", 6, 999, 0) is None


def test_the_candidate_window_skips_rows_outside_the_domain() -> None:
    """The polar cell has no row above it to search."""
    edges = topology._edge_neighbors("NE(0000/1000)")
    assert edges == ("NE(0001/0999)", "NW(0001/0999)")
    for neighbour in edges:
        assert share_an_edge("NE(0000/1000)", neighbour)


def test_classify_answers_none_when_no_allowed_direction_fits() -> None:
    """Its contract, exercised where the caller never lets it happen.

    ``_neighbor_of_atom`` always passes the whole cardinal or diagonal set, so
    the miss cannot occur there. Pinning it here keeps the helper safe to
    reuse with a narrower set.
    """
    origin, corner = INTERIOR, "NE(0500/0301)"
    assert topology._classify(origin, corner, topology._DIAGONALS) == "NE"
    assert topology._classify(origin, corner, ("SW",)) is None


def test_shared_length_is_zero_for_cells_that_do_not_share_an_edge() -> None:
    """Both ways of not sharing one: far apart, and touching at a point.

    The second is the one worth pinning. A vertex contact is not empty — the
    intersection is a point — so a test written against emptiness would call
    it an edge.
    """
    assert topology._shared_length(INTERIOR, "NE(0600/0300)") == 0.0
    assert topology._shared_length(INTERIOR, "NE(0500/0301)") == 0.0
    assert share_only_a_vertex(INTERIOR, "NE(0500/0301)")


def test_the_neighbour_set_never_holds_a_cell_that_is_not_touched() -> None:
    """Regression: the index route is verified in both callers, or neither.

    ``NE(0886/0707)`` sits one row below a border trapezoid. Its poleward step
    spells a column past the end of the row above, so it deflects across the
    seam and lands on ``NW(0884/0708)`` — a real cell that this one does not
    touch. ``get_neighbor`` checked for contact and rejected it;
    ``_contacts`` did not, and so held a phantom while missing the trapezoid
    that is actually there.

    The perimeter is what catches it: with the phantom in and the trapezoid
    out, the neighbours covered 0.959 of a perimeter of 1.162.
    """
    cell = "NE(0886/0707)"
    trapezoid = "NE(0884/0708)"
    neighbours = topology._edge_neighbors(cell)

    assert trapezoid in neighbours
    assert "NW(0884/0708)" not in neighbours
    for neighbour in neighbours:
        assert share_an_edge(cell, neighbour)

    perimeter = _polygon(cell).exterior.length
    covered = sum(shared_edge_length(cell, other) for other in neighbours)
    assert abs(perimeter - covered) < 1e-9

    assert topology.get_neighbor(cell, "N") == trapezoid
    assert cell in topology.grid_disk(trapezoid, 1, "manhattan")
    assert trapezoid in topology.grid_disk(cell, 1, "manhattan")


@pytest.mark.slow
def test_unit_disk_membership_is_symmetric_across_the_border_families() -> None:
    """Criterion 8 where the opening package expected it to fall.

    Every family at once: both zone limits in both zones, an ordinary seam, the
    row under the pole, the polar cells and the meridian triangles. Both
    metrics, because the phantom above showed only in one of them.
    """
    border: list[str] = []
    for quadrant, row in (
        ("NE", 708),
        ("NE", 709),
        ("NW", 708),
        ("NW", 709),
        ("SE", 170),
        ("SE", 171),
        ("SW", 170),
        ("SW", 171),
        ("NE", 300),
        ("NE", 999),
        ("NE", 100),
        ("NE", 900),
    ):
        last = last_column(quadrant, row)
        border += [
            f"{quadrant}({column:04d}/{row:04d})"
            for column in range(max(0, last - 3), last + 1)
        ]
    border += ["NE(0000/1000)", "SE(0000/1000)", "NE(0000/0300)", "NE(0000/0000)"]

    for metric in ("manhattan", "chebyshev"):
        asymmetric = [
            (cell, other)
            for cell in border
            for other in topology.grid_disk(cell, 1, metric)
            if cell not in topology.grid_disk(other, 1, metric)
        ]
        assert asymmetric == []


def test_a_child_inherits_a_missing_parent_neighbour() -> None:
    """The descent has nothing to wrap into, so it answers ``None``.

    Near the equator the last columns of consecutive rows step by more than
    one, so a diagonal off ``NE(2002/0001)`` leaves the domain. A child on
    that side of the parent leaves with it; a child on the other side does
    not, which is what makes this a property of the step rather than of the
    cell.
    """
    assert topology.get_neighbor("NE(2002/0001)", "SE") is None
    assert topology.get_neighbor("NE(2002/0001(2))", "SE") is None
    assert topology.get_neighbor("NE(2002/0001(4))", "SE") == "NE(2003/0001(1))"


def test_classify_answers_none_for_cells_that_do_not_touch() -> None:
    """No shared boundary, so no side of the ring carries it."""
    assert topology._classify(INTERIOR, "NE(0600/0300)", topology._CARDINALS) is None


@pytest.mark.slow
def test_the_candidate_window_is_wider_than_the_zones_need() -> None:
    """The window constants are sized from measurement, not from argument.

    ``_geometric_candidates`` searches three bands per neighbouring row: around
    the cell's own column, at the prime meridian, and at the antemeridian seam.
    The width of the last two is a constant, chosen because it covered the two
    zones that exist. Nothing made it fail if a future zone reached further
    past the antemeridian than Chukotka does, and a window that is too narrow
    loses a neighbour in silence.

    So the sufficiency is measured. Both zones are enumerated at both latitude
    limits, in all four quadrants, and each cell's neighbour set is compared
    against the set a window four columns wider would find. Equal sets mean the
    current width is enough; a wider zone that needs more will make them
    differ.
    """
    cells: list[str] = []
    for first, last in ZONE_ROWS.values():
        for row in (first - 1, first, last, last + 1):
            for quadrant in ("NE", "NW", "SE", "SW"):
                edge = last_column(quadrant, row)
                cells += [
                    f"{quadrant}({column:04d}/{row:04d})"
                    for column in range(max(0, edge - 2), edge + 1)
                ]
    cells = [cell for cell in cells if itacart.is_valid_cell(cell)]
    assert len(cells) >= 90

    narrow = {cell: topology._contacts(cell) for cell in cells}

    window, seam = topology._WINDOW, topology._SEAM_REACH
    try:
        topology._WINDOW, topology._SEAM_REACH = window + 4, seam + 4
        topology._contacts.cache_clear()
        wide = {cell: topology._contacts(cell) for cell in cells}
    finally:
        topology._WINDOW, topology._SEAM_REACH = window, seam
        topology._contacts.cache_clear()

    missed = {
        cell: (narrow[cell], wide[cell]) for cell in cells if narrow[cell] != wide[cell]
    }
    assert missed == {}


@pytest.mark.slow
def test_the_poleward_step_composes_while_it_stays_inside_the_row_above() -> None:
    """The rule's global consequence, enumerated rather than sampled.

    ``north^k(X, Y) = (X - k, Y + k)`` holds for as long as the target column
    still fits the row above. A whole column is walked from the equator to
    prove where that stops, in both directions the walk can end.

    From column 500 the walk reaches the prime meridian first: five hundred
    steps of pure arithmetic, no deflection, ending on the triangle at column
    zero, which has no poleward neighbour. From column 1500 the walk meets the
    edge of the row before it meets the meridian, and from there on it rides
    the border and the arithmetic no longer describes it — but it does reach
    the pole, crossing quadrants twice on the way.
    """

    def parts(cell: str) -> tuple[str, int, int]:
        return cell[:2], int(cell[3:7]), int(cell[8:12])

    def walk(start: str) -> tuple[str, int, list[int], list[int]]:
        cell, departures, crossings = start, [], []
        for step in range(1000):
            following = topology.get_neighbor(cell, "N")
            if following is None:
                return cell, step, departures, crossings
            quadrant, column, row = parts(cell)
            if parts(following) != (quadrant, column - 1, row + 1):
                departures.append(step + 1)
            if following[:2] != cell[:2]:
                crossings.append(step + 1)
            cell = following
        return cell, 1000, departures, crossings

    ended, steps, departures, crossings = walk("NE(0500/0000)")
    assert (ended, steps) == ("NE(0000/0500)", 500)
    assert departures == []
    assert crossings == []
    assert topology.get_neighbor(ended, "N") is None

    ended, steps, departures, crossings = walk("NE(1500/0000)")
    assert (ended, steps) == ("NE(0000/1000)", 1000)
    assert departures[0] == 778
    assert crossings == [800, 1000]

    first = "NE(0723/0777)"
    quadrant, column, row = parts(first)
    assert last_column(quadrant, row + 1) < column - 1


@pytest.mark.slow
def test_a_directed_step_is_not_symmetric_at_the_boundary() -> None:
    """Declared, because composing steps is a reasonable thing to want to do.

    Adjacency is symmetric and the disks are symmetric. A *named direction* is
    not, and two separate mechanisms break it.

    The first is mirroring. Direction labels are read in the origin's own
    quadrant, so ``E`` means away from the prime meridian on both sides of it.
    Stepping west off the meridian triangle lands in the western quadrant, and
    stepping east from there goes further west still, not back. The opposite
    label is simply not the inverse of a step that changes quadrant.

    The second is the tie-break. Where a side carries several cells the step
    names the one sharing the longest boundary; that cell is ordinary, has
    four directions, and spends them where its own lattice says.

    Every pair below is genuinely adjacent, which is what makes this a
    statement about directions and not about the tessellation.
    """
    opposite = {
        "N": "S",
        "S": "N",
        "E": "W",
        "W": "E",
        "NE": "SW",
        "SW": "NE",
        "NW": "SE",
        "SE": "NW",
    }

    cells: list[str] = []
    for quadrant, row in (
        ("NE", 708),
        ("NE", 709),
        ("NW", 708),
        ("NW", 709),
        ("SE", 170),
        ("SE", 171),
        ("NE", 300),
        ("NE", 999),
    ):
        edge = last_column(quadrant, row)
        cells += [
            f"{quadrant}({column:04d}/{row:04d})"
            for column in range(max(0, edge - 3), edge + 1)
        ]
        cells.append(
            f"{quadrant}(0000/{row:04d})"
            if quadrant[1] == "E"
            else f"{quadrant}(0001/{row:04d})"
        )
    cells = [cell for cell in cells if itacart.is_valid_cell(cell)]

    crossing, inside = 0, 0
    for cell in cells:
        for direction, back in opposite.items():
            neighbour = topology.get_neighbor(cell, direction)
            if neighbour is None or topology.get_neighbor(neighbour, back) == cell:
                continue
            if neighbour[:2] != cell[:2]:
                crossing += 1
            else:
                inside += 1
            if direction in ("N", "S", "E", "W"):
                assert topology.are_neighbor_cells(cell, neighbour)

    assert crossing > 0, "mirroring should break the inverse across a quadrant"
    assert inside > 0, "the tie-break should break it inside one"

    # The two named cases, so the mechanisms are pinned and not just counted.
    assert topology.get_neighbor("NE(0000/0708)", "W") == "NW(0001/0708)"
    assert topology.get_neighbor("NW(0001/0708)", "E") == "NW(0002/0708)"
    assert topology.are_neighbor_cells("NE(0000/0708)", "NW(0001/0708)")

    assert topology.get_neighbor("NE(0932/0709)", "N") == "NE(0930/0710)"
    assert topology.get_neighbor("NE(0930/0710)", "S") == "NE(0931/0709)"
    assert topology.are_neighbor_cells("NE(0932/0709)", "NE(0930/0710)")


def test_the_index_route_stays_lexical_and_the_exception_set_does_not() -> None:
    """Criterion 9, with its scope stated instead of implied.

    The derivation consults nothing, anywhere: ``_lexical_target`` answers with
    every geometric function and both existence predicates replaced by raisers.

    Resolving a neighbour is a weaker claim, and the weakness is exactly the
    exception set. Outside it the step is arithmetic plus a validity question,
    so it survives with the rings unavailable. Inside it the neighbour is
    measured, so it does not, and saying otherwise would make the criterion
    vacuous rather than satisfied.
    """

    def refuse(*args: object, **kwargs: object) -> object:
        raise AssertionError("consulted geometry")

    lexical = ("NE(0500/0300)", "NE(0001/0300)", "NE(0500/0300(1(A1)))")
    measured = ("NE(0000/0300)", "NE(1784/0300)", "NE(0000/1000)", "NW(0883/0708)")

    saved = itacart.cells.cell_to_boundary
    try:
        topology._contacts.cache_clear()
        itacart.cells.cell_to_boundary = refuse  # type: ignore[assignment]
        topology.cell_to_boundary = refuse  # type: ignore[assignment]
        for cell in lexical:
            assert isinstance(topology.get_neighbor(cell, "E"), str)
        for cell in measured:
            with pytest.raises(AssertionError, match="consulted geometry"):
                topology.get_neighbor(cell, "E")
    finally:
        itacart.cells.cell_to_boundary = saved  # type: ignore[assignment]
        topology.cell_to_boundary = saved  # type: ignore[assignment]
        topology._contacts.cache_clear()


def test_a_step_inside_one_frame_is_never_checked_against_geometry() -> None:
    """The interior path must not pay for the boundary's problem.

    Only a frame change can make the arithmetic non-invertible, and every
    boundary the normalisation handles changes quadrant, apart from the step
    onto a polar cell. Testing instead whether the step left its resolution-1
    cell would be true of every step at resolution 1 and would put a ring
    computation on the interior path.
    """
    assert not topology._changes_frame(INTERIOR, "NE(0499/0301)")
    assert not topology._changes_frame("NE(0500/0300(1(A1)))", "NE(0499/0300(2(A5)))")
    assert topology._changes_frame("NE(0000/0300)", "NW(0001/0300)")
    assert topology._changes_frame("NE(0001/0999)", "NE(0000/1000)")
