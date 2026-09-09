"""The parent-to-child contract on cells that absorb the domain border.

The contract is one sentence: the effective rings of the children cover
the effective ring of the parent, and their interiors do not overlap.
``boundary.plane_ring`` is the authority for "effective" -- it answers
after the pole and after absorption -- and every claim here is measured
against it rather than against the nominal lattice.

Why the nominal lattice is not enough is itself measured, in
:func:`test_the_nominal_footprint_does_not_cover_an_absorbing_parent`.
On the polar cap it covers 78.77% of the parent, and the shortfall does
not shrink with depth, because the absorbed surface is annexed at one
level and dropped by the next descent. Anything that enumerates children
by descending the footprint loses that surface permanently.

Three parents carry the regression. ``A5`` and ``E1`` are the two
trapezoids of the cap's first sub-row, mirror images of one another.
``B2`` is the polar triangle, whose effective ring is 2.93 times its
nominal one and which is therefore the hardest case in the grid.
"""

from __future__ import annotations

import pytest
from shapely.geometry import Polygon

import itacart
from itacart import boundary
from itacart import hierarchy as hy
from itacart.constants import CELL_SIZE_M, refinement_alphabet
from itacart.resolutions import refinement_ratio

POLAR_CAP = "NE(0000/1000)"
SOUTH_CAP = "SE(0000/1000)"
CAP_SUBDIVISION = "NE(0000/1000(1))"
EAST_TRAPEZOID = "NE(0000/1000(1(A5)))"
WEST_TRAPEZOID = "NE(0000/1000(1(E1)))"
POLAR_TRIANGLE = "NE(0000/1000(1(B2)))"

ABSORBING_PARENTS = [
    POLAR_CAP,
    CAP_SUBDIVISION,
    EAST_TRAPEZOID,
    WEST_TRAPEZOID,
    POLAR_TRIANGLE,
    SOUTH_CAP,
    "NE(1414/0500)",
    "NE(0311/0900)",
    "SW(1845/0237)",
    "NE(0817/0747)",
]

# The share of a parent's area that may stay uncovered, and the share of
# a cell side squared that sibling interiors may share. Both match the
# tolerances the enumeration itself proves against, so a test that passes
# here is not passing on a slacker rule than the code applies.
COVERAGE_TOLERANCE = 1e-6
OVERLAP_TOLERANCE = 1e-6


def _body(cell: str) -> Polygon:
    """The effective ring of ``cell`` as a polygon."""
    _, ring = boundary.plane_ring(cell)
    assert ring, f"{cell} has no effective ring"
    body = Polygon(ring)
    assert body.is_valid, f"{cell} has a self-intersecting effective ring"
    return body


def _outside_share(child: Polygon, parent: Polygon) -> float:
    """How much of ``child`` lies outside ``parent``, as a share of itself.

    By areas rather than by ``difference``. Sibling and parent rings share
    edges exactly, and the difference of such rings is where the geometry
    engine reports a side location conflict; the areas ask the same
    question without building the shape.
    """
    return (child.area - child.intersection(parent).area) / child.area


def _every_absorbing_cell_of(quadrant: str, step: int) -> list[str]:
    """The last-column cell of every ``step``-th row of ``quadrant``."""
    side = CELL_SIZE_M[1]
    assert side is not None
    family = []
    for row in range(0, 1000, step):
        column = itacart.last_lattice_column(quadrant, row, side)
        if column < 0:
            continue
        cell = f"{quadrant}({column:04d}/{row:04d})"
        if boundary.is_valid_cell(cell):
            family.append(cell)
    return family


class TestTheContract:
    """The five properties that hold, asserted over every named parent."""

    @pytest.mark.parametrize("parent", ABSORBING_PARENTS)
    def test_every_child_is_a_cell(self, parent: str) -> None:
        children = hy._children_of(parent)
        assert children
        for child in children:
            assert boundary.is_valid_cell(child), child

    @pytest.mark.parametrize("parent", ABSORBING_PARENTS)
    def test_sibling_interiors_do_not_overlap(self, parent: str) -> None:
        rings = [_body(child) for child in hy._children_of(parent)]
        shared = sum(
            first.intersection(second).area
            for index, first in enumerate(rings)
            for second in rings[index + 1 :]
        )
        assert shared <= OVERLAP_TOLERANCE * _body(parent).area

    @pytest.mark.parametrize("parent", ABSORBING_PARENTS)
    def test_the_children_cover_the_whole_effective_ring(self, parent: str) -> None:
        """The half of the contract the nominal descent cannot meet."""
        body = _body(parent)
        covered = sum(
            _body(child).intersection(body).area for child in hy._children_of(parent)
        )
        assert body.area - covered <= COVERAGE_TOLERANCE * body.area

    @pytest.mark.parametrize("parent", ABSORBING_PARENTS)
    def test_every_child_resolves_back_to_this_parent(self, parent: str) -> None:
        """The relation is a function both ways, spelling notwithstanding.

        A child of an absorbing cell can be spelled under a neighbour of
        its parent, by a column of the resolution-1 lattice or by a
        sibling code at an intermediate level. Neither spelling changes
        who the parent is, and this is what says so.
        """
        for child in hy._children_of(parent):
            assert hy._parent_cell(child) == parent, child

    def test_an_ordinary_cell_still_refines_into_exactly_the_alphabet(self) -> None:
        """Property seven: nothing outside the border family moved."""
        for cell in (
            "NE(0100/0100)",
            "SW(0500/0250)",
            "NE(0000/0500)",
            "SE(0007/0007)",
        ):
            for level in (2, 3):
                parent = cell if level == 2 else hy._children_of(cell)[0]
                if boundary.absorbs_border(parent):
                    continue
                children = hy._children_of(parent)
                assert len(children) == refinement_ratio(level), parent
                assert children == [
                    hy._descend(parent, code) for code in refinement_alphabet(level)
                ], parent


class TestThePolarCap:
    """The cap, where the nominal footprint fails hardest."""

    def test_the_nominal_footprint_does_not_cover_an_absorbing_parent(self) -> None:
        """Why children cannot be the descendants of the nominal lattice.

        Measured on the cap's first subdivision: the nominal rings of its
        own children cover 78.77% of it. The remaining fifth is surface
        the parent holds by absorption, and a descent that starts from
        the footprint never proposes a cell for it. The shortfall is a
        property of the footprint, not of depth, so it does not shrink
        one level down.
        """
        body = _body(CAP_SUBDIVISION)
        nominal = [
            Polygon(boundary._nominal_ring(child)[1])
            for child in hy._children_of(CAP_SUBDIVISION)
        ]
        covered = sum(ring.intersection(body).area for ring in nominal if ring.is_valid)
        assert covered / body.area < 0.79
        assert covered / body.area > 0.78

    def test_the_polar_triangle_is_the_hardest_parent_and_still_closes(self) -> None:
        """``B2``: effective ring 2.93 times the nominal one, tiled exactly.

        Its children are spelled under three different intermediate
        codes, two of which name no cell of their own. Four further
        spellings under that parent carry folded rings and are not
        children: the tiling closes without them, which is what says
        that a folded ring is an artefact of absorption rather than a
        cell.
        """
        body = _body(POLAR_TRIANGLE)
        nominal = Polygon(boundary._nominal_ring(POLAR_TRIANGLE)[1])
        assert body.area / nominal.area > 2.9

        children = hy._children_of(POLAR_TRIANGLE)
        assert len(children) == 8
        prefixes = {child[: child.rindex("(")] for child in children}
        assert len(prefixes) == 3

        covered = sum(_body(child).intersection(body).area for child in children)
        assert body.area - covered <= COVERAGE_TOLERANCE * body.area

    def test_the_cap_enumerates_symmetrically_east_and_west(self) -> None:
        """Property six, asserted on extents rather than on names.

        The two sides cannot be checked by swapping a quadrant letter or
        transposing a code, because a child's spelling can change stem
        during refinement. What is symmetric is the geometry, so that is
        what is compared: the outer reach of the western children is the
        outer reach of the eastern ones, mirrored.
        """
        for parent in (CAP_SUBDIVISION, POLAR_TRIANGLE):
            west, east = [], []
            for child in hy._children_of(parent):
                xs = [x for x, _ in boundary.plane_ring(child)[1]]
                if max(xs) <= 1e-9:
                    west.append(round(-min(xs), 6))
                elif min(xs) >= -1e-9:
                    east.append(round(max(xs), 6))
            assert west, parent
            assert sorted(west) == sorted(east), parent

    def test_the_two_trapezoids_of_the_first_sub_row_are_mirrors(self) -> None:
        """``A5`` and ``E1``, the regression pair, closing identically."""
        counts, areas = [], []
        for parent in (EAST_TRAPEZOID, WEST_TRAPEZOID):
            body = _body(parent)
            children = hy._children_of(parent)
            covered = sum(_body(child).intersection(body).area for child in children)
            assert body.area - covered <= COVERAGE_TOLERANCE * body.area
            counts.append(len(children))
            areas.append(round(body.area, 3))
        assert counts[0] == counts[1]
        assert areas[0] == areas[1]


class TestTheFamilyAsAWhole:
    """Totality instruments: the family, not an example from it."""

    @pytest.mark.slow
    def test_the_contract_holds_over_every_seventh_row_of_every_quadrant(self) -> None:
        """Coverage and disjointness across the whole antemeridian family.

        Every seventh row rather than every row, because the property is
        being asserted over a family and the cost is polygonal. The
        stride is coprime with nothing in the grid's periods, so it does
        not sample one phase of the border.
        """
        checked = 0
        for quadrant in ("NE", "NW", "SE", "SW"):
            for parent in _every_absorbing_cell_of(quadrant, 7):
                body = _body(parent)
                rings = [_body(child) for child in hy._children_of(parent)]
                covered = sum(ring.intersection(body).area for ring in rings)
                shared = sum(
                    first.intersection(second).area
                    for index, first in enumerate(rings)
                    for second in rings[index + 1 :]
                )
                assert body.area - covered <= COVERAGE_TOLERANCE * body.area, parent
                assert shared <= OVERLAP_TOLERANCE * body.area, parent
                checked += 1
        assert checked > 500

    @pytest.mark.slow
    def test_the_child_count_of_the_family_runs_from_two_to_seven(self) -> None:
        """The bound, re-measured after the enumeration was corrected.

        Seven is not an instance. Over the four quadrants it occurs in 45
        cells of 4 000, and the old bound of six was measured against an
        enumeration that reached one column east and no further -- which
        missed children spelled two columns out.
        """
        counts: dict[int, int] = {}
        for quadrant in ("NE", "NW", "SE", "SW"):
            for parent in _every_absorbing_cell_of(quadrant, 1):
                count = len(hy._children_of(parent))
                counts[count] = counts.get(count, 0) + 1
        assert min(counts) == 2
        assert max(counts) == 7
        assert counts[7] == 45
        assert sum(counts.values()) == 4000


class TestContainmentIsNotYetAProperty:
    """Property two, which does not hold, measured rather than assumed.

    A child of a border-absorbing cell reaches past its parent. Not by an
    epsilon, and not by accident: re-measured over every seventh row of
    every quadrant, all 572 parents are affected and exactly two children
    of each leave, which is the two trapezoids. The share of a child's own
    area that falls outside runs to a median of 0.0381%, a mean of
    0.0569%, a ninetieth percentile of 0.1041% and a worst case of 2.4110%
    at ``SE(1930/0196)``.

    The cause is a chord, not a staircase. ``_x_border`` does not depend
    on the lattice row at all -- measuring that is what killed the first
    explanation. The border is a smooth convex curve on the plane, and an
    absorbing cell replaces its outer side by the chord of that curve
    across its own height. A parent spans one cell height and takes one
    chord; its children span half that each and take two, and a chord of a
    convex curve over a shorter span lies further out.

    The model is exact, which is what
    :func:`test_the_excess_is_the_chord_model_to_four_places` asserts: the
    excess area equals the sagitta times half the parent's height, to four
    decimal places, over parents whose sagitta runs from 0.95 to 6.12
    metres.

    None of this is the enumeration's doing. The child rings are what
    ``plane_ring`` returns and the child set is the one the contract
    proves. This is pinned, not asserted away: when the absorbing side
    stops being a per-cell chord, these numbers go to zero and this test
    is the one that should fail.
    """

    def test_a_child_reaches_past_its_parent_on_the_lateral_border(self) -> None:
        for parent in ("NE(1414/0500)", "SW(1845/0237)"):
            body = _body(parent)
            shares = [
                _outside_share(_body(child), body) for child in hy._children_of(parent)
            ]
            assert max(shares) > 1e-6, parent

    @pytest.mark.slow
    def test_it_happens_across_the_family_and_stays_under_three_percent(self) -> None:
        worst = 0.0
        parents = 0
        reaching = 0
        for quadrant in ("NE", "NW", "SE", "SW"):
            for parent in _every_absorbing_cell_of(quadrant, 7):
                body = _body(parent)
                shares = [
                    _outside_share(_body(child), body)
                    for child in hy._children_of(parent)
                ]
                parents += 1
                if max(shares) > 1e-6:
                    reaching += 1
                worst = max(worst, max(shares))
        assert reaching == parents
        assert worst < 0.03

    def test_the_excess_is_the_chord_model_to_four_places(self) -> None:
        """The closed form, which is what a fix has to remove.

        One chord over a full cell height against two over half heights
        encloses two triangles, each of base half the height and height
        the sagitta at the midpoint. Their total is the sagitta times half
        the parent's height, and that is the measured excess exactly.
        """
        for quadrant, row in (("SE", 196), ("NE", 500), ("NE", 100), ("SW", 237)):
            side = CELL_SIZE_M[1]
            assert side is not None
            column = itacart.last_lattice_column(quadrant, row, side)
            parent = f"{quadrant}({column:04d}/{row:04d})"
            if not boundary.is_valid_cell(parent):
                continue
            whole = _body(parent)
            measured = sum(
                _body(child).area - _body(child).intersection(whole).area
                for child in hy._children_of(parent)
            )
            ys = [y for _, y in boundary.plane_ring(parent)[1]]
            low, high = min(ys), max(ys)
            edges = boundary._x_border
            sagitta = abs(
                edges(quadrant, row, (low + high) / 2.0)
                - (edges(quadrant, row, low) + edges(quadrant, row, high)) / 2.0
            )
            predicted = sagitta * (high - low) / 2.0
            assert sagitta > 0.9, parent
            assert measured == pytest.approx(predicted, rel=1e-4), parent

    def test_the_measured_shares_stay_where_they_were_re_measured(self) -> None:
        """The four numbers the docstring quotes, over a smaller stride."""
        shares = []
        for quadrant in ("NE", "SE"):
            for parent in _every_absorbing_cell_of(quadrant, 29):
                whole = _body(parent)
                shares.extend(
                    _outside_share(_body(child), whole)
                    for child in hy._children_of(parent)
                    if _outside_share(_body(child), whole) > 1e-9
                )
        assert shares
        shares.sort()
        median = shares[len(shares) // 2]
        assert 0.0003 < median < 0.0005
        assert max(shares) < 0.03

    def test_the_polar_cap_does_not_have_this_defect(self) -> None:
        """The cap is clipped by the pole, not by the staircase border."""
        for parent in (CAP_SUBDIVISION, POLAR_TRIANGLE, EAST_TRAPEZOID, WEST_TRAPEZOID):
            body = _body(parent)
            for child in hy._children_of(parent):
                assert _outside_share(_body(child), body) < 1e-6, child


class TestFoldedSpellingsAreNotCells:
    """The regression for what a folded ring is, and is not.

    Absorbing a border steeper than the lattice shear can fold the
    quadrilateral the absorption builds. The shoelace area of such a ring
    is the difference of its two lobes rather than nothing, so an area
    test accepts it, and the existence predicate used to.

    Nothing repairs one. The pair below reports 104 724.4 square metres
    from rings that cover 322 621.8 once rebuilt, and the surface they
    seem to claim already belongs to the trapezoid that absorbed it.
    """

    FOLDED = ("NE(0000/1000(3(A1)))", "NE(0000/1000(2(A1)))")
    ABSORBERS = (WEST_TRAPEZOID, EAST_TRAPEZOID)

    @pytest.mark.parametrize("cell", FOLDED)
    def test_a_folded_ring_does_not_name_a_cell(self, cell: str) -> None:
        assert not boundary.is_valid_cell(cell)

    @pytest.mark.parametrize("cell", FOLDED)
    def test_the_area_alone_would_have_accepted_it(self, cell: str) -> None:
        """Which is why the predicate cannot rest on area."""
        ring = boundary.plane_ring(cell)[1]
        assert boundary.ring_area(ring) > 100_000.0
        assert not Polygon(ring).is_valid

    @pytest.mark.parametrize("absorber", ABSORBERS)
    def test_the_absorbing_twin_is_a_cell_and_holds_that_surface(
        self, absorber: str
    ) -> None:
        assert boundary.is_valid_cell(absorber)
        effective = _body(absorber)
        nominal = Polygon(boundary._nominal_ring(absorber)[1])
        assert effective.is_valid
        assert round(effective.area - nominal.area, 1) == 104_724.4

    def test_no_child_of_an_absorbing_parent_carries_a_folded_ring(self) -> None:
        for parent in ABSORBING_PARENTS:
            for child in hy._children_of(parent):
                assert _body(child).is_valid, child

    def test_a_folded_spelling_is_never_proposed_as_a_child(self) -> None:
        """``B2`` closes exactly, and none of the four folds is in the set."""
        folded = {
            "NE(0000/1000(1(D2(1))))",
            "NE(0000/1000(1(D2(2))))",
            "NE(0000/1000(1(B4(1))))",
            "NE(0000/1000(1(B4(2))))",
        }
        children = set(hy._children_of(POLAR_TRIANGLE))
        assert not children & folded
        for cell in folded:
            assert not boundary.is_valid_cell(cell)


class TestTheSimplicityClauseIsDifferential:
    """What the new clause in the existence predicate changed, enumerated.

    Adding "and the ring closes without crossing itself" to
    :func:`itacart.boundary.is_valid_cell` is a widening of a refusal, and
    a widening has to be measured against the population it acts on
    rather than against the examples that motivated it. Two things are
    asserted over that population: the change is one-way, and every
    spelling it refuses is refused for the stated reason.

    Measured over 32 656 spellings -- the full alphabet of both polar caps
    to resolution 5, and the lateral family every thirty-seventh row to
    resolution 3 -- 196 spellings changed verdict, none in the other
    direction, and none of them carries a simple ring.
    """

    @staticmethod
    def _population() -> list[str]:
        side = CELL_SIZE_M[1]
        assert side is not None
        population: list[str] = []
        for quadrant in ("NE", "SE"):
            level_cells = [f"{quadrant}(0000/1000)"]
            for level in (2, 3, 4, 5):
                level_cells = [
                    hy._descend(cell, code)
                    for cell in level_cells
                    for code in refinement_alphabet(level)
                ]
                population.extend(level_cells)
        for quadrant in ("NE", "NW", "SE", "SW"):
            for row in range(0, 1000, 37):
                column = itacart.last_lattice_column(quadrant, row, side)
                if column < 0:
                    continue
                level_cells = [f"{quadrant}({column:04d}/{row:04d})"]
                for level in (2, 3):
                    level_cells = [
                        hy._descend(cell, code)
                        for cell in level_cells
                        for code in refinement_alphabet(level)
                    ]
                    population.extend(level_cells)
        return population

    @pytest.mark.slow
    def test_the_clause_only_ever_refuses_and_only_ever_for_a_fold(self) -> None:
        population = self._population()
        assert len(population) > 30_000

        refused = []
        for cell in population:
            ring = boundary._safe_ring(cell)[1]
            had_area = boundary.ring_area(ring) > boundary._AREA_EPSILON_M2
            exists = bool(boundary.is_valid_cell(cell))
            # One-way: nothing the old rule refused is accepted now.
            assert not (exists and not had_area), cell
            if had_area and not exists:
                refused.append(cell)

        assert refused
        for cell in refused:
            assert not boundary._is_simple_ring(
                boundary._safe_ring(cell)[1]
            ), f"{cell} was refused although its ring is simple"

    def test_an_ordinary_ring_and_a_degenerate_one_are_told_apart(self) -> None:
        """The clause itself, on rings written here rather than found."""
        square = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        bowtie = [(0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0)]
        sliver = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]
        assert boundary._is_simple_ring(square)
        assert not boundary._is_simple_ring(bowtie)
        # Collinear and enclosing nothing is degenerate, not folded: the
        # area test is what refuses it, and this clause must not.
        assert boundary._is_simple_ring(sliver)
        assert not boundary._is_simple_ring(square[:2])
