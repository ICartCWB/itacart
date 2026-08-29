"""Tests for :mod:`itacart.index`.

The suite is organised by acceptance criterion of F2, then by surface.
The tests around the even and odd alphabets are permanent regression
guards: the prototype mapped even resolutions to ``A1``-``E5`` and odd
ones to ``1``-``4``, which is the inversion the paper's section 3.1
forbids.
"""

from __future__ import annotations

from itertools import combinations
from typing import Iterator

import pytest

from itacart import index as ix
from itacart.constants import (
    MAX_RESOLUTION,
    QUADRANTS,
    QUATERNARY_CODES,
    QUINARY_CODES,
    RES1_DIGITS,
    RES1_MAX_INDEX,
    refinement_alphabet,
)
from itacart.exceptions import (
    InvalidIndexError,
    InvalidQuadrantError,
    InvalidRefinementCodeError,
    ITACaRTError,
    NonAtomicIndexError,
)

# --------------------------------------------------------------------------
# Corpus
# --------------------------------------------------------------------------

ALL_QUATERNARY = ",".join(QUATERNARY_CODES)
ALL_QUINARY = ",".join(QUINARY_CODES)


def _every_resolution() -> Iterator[int]:
    """Resolutions addressed by a refinement code: 2 to 13."""
    return iter(range(2, MAX_RESOLUTION + 1))


def _chain_to(resolution: int, *, tail: str = "") -> str:
    """Build a single-path index descending to ``resolution``.

    The refinement code at each level is taken from
    :func:`refinement_alphabet`, so the helper cannot itself encode the
    even/odd rule the tests exist to check.
    """
    codes = [refinement_alphabet(level)[0] for level in range(2, resolution + 1)]
    if tail:
        codes.append(tail)
    body = "(".join(["0001/0002", *codes])
    return f"NE({body}{')' * (len(codes) + 1)}"


def _graft(index: str, subtree: str) -> str:
    """Hang ``subtree`` under the deepest node of a single-path ``index``."""
    closers = len(index) - len(index.rstrip(")"))
    return f"{index[:-closers]}({subtree}){')' * closers}"


#: Indices spanning every shape the grammar admits. Used by the property
#: tests, which must hold over all of them and not only over Figure 7.
CORPUS: tuple[str, ...] = (
    "NE",
    "SE(1400/0374)",
    "SE(1400/0374(3))",
    "SE(1400/0374(3(C2(3))))",
    "SE(1400/0374(1,2))",
    "SE(1400/0374(3(C2,C3)))",
    "SE(1400/0374,1401/0374)",
    _chain_to(MAX_RESOLUTION),
    "NW(0625/0451(1(E1(3(B2(4(A2,B2)))))))",
    "SW(0001/0002(4)),NE(0003/0004(2))",
)


# --------------------------------------------------------------------------
# Criterion 1 - Figure 7 round-trip
# --------------------------------------------------------------------------


def test_figure_7_roundtrips_through_decompose_and_compose(
    central_park_index: str,
) -> None:
    """parse -> decompose -> compose returns the published string itself."""
    assert ix.compose(ix.decompose(central_park_index)) == central_park_index


def test_figure_7_decomposes_to_114_atomic_cells(central_park_index: str) -> None:
    cells = ix.decompose(central_park_index)
    assert len(cells) == 114
    assert len(set(cells)) == 114


def test_figure_7_cells_are_all_atomic(central_park_index: str) -> None:
    assert all(ix.is_atomic(cell) for cell in ix.decompose(central_park_index))


def test_figure_7_count_cells_agrees_with_decompose(central_park_index: str) -> None:
    assert ix.count_cells(central_park_index) == len(ix.decompose(central_park_index))


def test_figure_7_is_already_in_canonical_form(central_park_index: str) -> None:
    """The published index needs no rewriting; canonicalisation is a no-op."""
    assert ix.normalize(central_park_index) == central_park_index


def test_figure_7_spans_resolutions_6_and_7(central_park_index: str) -> None:
    """The caption says resolutions 6 and 7; the string must agree."""
    depths = {
        len(ix.split_components(cell)) - 1 for cell in ix.decompose(central_park_index)
    }
    assert depths == {6, 7}


@pytest.mark.parametrize("index", CORPUS)
def test_compose_inverts_decompose_over_the_corpus(index: str) -> None:
    assert ix.compose(ix.decompose(index)) == index


def test_paper_example_roundtrips(paper_example_index: str) -> None:
    assert ix.compose(ix.decompose(paper_example_index)) == paper_example_index
    assert ix.split_components(paper_example_index) == [
        "SE",
        "1400/0374",
        "3",
        "C2",
        "3",
    ]


def test_sydney_cell_is_a_resolution_7_single_path(sydney_cell: str) -> None:
    """The fixture docstring claims resolution 9; the string is 7."""
    assert ix.is_atomic(sydney_cell)
    assert len(ix.split_components(sydney_cell)) - 1 == 7


# --------------------------------------------------------------------------
# Criterion 2 - completeness collapse
# --------------------------------------------------------------------------


def test_normalize_collapses_a_complete_quaternary_set() -> None:
    assert ix.normalize(f"SE(1400/0374({ALL_QUATERNARY}))") == "SE(1400/0374)"


def test_normalize_collapses_a_complete_quaternary_set_under_a_quinary_parent() -> None:
    """The realisable form of criterion 2's ``4(1,2,3,4)`` example.

    A quaternary child list can only hang off a resolution-1 node or off a
    quinary one, never off a quaternary one — see the test below.
    """
    assert ix.normalize(f"SE(1400/0374(3(C2({ALL_QUATERNARY}))))") == (
        "SE(1400/0374(3(C2)))"
    )


@pytest.mark.parametrize("resolution", [2, 4, 6, 8, 10, 12])
def test_criterion_2_literal_example_is_unrealisable(resolution: int) -> None:
    """``4(1,2,3,4)`` is not a valid fragment at any level.

    ``4`` is a quaternary code, so it lives at an even resolution, so its
    children are odd and therefore quinary. The alternating refinement
    makes a quaternary node with quaternary children impossible. The
    briefing's criterion 2 quotes this string; the property it stands for
    is closed by the two collapse tests above, not by this one.
    """
    node_four = _chain_to(resolution, tail="")
    assert ix.is_valid_index(node_four)
    with pytest.raises(InvalidRefinementCodeError):
        ix.parse(_graft(node_four, ALL_QUATERNARY))


def test_normalize_collapses_a_complete_quinary_set() -> None:
    index = f"SE(1400/0374(3({ALL_QUINARY})))"
    assert ix.normalize(index) == "SE(1400/0374(3))"


def test_normalize_collapses_recursively_through_two_levels() -> None:
    """A child that collapses to a whole cell lets its parent collapse too."""
    inner = ",".join(f"{code}({ALL_QUINARY})" for code in QUATERNARY_CODES)
    assert ix.normalize(f"SE(1400/0374({inner}))") == "SE(1400/0374)"


def test_normalize_does_not_collapse_an_incomplete_set() -> None:
    index = "SE(1400/0374(1,2,3))"
    assert ix.normalize(index) == index


def test_normalize_does_not_collapse_a_complete_set_with_a_partial_child() -> None:
    """All four codes present, but one is subdivided, so coverage is partial."""
    index = "SE(1400/0374(1,2,3,4(A1)))"
    assert ix.normalize(index) == index


def test_normalize_does_not_collapse_quadrants_into_the_globe() -> None:
    """Resolution 0 has no code above it, so a full set of quadrants stays."""
    index = ",".join(f"{quadrant}(0001/0001)" for quadrant in QUADRANTS)
    assert ix.normalize(index) == index


def test_normalize_merges_repeated_siblings() -> None:
    assert ix.normalize("SE(1400/0374(1),1400/0374(2))") == "SE(1400/0374(1,2))"


def test_normalize_lets_a_whole_cell_absorb_its_own_descendants() -> None:
    assert ix.normalize("SE(1400/0374(4,4(A1)))") == "SE(1400/0374(4))"


def test_normalize_pads_resolution_1_components() -> None:
    assert ix.normalize("SE(1/2)") == f"SE({'1'.rjust(RES1_DIGITS, '0')}/0002)"


def test_normalize_orders_siblings_by_alphabet_position() -> None:
    assert ix.normalize("SE(1400/0374(3(E5,A1,C3)))") == "SE(1400/0374(3(A1,C3,E5)))"


def test_normalize_orders_base_cells_numerically_not_lexically() -> None:
    assert ix.normalize("SE(0010/0001,0002/0001)") == "SE(0002/0001,0010/0001)"


def test_normalize_orders_quadrants_by_the_table_order() -> None:
    assert ix.normalize("SW(0001/0001),NE(0001/0001)") == (
        "NE(0001/0001),SW(0001/0001)"
    )


# --------------------------------------------------------------------------
# Criterion 3 - idempotence
# --------------------------------------------------------------------------


@pytest.mark.parametrize("index", CORPUS)
def test_normalize_is_idempotent(index: str) -> None:
    once = ix.normalize(index)
    assert ix.normalize(once) == once


@pytest.mark.parametrize(
    "index",
    [
        f"SE(1400/0374({ALL_QUATERNARY}))",
        "SE(1400/0374(1),1400/0374(2))",
        "SE(1/2)",
        "SE(1400/0374(3(E5,A1,C3)))",
        "SE(1400/0374(4,4(A1)))",
    ],
)
def test_normalize_is_idempotent_where_it_actually_rewrites(index: str) -> None:
    once = ix.normalize(index)
    assert once != index
    assert ix.normalize(once) == once


def test_normalize_is_idempotent_on_figure_7(central_park_index: str) -> None:
    once = ix.normalize(central_park_index)
    assert ix.normalize(once) == once


# --------------------------------------------------------------------------
# Criterion 4 - OGC Req 13, one canonical form per region
# --------------------------------------------------------------------------

#: Groups of spellings. Within a group every spelling denotes the same
#: region; across groups the regions differ.
EQUIVALENCE_GROUPS: tuple[tuple[str, ...], ...] = (
    (
        "SE(1400/0374)",
        f"SE(1400/0374({ALL_QUATERNARY}))",
        f"SE(1400/0374({','.join(reversed(QUATERNARY_CODES))}))",
        "SE(1400/374)",
        f"SE(1400/0374(1,2,3,4({ALL_QUINARY})))",
    ),
    (
        "SE(1400/0374(1,2))",
        "SE(1400/0374(2,1))",
        "SE(1400/0374(1),1400/0374(2))",
        "SE(1400/0374(1)),SE(1400/0374(2))",
    ),
    (
        "SE(1400/0374(3(A1)))",
        "SE(1400/0374(3(A1),3(A1)))",
    ),
    ("SE(1400/0374(1))",),
    ("NE(1400/0374)",),
    ("SE(1401/0374)",),
)


@pytest.mark.parametrize("group", EQUIVALENCE_GROUPS)
def test_req_13_equivalent_spellings_share_one_canonical_form(
    group: tuple[str, ...],
) -> None:
    canonical = {ix.normalize(spelling) for spelling in group}
    assert len(canonical) == 1


def test_req_13_distinct_regions_have_distinct_canonical_forms() -> None:
    representatives = [group[0] for group in EQUIVALENCE_GROUPS]
    canonical = [ix.normalize(index) for index in representatives]
    assert len(set(canonical)) == len(canonical)


def test_req_13_canonical_form_partitions_the_corpus_by_region() -> None:
    """Equal canonical form if and only if the same region, pairwise.

    Region identity is the declared grouping above, not something derived
    from ``normalize`` — deriving it would let the test agree with a
    broken canonicaliser. The "only if" half is what a weaker test misses:
    it is not enough that equivalent spellings agree, two different
    regions must also not be conflated.
    """
    labelled = [
        (group_id, spelling)
        for group_id, group in enumerate(EQUIVALENCE_GROUPS)
        for spelling in group
    ]
    for (left_group, left), (right_group, right) in combinations(labelled, 2):
        same_canonical = ix.normalize(left) == ix.normalize(right)
        assert same_canonical == (left_group == right_group), (left, right)


def test_req_13_canonical_form_survives_a_decompose_compose_detour() -> None:
    """Canonicalising before or after a round-trip gives the same string."""
    for index in CORPUS:
        detour = ix.compose(ix.decompose(index))
        assert ix.normalize(detour) == ix.normalize(index)


# --------------------------------------------------------------------------
# Criterion 5 -- the even/odd alphabets
# --------------------------------------------------------------------------


@pytest.mark.parametrize("resolution", list(_every_resolution()))
def test_every_resolution_accepts_exactly_its_own_alphabet(resolution: int) -> None:
    """The parser's alphabet at each level is refinement_alphabet's, exactly.

    Both directions are checked at once: every code of the level's own
    alphabet parses, and every code of the other alphabet is refused. That
    is the whole content of the inversion guard, stated once for all
    twelve levels.
    """
    own = refinement_alphabet(resolution)
    other = QUINARY_CODES if own is QUATERNARY_CODES else QUATERNARY_CODES
    for code in own:
        assert ix.is_valid_index(_chain_to(resolution - 1, tail=code))
    for code in other:
        assert not ix.is_valid_index(_chain_to(resolution - 1, tail=code))


@pytest.mark.parametrize("resolution", [2, 4, 6, 8, 10, 12])
def test_even_resolutions_are_quaternary(resolution: int) -> None:
    """Even resolutions take 1-4. The prototype had A1-E5 here."""
    assert refinement_alphabet(resolution) == QUATERNARY_CODES
    assert ix.is_valid_index(_chain_to(resolution - 1, tail="4"))


@pytest.mark.parametrize("resolution", [3, 5, 7, 9, 11, 13])
def test_odd_resolutions_are_quinary(resolution: int) -> None:
    """Odd resolutions from 3 take A1-E5. The prototype had 1-4 here."""
    assert refinement_alphabet(resolution) == QUINARY_CODES
    assert ix.is_valid_index(_chain_to(resolution - 1, tail="E5"))


def test_even_resolution_rejects_a_quinary_code() -> None:
    with pytest.raises(InvalidRefinementCodeError):
        ix.parse("SE(1400/0374(C2))")


def test_odd_resolution_rejects_a_quaternary_code() -> None:
    with pytest.raises(InvalidRefinementCodeError):
        ix.parse("SE(1400/0374(3(2)))")


def test_rejection_names_the_expected_alphabet() -> None:
    with pytest.raises(InvalidRefinementCodeError, match="A1..E5"):
        ix.parse("SE(1400/0374(3(2)))")


def test_invalid_refinement_code_is_an_invalid_index_error() -> None:
    """F2 phases downstream guard on the parent class, so the link matters."""
    assert issubclass(InvalidRefinementCodeError, InvalidIndexError)
    assert issubclass(InvalidIndexError, ITACaRTError)


@pytest.mark.parametrize("code", ["A6", "F1", "a1", "A0", "5", "0", "AA", ""])
def test_codes_outside_every_alphabet_are_refused(code: str) -> None:
    assert not ix.is_valid_index(f"SE(1400/0374(3({code})))")


# --------------------------------------------------------------------------
# Criterion 6 - deterministic decompose order
# --------------------------------------------------------------------------


def test_decompose_is_depth_first_left_to_right() -> None:
    index = "SE(1400/0374(1(A2,A1),2))"
    assert ix.decompose(index) == [
        "SE(1400/0374(1(A2)))",
        "SE(1400/0374(1(A1)))",
        "SE(1400/0374(2))",
    ]


def test_decompose_keeps_source_order_and_does_not_sort() -> None:
    """Ordering is normalize's business; decompose reports what is written."""
    index = "SE(1400/0374(3(E5,A1)))"
    assert ix.decompose(index) == [
        "SE(1400/0374(3(E5)))",
        "SE(1400/0374(3(A1)))",
    ]


@pytest.mark.parametrize("index", CORPUS)
def test_decompose_is_stable_across_calls(index: str) -> None:
    assert ix.decompose(index) == ix.decompose(index)


@pytest.mark.parametrize("index", CORPUS)
def test_iter_cells_yields_the_decompose_order(index: str) -> None:
    assert list(ix.iter_cells(index)) == ix.decompose(index)


def test_iter_cells_is_lazy(central_park_index: str) -> None:
    stream = ix.iter_cells(central_park_index)
    assert isinstance(stream, Iterator)
    assert next(stream) == "NW(0625/0451(1(E1(3(B2(4(A2)))))))"


def test_decompose_order_survives_a_compose_roundtrip(
    central_park_index: str,
) -> None:
    cells = ix.decompose(central_park_index)
    assert ix.decompose(ix.compose(cells)) == cells


# --------------------------------------------------------------------------
# parse
# --------------------------------------------------------------------------


def test_parse_exposes_the_documented_root_keys(paper_example_index: str) -> None:
    tree = ix.parse(paper_example_index)
    assert tree["quadrant"] == "SE"
    assert tree["base_cell"] == "1400/0374"
    assert isinstance(tree["children"], list)


def test_parse_root_is_the_globe_one_level_above_the_quadrants() -> None:
    tree = ix.parse("NE(0001/0001)")
    assert tree["resolution"] == ix.GLOBE_RESOLUTION == -1
    assert tree["code"] is None


def test_parse_labels_every_node_with_its_resolution(
    paper_example_index: str,
) -> None:
    node = ix.parse(paper_example_index)
    seen = []
    while True:
        children = node["children"]
        assert isinstance(children, list)
        if not children:
            break
        node = children[0]
        seen.append((node["code"], node["resolution"]))
    assert seen == [
        ("SE", 0),
        ("1400/0374", 1),
        ("3", 2),
        ("C2", 3),
        ("3", 4),
    ]


def test_parse_reports_no_single_quadrant_when_several_are_addressed() -> None:
    tree = ix.parse("NE(0001/0001),SW(0002/0002)")
    assert tree["quadrant"] is None
    assert tree["base_cell"] is None


def test_parse_reports_no_single_base_cell_when_several_are_addressed() -> None:
    tree = ix.parse("NE(0001/0001,0002/0002)")
    assert tree["quadrant"] == "NE"
    assert tree["base_cell"] is None


# --------------------------------------------------------------------------
# Syntax errors
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "index",
    [
        "",
        "   ",
        "XX(0001/0001)",
        "N(0001/0001)",
        "SE(0001/0001",
        "SE0001/0001)",
        "SE(0001/0001))",
        "SE()",
        "SE(0001/0001(,1))",
        "SE(0001/0001)trailing",
        "SE(0001/0001),",
        "SE(0001)",
        "SE(0001/)",
        "SE(/0001)",
        "SE(0001/0001/0002)",
        "SE(00a1/0001)",
        "SE(0001/0001( 1))",
    ],
)
def test_malformed_strings_are_refused(index: str) -> None:
    assert not ix.is_valid_index(index)
    with pytest.raises(InvalidIndexError):
        ix.parse(index)


def test_unknown_quadrant_raises_the_specific_subclass() -> None:
    with pytest.raises(InvalidQuadrantError):
        ix.parse("XX(0001/0001)")


def test_x_index_admits_one_column_past_the_table_1_maximum() -> None:
    """The grammar reaches 2004, and stops there.

    Table 1 caps the column at 2003, and every resolution-1 cell obeys
    that. The grammar stands one column further out because refining a
    trapezoidal cell of column 2003 produces a child whose resolution-1
    prefix is 2004 -- the only place in the system where a child's prefix
    is not its parent. That prefix has to be spellable.

    The ceiling is syntactic. It says nothing about which cells exist,
    which ``boundary.is_valid_cell`` decides, and 2004 is never a
    resolution-1 cell in its own right.
    """
    table_1_max = int(RES1_MAX_INDEX.split("/")[0])
    assert ix.is_valid_index(f"SE({table_1_max}/0001)")
    assert ix.is_valid_index(f"SE({table_1_max + 1}/0001)")
    assert not ix.is_valid_index(f"SE({table_1_max + 2}/0001)")


def test_the_x_ceiling_stands_exactly_one_column_past_table_1() -> None:
    """A regression guard on the size of the concession.

    Loosening the ceiling is a deliberate exception, and an exception that
    can drift is not an exception. Nothing measured anywhere on the globe
    reaches 2005: the widest parallel is the equator at 2003.7508 cell
    sides, and the extension zones never approach it because ``cos φ``
    withdraws faster than the extension adds.
    """
    assert ix._RES1_MAX_X == int(RES1_MAX_INDEX.split("/")[0]) + 1


def test_y_index_beyond_the_table_1_maximum_is_refused() -> None:
    max_y = int(RES1_MAX_INDEX.split("/")[1])
    assert ix.is_valid_index(f"SE(0001/{max_y})")
    assert not ix.is_valid_index(f"SE(0001/{max_y + 1})")


def test_unicode_digits_are_not_accepted_as_base_indices() -> None:
    """``"\u00b2".isdigit()`` is True but ``int()`` refuses it."""
    assert not ix.is_valid_index("SE(\u00b2/0001)")


def test_descending_past_the_finest_resolution_is_refused() -> None:
    assert ix.is_valid_index(_chain_to(MAX_RESOLUTION))
    with pytest.raises(InvalidIndexError, match="past resolution"):
        ix.parse(_chain_to(MAX_RESOLUTION, tail="1"))


def test_a_non_string_is_refused_without_a_bare_type_error() -> None:
    with pytest.raises(InvalidIndexError):
        ix.parse(None)  # type: ignore[arg-type]


def test_surrounding_whitespace_is_tolerated(paper_example_index: str) -> None:
    assert ix.normalize(f"  {paper_example_index}\n") == paper_example_index


# --------------------------------------------------------------------------
# is_valid_index / is_atomic / count_cells
# --------------------------------------------------------------------------


@pytest.mark.parametrize("index", CORPUS)
def test_corpus_is_structurally_valid(index: str) -> None:
    assert ix.is_valid_index(index)


def test_a_bare_quadrant_addresses_the_whole_quadrant() -> None:
    """A node without a subtree is the whole cell, resolution 0 included."""
    assert ix.is_valid_index("NE")
    assert ix.is_atomic("NE")
    assert ix.decompose("NE") == ["NE"]
    assert ix.split_components("NE") == ["NE"]


def test_is_atomic_is_false_for_a_region(central_park_index: str) -> None:
    assert not ix.is_atomic(central_park_index)


def test_is_atomic_raises_rather_than_calling_a_malformed_string_a_region() -> None:
    with pytest.raises(InvalidIndexError):
        ix.is_atomic("SE(0001/0001")


def test_count_cells_matches_decompose_over_the_corpus() -> None:
    for index in CORPUS:
        assert ix.count_cells(index) == len(ix.decompose(index))


# --------------------------------------------------------------------------
# compose
# --------------------------------------------------------------------------


def test_compose_factors_out_shared_ancestry() -> None:
    cells = ["SE(1400/0374(1(A1)))", "SE(1400/0374(1(A2)))"]
    assert ix.compose(cells) == "SE(1400/0374(1(A1,A2)))"


def test_compose_keeps_first_appearance_order() -> None:
    cells = ["SE(1400/0374(1(A2)))", "SE(1400/0374(1(A1)))"]
    assert ix.compose(cells) == "SE(1400/0374(1(A2,A1)))"


def test_compose_does_not_collapse_a_complete_set() -> None:
    """Collapse belongs to normalize; compose must not silently rewrite."""
    cells = [f"SE(1400/0374({code}))" for code in QUATERNARY_CODES]
    assert ix.compose(cells) == f"SE(1400/0374({ALL_QUATERNARY}))"


def test_compose_accepts_cells_of_different_resolutions() -> None:
    cells = ["SE(1400/0374(1))", "SE(1400/0374(2(A1)))"]
    assert ix.compose(cells) == "SE(1400/0374(1,2(A1)))"


def test_compose_accepts_cells_of_different_quadrants() -> None:
    cells = ["NE(0001/0001)", "SW(0002/0002)"]
    composed = ix.compose(cells)
    assert composed == "NE(0001/0001),SW(0002/0002)"
    assert ix.decompose(composed) == cells


def test_compose_accepts_a_generator(central_park_index: str) -> None:
    """iter_cells exists to avoid materialising; compose must take it."""
    assert ix.compose(ix.iter_cells(central_park_index)) == central_park_index


def test_compose_accepts_a_compositional_entry() -> None:
    cells = ["SE(1400/0374(1(A1,A2)))", "SE(1400/0374(2))"]
    assert ix.compose(cells) == "SE(1400/0374(1(A1,A2),2))"


def test_compose_lets_a_whole_cell_absorb_a_descendant() -> None:
    assert ix.compose(["SE(1400/0374(1))", "SE(1400/0374(1(A1)))"]) == (
        "SE(1400/0374(1))"
    )


def test_compose_of_nothing_is_an_error_not_an_empty_string() -> None:
    with pytest.raises(InvalidIndexError):
        ix.compose([])


def test_compose_propagates_a_malformed_entry() -> None:
    with pytest.raises(InvalidIndexError):
        ix.compose(["SE(1400/0374(1))", "not an index"])


# --------------------------------------------------------------------------
# split_components / quadrant_of / base_cell_of
# --------------------------------------------------------------------------


def test_split_components_refuses_a_region(central_park_index: str) -> None:
    with pytest.raises(NonAtomicIndexError):
        ix.split_components(central_park_index)


def test_split_components_of_a_two_cell_index_names_the_count() -> None:
    with pytest.raises(NonAtomicIndexError, match="2 cells"):
        ix.split_components("SE(1400/0374(1,2))")


def test_quadrant_of_a_single_cell_is_a_scalar(paper_example_index: str) -> None:
    assert ix.quadrant_of(paper_example_index) == "SE"


def test_quadrant_of_a_region_is_positionally_aligned() -> None:
    index = "SE(1400/0374(1,2)),NE(0001/0001)"
    assert ix.quadrant_of(index) == ["SE", "SE", "NE"]
    assert len(ix.quadrant_of(index)) == ix.count_cells(index)


def test_base_cell_of_is_fully_qualified(paper_example_index: str) -> None:
    assert ix.base_cell_of(paper_example_index) == "SE(1400/0374)"


def test_base_cell_of_a_region_is_positionally_aligned() -> None:
    index = "SE(1400/0374(1),1401/0374(2))"
    assert ix.base_cell_of(index) == ["SE(1400/0374)", "SE(1401/0374)"]


def test_base_cell_of_aligns_with_decompose_on_figure_7(
    central_park_index: str,
) -> None:
    bases = ix.base_cell_of(central_park_index)
    assert isinstance(bases, list)
    assert len(bases) == len(ix.decompose(central_park_index))
    assert set(bases) == {"NW(0625/0451)"}


def test_base_cell_of_a_whole_quadrant_has_no_answer() -> None:
    with pytest.raises(InvalidIndexError, match="depth 1"):
        ix.base_cell_of("NE")


# --------------------------------------------------------------------------
# Constants exposed by the module
# --------------------------------------------------------------------------


def test_level_constants_agree_with_the_specification() -> None:
    assert ix.QUADRANT_RESOLUTION == 0
    assert ix.BASE_CELL_RESOLUTION == 1
    assert ix.GLOBE_RESOLUTION == -1


def test_public_surface_matches_the_contract() -> None:
    assert set(ix.__all__) == {
        "parse",
        "is_valid_index",
        "is_atomic",
        "decompose",
        "compose",
        "normalize",
        "count_cells",
        "iter_cells",
        "split_components",
        "quadrant_of",
        "base_cell_of",
    }
    for name in ix.__all__:
        assert callable(getattr(ix, name))
