"""Unit tests for :mod:`itacart.resolutions`.

Covers acceptance criteria 6 and 7 of F3.

Table 1 is checked by two independent routes, as criterion 7 requires.
The literal transcription below is a *second witness*: it is typed from
the published table rather than read from :mod:`itacart.constants`, so a
mistyped constant fails here instead of agreeing with itself. The closed
forms are a third route, and the fixed ratio of 5 between the two scale
families is a fourth.
"""

from __future__ import annotations

import math

import pytest

from itacart import resolutions as res
from itacart.exceptions import ITACaRTError, ResolutionError

# --------------------------------------------------------------------------
# Table 1, transcribed from the paper. Not imported from constants.
# (resolution, side_m, area_m2, children, visualization_scale, analysis_scale)
# --------------------------------------------------------------------------

TABLE_1 = (
    (1, 10_000.0, 1e8, None, 100_000_000, 20_000_000),
    (2, 5_000.0, 2.5e7, 4, 50_000_000, 10_000_000),
    (3, 1_000.0, 1e6, 25, 10_000_000, 2_000_000),
    (4, 500.0, 250_000.0, 4, 5_000_000, 1_000_000),
    (5, 100.0, 10_000.0, 25, 1_000_000, 200_000),
    (6, 50.0, 2_500.0, 4, 500_000, 100_000),
    (7, 10.0, 100.0, 25, 100_000, 20_000),
    (8, 5.0, 25.0, 4, 50_000, 10_000),
    (9, 1.0, 1.0, 25, 10_000, 2_000),
    (10, 0.5, 0.25, 4, 5_000, 1_000),
    (11, 0.1, 0.01, 25, 1_000, 200),
    (12, 0.05, 0.0025, 4, 500, 100),
    (13, 0.01, 0.0001, 25, 100, 20),
)

METRIC_RESOLUTIONS = tuple(row[0] for row in TABLE_1)
REFINED_RESOLUTIONS = tuple(r for r in METRIC_RESOLUTIONS if r >= 2)

MINIMUM_VISIBLE_LINE_M = 1e-4  # Jenny et al. (2008), 0.1 mm on paper
SAMPLING_INTERVAL_M = 5e-4  # Tobler (1987), 0.5 mm on paper


# ==========================================================================
# Criterion 6 -- nominal_cell_area against Table 1
# ==========================================================================


@pytest.mark.parametrize("resolution,_side,area,_c,_v,_a", TABLE_1)
def test_criterion_6_nominal_area_matches_the_published_table(
    resolution: int, _side: float, area: float, _c: int | None, _v: int, _a: int
) -> None:
    """Route 1: literal transcription of the 13 metric rows."""
    assert res.nominal_cell_area(resolution) == area


@pytest.mark.parametrize("resolution", METRIC_RESOLUTIONS)
def test_criterion_6_area_equals_side_squared(resolution: int) -> None:
    """Route 2: the equal-area property stated as a table invariant.

    The cell is a parallelogram whose base and height are equal, so its
    area is the square of its side at every resolution. Holding on all 13
    rows is what makes the area column an assertion about the geometry
    rather than about the typing.
    """
    side = res.cell_size(resolution)
    assert res.nominal_cell_area(resolution) == pytest.approx(side * side, rel=1e-12)


@pytest.mark.parametrize("resolution,side,_a,_c,_v,_an", TABLE_1)
def test_cell_size_matches_the_published_table(
    resolution: int, side: float, _a: float, _c: int | None, _v: int, _an: int
) -> None:
    assert res.cell_size(resolution) == side


def test_resolution_zero_has_no_size_or_area() -> None:
    """Table 1 prints a dash for the quadrant row, not a number."""
    with pytest.raises(ResolutionError, match="quadrant"):
        res.cell_size(0)
    with pytest.raises(ResolutionError, match="quadrant"):
        res.nominal_cell_area(0)


@pytest.mark.parametrize("resolution", [-1, 14, 99])
def test_size_and_area_reject_out_of_range_resolutions(resolution: int) -> None:
    with pytest.raises(ResolutionError, match="outside"):
        res.cell_size(resolution)
    with pytest.raises(ResolutionError, match="outside"):
        res.nominal_cell_area(resolution)


@pytest.mark.parametrize("bogus", [True, False, 1.0, "9", None])
def test_size_rejects_non_integer_resolutions(bogus: object) -> None:
    """``bool`` is an ``int`` and would silently address resolution 1."""
    with pytest.raises(ResolutionError, match="must be an int"):
        res.cell_size(bogus)  # type: ignore[arg-type]


# ==========================================================================
# Criterion 7 -- scale_for_resolution against both columns, two routes
# ==========================================================================


@pytest.mark.parametrize("resolution,_s,_a,_c,visualization,analysis", TABLE_1)
def test_criterion_7_scales_match_the_published_table(
    resolution: int,
    _s: float,
    _a: float,
    _c: int | None,
    visualization: int,
    analysis: int,
) -> None:
    """Route 1: literal transcription of both scale columns, 13 rows."""
    assert res.scale_for_resolution(resolution, "visualization") == visualization
    assert res.scale_for_resolution(resolution, "analysis") == analysis


@pytest.mark.parametrize("resolution", METRIC_RESOLUTIONS)
def test_criterion_7_scales_match_their_closed_forms(resolution: int) -> None:
    """Route 2: the closed forms the two cited rules imply.

    Visualization divides the side by the 0.1 mm minimum visible line;
    analysis divides it by the 0.5 mm sampling interval. Agreeing with
    route 1 on all 13 rows is an assertion about the table; either route
    alone is an assertion about the typing.
    """
    side = res.cell_size(resolution)
    assert res.scale_for_resolution(resolution, "visualization") == pytest.approx(
        side / MINIMUM_VISIBLE_LINE_M, rel=1e-12
    )
    assert res.scale_for_resolution(resolution, "analysis") == pytest.approx(
        side / SAMPLING_INTERVAL_M, rel=1e-12
    )


@pytest.mark.parametrize("resolution", METRIC_RESOLUTIONS)
def test_criterion_7_the_two_scale_families_differ_by_a_factor_of_five(
    resolution: int,
) -> None:
    """Route 3, free: 0.5 mm over 0.1 mm is 5 at every resolution."""
    visualization = res.scale_for_resolution(resolution, "visualization")
    analysis = res.scale_for_resolution(resolution, "analysis")
    assert visualization == 5 * analysis


def test_visualization_is_the_default_scale_family() -> None:
    assert res.scale_for_resolution(9) == res.scale_for_resolution(9, "visualization")


def test_scale_rejects_an_unknown_family() -> None:
    with pytest.raises(ResolutionError, match="unknown scale kind"):
        res.scale_for_resolution(9, "cadastral")  # type: ignore[arg-type]


def test_resolution_zero_has_no_scale() -> None:
    with pytest.raises(ResolutionError, match="quadrant"):
        res.scale_for_resolution(0)


# --------------------------------------------------------------------------
# resolution_for_scale -- the inverse
# --------------------------------------------------------------------------


@pytest.mark.parametrize("resolution", METRIC_RESOLUTIONS)
@pytest.mark.parametrize("kind", ["visualization", "analysis"])
def test_resolution_for_scale_inverts_scale_for_resolution(
    resolution: int, kind: res.ScaleKind
) -> None:
    """Round trip on every printed denominator, in both families."""
    denominator = res.scale_for_resolution(resolution, kind)
    assert res.resolution_for_scale(denominator, kind) == resolution


def test_between_rows_the_two_families_answer_differently() -> None:
    """The two scale families, stated as the example that motivates them.

    At 1:2 000 a resolution-11 cell would draw at 0.05 mm and vanish, so
    the drawable answer is 10. At 1:3 000 the analysis rule asks for cells
    at or below 1.5 m on the ground, so the coarsest adequate level is 9.
    """
    assert res.resolution_for_scale(2_000, "visualization") == 10
    assert res.resolution_for_scale(3_000, "analysis") == 9


def test_resolution_for_scale_rejects_targets_no_level_satisfies() -> None:
    with pytest.raises(ResolutionError, match="drawable"):
        res.resolution_for_scale(200_000_000, "visualization")
    with pytest.raises(ResolutionError, match="fine enough"):
        res.resolution_for_scale(10, "analysis")


@pytest.mark.parametrize("bogus", [0, -1])
def test_resolution_for_scale_rejects_non_positive_denominators(bogus: int) -> None:
    with pytest.raises(ResolutionError, match="positive"):
        res.resolution_for_scale(bogus)


@pytest.mark.parametrize("bogus", [True, 1_000.0, "1000", None])
def test_resolution_for_scale_rejects_non_integer_denominators(bogus: object) -> None:
    with pytest.raises(ResolutionError, match="must be an int"):
        res.resolution_for_scale(bogus)  # type: ignore[arg-type]


def test_resolution_for_scale_rejects_an_unknown_family() -> None:
    with pytest.raises(ResolutionError, match="unknown scale kind"):
        res.resolution_for_scale(1_000, "cadastral")  # type: ignore[arg-type]


# ==========================================================================
# Refinement ratios -- the area/linear distinction
# ==========================================================================


@pytest.mark.parametrize("resolution,_s,_a,children,_v,_an", TABLE_1)
def test_refinement_ratio_matches_the_published_table(
    resolution: int,
    _s: float,
    _a: float,
    children: int | None,
    _v: int,
    _an: int,
) -> None:
    if children is None:
        pytest.skip("resolution 1 is the Cartesian base grid, not a refinement")
    assert res.refinement_ratio(resolution) == children


@pytest.mark.parametrize("resolution", REFINED_RESOLUTIONS)
def test_linear_ratio_is_the_square_root_of_the_area_ratio(resolution: int) -> None:
    """The property, stated once, rather than a second even/odd table."""
    linear = res.linear_refinement_ratio(resolution)
    assert linear * linear == res.refinement_ratio(resolution)


@pytest.mark.parametrize("resolution", REFINED_RESOLUTIONS)
def test_linear_ratio_is_the_actual_side_ratio_of_the_table(resolution: int) -> None:
    """The area/linear distinction measured, not asserted.

    This is the check that would catch a descent written against the area
    ratio: the side of a child really does shrink by 2 or by 5, and the
    published side column says so on all 12 refined levels.
    """
    parent_side = res.cell_size(resolution - 1)
    child_side = res.cell_size(resolution)
    assert parent_side / child_side == pytest.approx(
        res.linear_refinement_ratio(resolution), rel=1e-12
    )


def test_the_two_ratios_are_never_equal() -> None:
    """Substituting one for the other is always a real error, never a no-op."""
    for resolution in REFINED_RESOLUTIONS:
        assert res.refinement_ratio(resolution) != res.linear_refinement_ratio(
            resolution
        )


@pytest.mark.parametrize("resolution", [0, 1])
def test_refinement_ratios_reject_the_unrefined_levels(resolution: int) -> None:
    with pytest.raises(ResolutionError, match="not produced by a refinement"):
        res.refinement_ratio(resolution)
    with pytest.raises(ResolutionError, match="not produced by a refinement"):
        res.linear_refinement_ratio(resolution)


@pytest.mark.parametrize("resolution", REFINED_RESOLUTIONS)
def test_even_levels_are_quaternary_and_odd_levels_quinary(resolution: int) -> None:
    expected = 4 if resolution % 2 == 0 else 25
    assert res.refinement_ratio(resolution) == expected


# ==========================================================================
# get_resolution
# ==========================================================================


def test_get_resolution_of_a_quadrant_is_zero() -> None:
    assert res.get_resolution("NE") == 0


def test_get_resolution_of_a_base_cell_is_one() -> None:
    assert res.get_resolution("NE(0001/0002)") == 1


def test_get_resolution_of_the_paper_example(paper_example_index: str) -> None:
    """Section 3.1 of the paper: quadrant, base pair, then three codes."""
    assert res.get_resolution(paper_example_index) == 4


def test_get_resolution_accepts_a_uniform_compositional_index() -> None:
    assert res.get_resolution("NE(0001/0002(1(A1,A2,B3)))") == 3


def test_get_resolution_rejects_a_mixed_resolution_index() -> None:
    """decompose() does not uniformize; neither does this (invariant F2)."""
    with pytest.raises(ResolutionError, match="mixed-resolution"):
        res.get_resolution("NE(0001/0002(1(A1,A2),2))")


def test_get_resolution_of_figure_7_reports_the_mixture(
    central_park_index: str,
) -> None:
    """Figure 7(a) has leaves at resolutions 6 and 7 simultaneously."""
    with pytest.raises(ResolutionError, match="mixed-resolution"):
        res.get_resolution(central_park_index)


# ==========================================================================
# effective_cell_area
# ==========================================================================


def test_effective_cell_area_reads_the_table_for_unclipped_cells() -> None:
    """Parallelograms and meridian triangles carry the nominal area exactly.

    The triangle is the interesting half: base twice height, halved, is
    the area of the parallelogram it replaces, which is why the prime
    meridian does not cost the grid its equal-area property.
    """
    for cell in ("NE(1400/0374)", "SE(0900/0200(3))", "NE(0000/0500)"):
        assert res.effective_cell_area(cell) == res.nominal_cell_area(
            res.get_resolution(cell)
        )


def test_effective_cell_area_measures_a_clipped_cell() -> None:
    """A trapezoid is measured from its own vertices, and differs."""
    from itacart import boundary

    cell = "NE(2003/0000)"
    assert boundary.cell_shape(cell) == "trapezoid"
    effective = res.effective_cell_area(cell)
    assert effective != res.nominal_cell_area(1)
    assert effective > 0.0


def test_effective_cell_area_is_positionally_aligned() -> None:
    """A composed index answers a list, one entry per terminal cell."""
    areas = res.effective_cell_area("NE(1400/0374(1,2,3,4))")
    assert areas == [res.nominal_cell_area(2)] * 4


def test_effective_cell_area_does_not_masquerade_as_a_domain_error() -> None:
    """An unbuilt function must not be swallowed by ``except ITACaRTError``.

    Same family of mistake as the earlier transcription errors, in the
    other direction:
    there the package leaked a bare exception, here it must refuse to
    capture one.
    """
    assert not issubclass(NotImplementedError, ITACaRTError)


# ==========================================================================
# resolution_table and tokenization
# ==========================================================================


def test_resolution_table_has_one_record_per_level_in_order() -> None:
    table = res.resolution_table()
    assert [row["resolution"] for row in table] == list(range(14))


def test_resolution_table_carries_the_documented_keys() -> None:
    expected = {
        "resolution",
        "cell_size_m",
        "cell_area_m2",
        "refinement",
        "index_alphabet",
        "visualization_scale",
        "analysis_scale",
    }
    assert all(set(row) == expected for row in res.resolution_table())


@pytest.mark.parametrize(
    "resolution,side,area,children,visualization,analysis", TABLE_1
)
def test_resolution_table_agrees_with_the_published_table(
    resolution: int,
    side: float,
    area: float,
    children: int | None,
    visualization: int,
    analysis: int,
) -> None:
    row = res.resolution_table()[resolution]
    assert row["cell_size_m"] == side
    assert row["cell_area_m2"] == area
    assert row["refinement"] == children
    assert row["visualization_scale"] == visualization
    assert row["analysis_scale"] == analysis


def test_resolution_table_row_zero_is_the_quadrant_row() -> None:
    row = res.resolution_table()[0]
    assert row["cell_size_m"] is None
    assert row["cell_area_m2"] is None
    assert row["index_alphabet"] == ("NE", "NW", "SE", "SW")


def test_resolution_table_is_freshly_built_on_each_call() -> None:
    """Mutating the result must not corrupt the constants behind it."""
    first = res.resolution_table()
    first[9]["cell_size_m"] = -1.0
    assert res.resolution_table()[9]["cell_size_m"] == 1.0


@pytest.mark.parametrize("resolution", METRIC_RESOLUTIONS)
def test_tokenizable_levels_are_exactly_those_with_decimal_areas(
    resolution: int,
) -> None:
    """The property behind the constant, not a second copy of the list.

    A level is tokenizable when its side is an exact power of ten metres,
    which makes its area an exact power of ten square metres.
    """
    exponent = math.log10(res.cell_size(resolution))
    decimal_side = math.isclose(exponent, round(exponent), abs_tol=1e-12)
    assert res.is_tokenizable_resolution(resolution) is decimal_side


def test_resolution_thirteen_is_exactly_one_square_centimetre() -> None:
    assert res.is_tokenizable_resolution(13)
    assert res.nominal_cell_area(13) == pytest.approx(1e-4, rel=1e-12)


def test_tokenizable_rejects_the_quadrant_level() -> None:
    with pytest.raises(ResolutionError, match="quadrant"):
        res.is_tokenizable_resolution(0)


# ==========================================================================
# Fidelity to itacart_core, and the two deliberate divergences
# ==========================================================================

# itacart_core/res.py stores sides as exact integer centimetres.
# Transcribed here so the port is measured against the origin as well as
# against the paper -- a third witness, independent of both.
ORIGIN_BASE_LENGTH_CM = {
    1: 1_000_000,
    2: 500_000,
    3: 100_000,
    4: 50_000,
    5: 10_000,
    6: 5_000,
    7: 1_000,
    8: 500,
    9: 100,
    10: 50,
    11: 10,
    12: 5,
    13: 1,
}


@pytest.mark.parametrize("resolution", METRIC_RESOLUTIONS)
def test_sides_and_areas_reproduce_the_origin_exactly(resolution: int) -> None:
    """Zero ulps, not a tolerance.

    The origin derives every figure from an exact integer centimetre
    count, so agreement here is exact by construction and any drift is a
    real transcription error rather than a libm difference.
    """
    side_cm = ORIGIN_BASE_LENGTH_CM[resolution]
    assert res.cell_size(resolution) == side_cm / 100.0
    assert res.nominal_cell_area(resolution) == (side_cm * side_cm) / 10_000.0


def test_divergence_resolution_one_is_tokenizable_here_but_not_in_the_origin() -> None:
    """Pinned so the divergence cannot drift back unnoticed.

    ``itacart_core`` enumerates only the odd levels. Resolution 1 has a
    side of 10^4 m and an area of 10^8 m2, so it meets the stated
    property; the paper enumerates neither set. F0 owns this reading.
    """
    assert res.is_tokenizable_resolution(1) is True
    assert 1 not in {3, 5, 7, 9, 11, 13}


def test_divergence_the_linear_ratio_refuses_resolution_one() -> None:
    """The origin returns 1 here; returning 1 is the trap.

    ``subdivisions_per_axis(1) == 1`` lets a descent loop divide the side
    by one and emit a cell one level shallower than asked, with nothing
    raised. Refusing is the only answer that fails loudly.
    """
    with pytest.raises(ResolutionError, match="not produced by a refinement"):
        res.linear_refinement_ratio(1)
