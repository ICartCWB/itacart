"""Tests for :mod:`itacart.constants`.

The point of this file is auditability. Table 1 of the paper is restated
here as literals, independently of the module, so that the test fails if
either copy drifts. Everything else is checked as an invariant that the
table must satisfy, which catches a typo the literal comparison would
happily reproduce if both were wrong in the same way.
"""

from __future__ import annotations

import ast
import math
import pathlib
from typing import cast

import pytest

from itacart import constants as c
from itacart.exceptions import ResolutionError

# --------------------------------------------------------------------------
# Table 1 of the paper, transcribed independently of the module.
# (resolution, side_m, area_m2, alphabet_label, visualization, analysis)
# --------------------------------------------------------------------------

Row = tuple[int, float | None, float | None, str, int | None, int | None]
MetricRow = tuple[int, float, float, str, int, int]

TABLE_1: list[Row] = [
    (0, None, None, "quadrant", None, None),
    (1, 10_000.0, 100e6, "xxxx/yyyy", 100_000_000, 20_000_000),
    (2, 5_000.0, 25e6, "1..4", 50_000_000, 10_000_000),
    (3, 1_000.0, 1e6, "A1..E5", 10_000_000, 2_000_000),
    (4, 500.0, 250_000.0, "1..4", 5_000_000, 1_000_000),
    (5, 100.0, 10_000.0, "A1..E5", 1_000_000, 200_000),
    (6, 50.0, 2_500.0, "1..4", 500_000, 100_000),
    (7, 10.0, 100.0, "A1..E5", 100_000, 20_000),
    (8, 5.0, 25.0, "1..4", 50_000, 10_000),
    (9, 1.0, 1.0, "A1..E5", 10_000, 2_000),
    (10, 0.5, 0.25, "1..4", 5_000, 1_000),
    (11, 0.1, 0.01, "A1..E5", 1_000, 200),
    (12, 0.05, 0.0025, "1..4", 500, 100),
    (13, 0.01, 0.0001, "A1..E5", 100, 20),
]

METRIC_ROWS: list[MetricRow] = [cast(MetricRow, row) for row in TABLE_1 if row[0] > 0]
"""The 13 rows carrying a metric size; resolution 0 has none."""


def _alphabet_label(resolution: int) -> str:
    if resolution == 0:
        return "quadrant"
    if resolution == 1:
        return "xxxx/yyyy"
    if c.REFINEMENT_ALPHABET[resolution] is c.QUATERNARY_CODES:
        return "1..4"
    return "A1..E5"


# --------------------------------------------------------------------------
# Criterion 4 -- every constant of Table 1 checks out, value by value
# --------------------------------------------------------------------------


@pytest.mark.parametrize("row", TABLE_1, ids=[f"res{r[0]}" for r in TABLE_1])
def test_table1_row_by_row(row: Row) -> None:
    """Each of the 14 rows of Table 1, checked field by field."""
    resolution, side, area, alphabet, visualization, analysis = row
    assert c.CELL_SIZE_M[resolution] == side
    assert c.CELL_AREA_M2[resolution] == area
    assert c.VISUALIZATION_SCALE[resolution] == visualization
    assert c.ANALYSIS_SCALE[resolution] == analysis
    assert _alphabet_label(resolution) == alphabet


def test_table1_has_fourteen_rows() -> None:
    assert c.RESOLUTION_COUNT == 14
    assert (c.MIN_RESOLUTION, c.MAX_RESOLUTION) == (0, 13)
    for table in (
        c.CELL_SIZE_M,
        c.CELL_AREA_M2,
        c.VISUALIZATION_SCALE,
        c.ANALYSIS_SCALE,
        c.REFINEMENT_RATIO,
        c.REFINEMENT_ALPHABET,
        c.RESOLUTION_TABLE,
    ):
        assert len(table) == 14


def test_resolution_zero_has_no_metric_size() -> None:
    """Resolution 0 is a global quadrant: no side, no area, no scale."""
    assert c.CELL_SIZE_M[0] is None
    assert c.CELL_AREA_M2[0] is None
    assert c.VISUALIZATION_SCALE[0] is None
    assert c.ANALYSIS_SCALE[0] is None
    assert c.RESOLUTION_TABLE[0].alphabet == c.QUADRANTS


# --------------------------------------------------------------------------
# Cross-invariants of Table 1
# --------------------------------------------------------------------------


def test_i1_refinement_ratio_alternates() -> None:
    """I-1: side divides by 2 into even resolutions, by 5 into odd ones."""
    for resolution in range(2, 14):
        previous = c.CELL_SIZE_M[resolution - 1]
        current = c.CELL_SIZE_M[resolution]
        assert previous is not None and current is not None
        expected = 2.0 if resolution % 2 == 0 else 5.0
        assert previous / current == pytest.approx(expected, rel=1e-12)


def test_i1_decimal_decade_every_two_levels() -> None:
    """Consequence of I-1: an exact decimal decade every two levels."""
    for resolution in range(3, 14):
        two_up = c.CELL_SIZE_M[resolution - 2]
        current = c.CELL_SIZE_M[resolution]
        assert two_up is not None and current is not None
        assert two_up / current == pytest.approx(10.0, rel=1e-12)


def test_i2_visualization_scale_from_minimum_visible_line() -> None:
    """I-2: scale = side / 1e-4 (Jenny et al., 2008; 0.1 mm line)."""
    for resolution, side, _, _, visualization, _ in METRIC_ROWS:
        assert c.VISUALIZATION_SCALE[resolution] == pytest.approx(side / 1e-4)
        assert visualization == pytest.approx(side / 1e-4)


def test_i3_analysis_scale_from_sampling_theory() -> None:
    """I-3: scale = side * 2 * 1000 (Tobler, 1987)."""
    for resolution, side, _, _, _, analysis in METRIC_ROWS:
        assert c.ANALYSIS_SCALE[resolution] == pytest.approx(side * 2 * 1000)
        assert analysis == pytest.approx(side * 2 * 1000)


def test_i4_area_is_the_square_of_the_side() -> None:
    """Base and height are equal, so nominal area is the side squared.

    Not listed in the opening package, but it holds on all 13 metric rows
    and it is the table-level statement of the equal-area property.
    """
    for resolution, side, area, _, _, _ in METRIC_ROWS:
        assert c.CELL_AREA_M2[resolution] == pytest.approx(side * side, rel=1e-12)
        assert area == pytest.approx(side * side, rel=1e-12)


def test_resolution_table_agrees_with_the_component_tuples() -> None:
    """RESOLUTION_TABLE is assembled, never restated."""
    for spec in c.RESOLUTION_TABLE:
        resolution = spec.resolution
        assert spec.cell_size_m == c.CELL_SIZE_M[resolution]
        assert spec.cell_area_m2 == c.CELL_AREA_M2[resolution]
        assert spec.refinement_ratio == c.REFINEMENT_RATIO[resolution]
        assert spec.visualization_scale == c.VISUALIZATION_SCALE[resolution]
        assert spec.analysis_scale == c.ANALYSIS_SCALE[resolution]
    assert [spec.resolution for spec in c.RESOLUTION_TABLE] == list(range(14))


# --------------------------------------------------------------------------
# Criterion 5 -- QUINARY_CODES, and the inversion trap
# --------------------------------------------------------------------------


def test_quinary_codes_are_25_row_major_from_a1_to_e5() -> None:
    assert len(c.QUINARY_CODES) == 25
    assert len(set(c.QUINARY_CODES)) == 25
    assert c.QUINARY_CODES[0] == "A1"
    assert c.QUINARY_CODES[-1] == "E5"
    expected = tuple(f"{letter}{digit}" for letter in "ABCDE" for digit in "12345")
    assert c.QUINARY_CODES == expected


def test_quinary_codes_row_major_not_column_major() -> None:
    """Row-major: the sixth code opens a new row, it does not continue one.

    Column-major would give ``B1`` at position 1 instead of position 5.
    """
    assert c.QUINARY_CODES[1] == "A2"
    assert c.QUINARY_CODES[5] == "B1"
    assert c.QUINARY_CODES[24] == "E5"


def test_quaternary_codes() -> None:
    assert c.QUATERNARY_CODES == ("1", "2", "3", "4")
    assert c.QUATERNARY_GRID_SIZE == 2
    assert c.QUINARY_GRID_SIZE == 5


def test_even_is_quaternary_odd_is_quinary() -> None:
    """Regression guard: the two alphabets must not be inverted."""
    for resolution in range(2, 14):
        alphabet = c.refinement_alphabet(resolution)
        if resolution % 2 == 0:
            assert alphabet == c.QUATERNARY_CODES, f"res {resolution}"
            assert c.REFINEMENT_RATIO[resolution] == 4
        else:
            assert alphabet == c.QUINARY_CODES, f"res {resolution}"
            assert c.REFINEMENT_RATIO[resolution] == 25


def test_refinement_ratio_matches_alphabet_length() -> None:
    for resolution in range(2, 14):
        ratio = c.REFINEMENT_RATIO[resolution]
        assert ratio == len(c.refinement_alphabet(resolution))


@pytest.mark.parametrize("resolution", [0, 1])
def test_refinement_alphabet_rejects_non_refined_resolutions(
    resolution: int,
) -> None:
    """Resolutions 0 and 1 are addressed, not refined."""
    assert c.REFINEMENT_ALPHABET[resolution] is None
    with pytest.raises(ResolutionError):
        c.refinement_alphabet(resolution)


@pytest.mark.parametrize("resolution", [-1, 14, 100])
def test_refinement_alphabet_rejects_out_of_range(resolution: int) -> None:
    with pytest.raises(ResolutionError):
        c.refinement_alphabet(resolution)


@pytest.mark.parametrize("resolution", [3.0, "3", None, True])
def test_refinement_alphabet_rejects_non_int(resolution: object) -> None:
    """Including ``True``, which is an int and would silently mean 1."""
    with pytest.raises(ResolutionError):
        c.refinement_alphabet(resolution)  # type: ignore[arg-type]


# --------------------------------------------------------------------------
# WGS84 -- derived at import time and audited against the published values
# --------------------------------------------------------------------------


def test_wgs84_defining_parameters() -> None:
    assert c.WGS84_A == 6378137.0
    assert c.WGS84_INV_F == 298.257223563
    assert c.WGS84_F == 1.0 / c.WGS84_INV_F


def test_wgs84_derived_parameters_match_published_values() -> None:
    """Derived at import time, audited here against the published figures."""
    assert c.WGS84_B == pytest.approx(6356752.314245, abs=1e-6)
    assert c.WGS84_E2 == pytest.approx(6.694379990141e-03, rel=1e-12)
    assert c.WGS84_E == pytest.approx(8.1819190842622e-02, rel=1e-12)
    assert c.WGS84_E == pytest.approx(math.sqrt(c.WGS84_E2), rel=1e-15)
    assert c.WGS84_EP2 == pytest.approx(c.WGS84_E2 / (1.0 - c.WGS84_E2), rel=1e-15)


def test_meridian_quadrant_is_bounded_by_the_two_axes() -> None:
    """Sanity bracket: b < quadrant arc / (pi/2) < a, and closer to a."""
    mean_radius = c.MERIDIAN_QUADRANT / (math.pi / 2)
    assert c.WGS84_B < mean_radius < c.WGS84_A
    assert c.MERIDIAN_QUADRANT == pytest.approx(10001965.729, abs=1e-3)


def test_equator_quadrant_is_pi_times_a() -> None:
    assert c.EQUATOR_QUADRANT == pytest.approx(20037508.3427892, abs=1e-6)


def test_res1_cell_counts_follow_from_the_projection() -> None:
    """2003 along the equator and 1000 along the meridian, both partial.

    The fractional part is what produces the partial last column and row
    dealt with in F4; it is recorded here so the intent is not lost.
    """
    along_equator = c.EQUATOR_QUADRANT / 10_000.0
    along_meridian = c.MERIDIAN_QUADRANT / 10_000.0
    assert math.floor(along_equator) == c.RES1_CELLS_X == 2003
    assert math.floor(along_meridian) == c.RES1_CELLS_Y == 1000
    assert along_equator % 1 == pytest.approx(0.7508, abs=1e-4)
    assert along_meridian % 1 == pytest.approx(0.1966, abs=1e-4)


def test_res1_full_cell_count_reconciles_with_the_max_index() -> None:
    """N full cells occupy indices 0..N-1, so the partial one is index N.

    That is what makes the paper's "approximately 2,003 cells" (section 3.1)
    and Table 1's ``2003/1000`` upper bound the same statement rather than
    an off-by-one between them.
    """
    max_x, max_y = c.RES1_MAX_INDEX.split(c.RES1_SEPARATOR)
    assert int(max_x) == c.RES1_CELLS_X
    assert int(max_y) == c.RES1_CELLS_Y
    # the last full cell sits one index below the partial one
    assert int(max_x) - 1 == c.RES1_CELLS_X - 1
    assert c.EQUATOR_QUADRANT / 10_000.0 > c.RES1_CELLS_X
    assert c.MERIDIAN_QUADRANT / 10_000.0 > c.RES1_CELLS_Y


def test_res1_index_bounds_match_table1() -> None:
    assert c.RES1_MIN_INDEX == "0000/0000"
    assert c.RES1_MAX_INDEX == "2003/1000"
    assert c.RES1_DIGITS == 4
    for bound in (c.RES1_MIN_INDEX, c.RES1_MAX_INDEX):
        x, y = bound.split(c.RES1_SEPARATOR)
        assert len(x) == len(y) == c.RES1_DIGITS


# --------------------------------------------------------------------------
# Quadrantes, sintaxe do indice e formas de celula
# --------------------------------------------------------------------------


def test_quadrants() -> None:
    assert set(c.QUADRANTS) == {"NE", "NW", "SE", "SW"}
    assert len(c.QUADRANTS) == 4
    assert all(len(q) == c.QUADRANT_CODE_LENGTH for q in c.QUADRANTS)


def test_index_syntax_atoms() -> None:
    assert (c.DESCENT_OPEN, c.DESCENT_CLOSE) == ("(", ")")
    assert c.SIBLING_SEPARATOR == ","
    assert c.RES1_SEPARATOR == "/"


def test_canonical_index_example_is_well_formed() -> None:
    """``SE(1400/0374(3(C2(3))))`` parses structurally against the atoms."""
    example = c.INDEX_EXAMPLE_ATOMIC
    assert example.startswith("SE")
    assert example.count(c.DESCENT_OPEN) == example.count(c.DESCENT_CLOSE)
    assert "-" not in example  # no negative index anywhere
    # resolutions 1, 2, 3, 4: base cell, quaternary, quinary, quaternary
    assert "1400/0374" in example
    assert "C2" in c.QUINARY_CODES
    assert "3" in c.QUATERNARY_CODES


def test_cell_shapes() -> None:
    assert c.CELL_SHAPES == ("parallelogram", "triangle", "trapezoid")
    assert c.PARALLELOGRAM_BASE_ANGLE_DEG == 45.0
    assert c.TRIANGLE_BASE_TO_HEIGHT_RATIO == 2.0


def test_triangle_preserves_the_parallelogram_area() -> None:
    """Isosceles cell with base = 2h has area h^2, same as the parallelogram."""
    height = 1.0
    base = c.TRIANGLE_BASE_TO_HEIGHT_RATIO * height
    assert 0.5 * base * height == pytest.approx(height * height)


# --------------------------------------------------------------------------
# Extension zones (paper, Figure 5)
# --------------------------------------------------------------------------


def test_extension_zone_keys_are_the_load_bearing_strings() -> None:
    assert set(c.EXTENSION_ZONES) == {"FIJI", "CHUKOTKA"}


def test_fiji_bounds_match_figure_5() -> None:
    zone = c.EXTENSION_ZONES["FIJI"]
    assert zone.quadrant == "SE"
    assert zone.lon_limit == -178.0
    assert (zone.lat_min, zone.lat_max) == (-21.5, -15.5)


def test_chukotka_bounds_match_figure_5() -> None:
    zone = c.EXTENSION_ZONES["CHUKOTKA"]
    assert zone.quadrant == "NE"
    assert zone.lon_limit == -169.5
    assert (zone.lat_min, zone.lat_max) == (64.0, 72.0)


def test_extension_zones_are_self_consistent() -> None:
    for name, zone in c.EXTENSION_ZONES.items():
        assert zone.name == name
        assert zone.lat_min < zone.lat_max
        assert -180.0 < zone.lon_limit < 0.0  # western longitudes
        assert zone.lon_limit % c.EXTENSION_LON_PRECISION == 0.0
        assert zone.quadrant in c.QUADRANTS
        # eastern quadrants only: the extension pushes east across 180
        assert zone.quadrant.endswith("E")
        hemisphere = "N" if zone.lat_min >= 0 else "S"
        assert zone.quadrant.startswith(hemisphere)


def test_antarctica_is_not_an_extension_zone() -> None:
    """Deliberate exclusion, not an omission to be fixed later."""
    for zone in c.EXTENSION_ZONES.values():
        assert zone.lat_min > -60.0


def _sinusoidal_x(lon_deg: float, lat_deg: float) -> float:
    """Eq. (1) of the paper, inlined so F0 does not depend on F1.

    ``x = lambda * a cos(phi) / sqrt(1 - e^2 sin^2 phi)``. Once F1 lands,
    this should be replaced by ``geodetic_to_sinusoidal``.
    """
    lam = math.radians(lon_deg)
    phi = math.radians(lat_deg)
    return (
        lam * c.WGS84_A * math.cos(phi) / math.sqrt(1 - c.WGS84_E2 * math.sin(phi) ** 2)
    )


def test_extension_zones_stay_within_the_res1_x_bound() -> None:
    """No extension may push X past ``RES1_MAX_INDEX``.

    The zones reach east of 180 deg, which raises x, but they sit at
    latitudes where ``cos(phi)`` lowers it more. Fiji adds 2 deg of
    longitude (+1.11%) and loses 3.6% to the cosine, topping out near
    X = 1953; Chukotka tops out near X = 932. The bound holds — but only
    because of where the zones are. A zone declared near the equator with
    the same reach would break it, and this test is what would say so.
    """
    max_x = int(c.RES1_MAX_INDEX.split(c.RES1_SEPARATOR)[0])
    for name, zone in c.EXTENSION_ZONES.items():
        # eastward extent, in degrees from the prime meridian
        lon_extent = 360.0 + zone.lon_limit
        # x is largest at whichever bound lies closer to the equator
        nearest_lat = min(abs(zone.lat_min), abs(zone.lat_max))
        reached = _sinusoidal_x(lon_extent, nearest_lat) / 10_000.0
        assert reached < max_x, f"{name} reaches X = {reached:.2f}"


def test_the_x_bound_is_attained_at_the_equator() -> None:
    """2003.75 is a global maximum of x, reached on the antemeridian at phi = 0."""
    at_equator = _sinusoidal_x(180.0, 0.0)
    assert at_equator == pytest.approx(c.EQUATOR_QUADRANT, rel=1e-12)
    for latitude in (5.0, 15.5, 30.0, 45.0, 64.0, 80.0):
        assert _sinusoidal_x(180.0, latitude) < at_equator


# --------------------------------------------------------------------------
# Imutabilidade
# --------------------------------------------------------------------------


def test_tables_are_immutable() -> None:
    for table in (
        c.CELL_SIZE_M,
        c.CELL_AREA_M2,
        c.VISUALIZATION_SCALE,
        c.ANALYSIS_SCALE,
        c.QUADRANTS,
        c.QUATERNARY_CODES,
        c.QUINARY_CODES,
        c.REFINEMENT_RATIO,
        c.TOKENIZABLE_RESOLUTIONS,
    ):
        assert isinstance(table, tuple)


def test_extension_zones_mapping_is_read_only() -> None:
    zones = c.EXTENSION_ZONES
    with pytest.raises(TypeError):
        zones["ANTARCTICA"] = zones["FIJI"]  # type: ignore[index]


# --------------------------------------------------------------------------
# Tokenizable resolutions
# --------------------------------------------------------------------------


def test_tokenizable_resolutions_have_power_of_ten_sides() -> None:
    for resolution in c.TOKENIZABLE_RESOLUTIONS:
        side = c.CELL_SIZE_M[resolution]
        assert side is not None
        exponent = math.log10(side)
        assert exponent == pytest.approx(round(exponent), abs=1e-12)


def test_non_tokenizable_resolutions_do_not() -> None:
    others = [r for r in range(1, 14) if r not in c.TOKENIZABLE_RESOLUTIONS]
    assert others == [2, 4, 6, 8, 10, 12]
    for resolution in others:
        side = c.CELL_SIZE_M[resolution]
        assert side is not None
        exponent = math.log10(side)
        assert exponent != pytest.approx(round(exponent), abs=1e-9)


# --------------------------------------------------------------------------
# Contrato de API - os nomes que __init__.py importa precisam existir
# --------------------------------------------------------------------------


def _names_imported_from(module_name: str) -> list[str]:
    init = pathlib.Path(c.__file__).with_name("__init__.py")
    tree = ast.parse(init.read_text(encoding="utf-8"))
    names: list[str] = []
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.ImportFrom)
            and node.level == 1
            and node.module == module_name
        ):
            names.extend(alias.name for alias in node.names)
    return names


def test_init_imports_from_constants_all_exist() -> None:
    """A constant renamed here breaks ``import itacart`` at package level."""
    expected = _names_imported_from("constants")
    assert expected, "no names imported from .constants - parse failed"
    missing = [name for name in expected if not hasattr(c, name)]
    assert missing == []


def test_init_imports_from_exceptions_all_exist() -> None:
    from itacart import exceptions

    expected = _names_imported_from("exceptions")
    assert expected, "no names imported from .exceptions - parse failed"
    missing = [name for name in expected if not hasattr(exceptions, name)]
    assert missing == []
