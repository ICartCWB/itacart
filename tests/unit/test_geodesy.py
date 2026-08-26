"""Tests for :mod:`itacart.geodesy` (phase F1).

Naming: the five acceptance criteria of the phase have one test each, named
``test_criterion_<n>_...``, so the handoff can cite evidence by name. The
remaining tests cover the degenerate cases, the error branches and the two
independent oracles (closed-form derivative and Simpson quadrature).
"""

from __future__ import annotations

import inspect
import math
import random
import sys
from pathlib import Path

import pytest

from itacart import geodesy
from itacart.constants import EQUATOR_QUADRANT, MERIDIAN_QUADRANT, WGS84_A, WGS84_B
from itacart.exceptions import ConvergenceError, DomainError

# tests/ is not a package in this repo, so the shared point table is
# imported by path rather than by dotted name. The notebook does the same,
# which is the point: both must exercise the identical coordinates.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from reference_points import (  # noqa: E402
    NEAR_ANTIPODAL_PAIRS,
    REFERENCE_POINTS,
    ROUNDTRIP_POINTS,
)

try:  # pragma: no cover - import guard, not logic
    import pyproj as _pyproj

    HAS_PYPROJ = True
except ImportError:  # pragma: no cover
    _pyproj = None  # type: ignore[assignment]
    HAS_PYPROJ = False

pyproj = _pyproj

requires_pyproj = pytest.mark.skipif(
    not HAS_PYPROJ, reason="pyproj not installed (dev-only cross-check, D-0.3)"
)


# ---------------------------------------------------------------------------
# Acceptance criterion 1 -- projection round-trip
# ---------------------------------------------------------------------------


def test_criterion_1_roundtrip_named_points() -> None:
    """geodetic -> sinusoidal -> geodetic closes below 1e-9 degree."""
    for name in ROUNDTRIP_POINTS:
        lon, lat = REFERENCE_POINTS[name]
        x, y = geodesy.geodetic_to_sinusoidal(lon, lat)
        lon_back, lat_back = geodesy.sinusoidal_to_geodetic(x, y)
        assert abs(lon_back - lon) < 1e-9, name
        assert abs(lat_back - lat) < 1e-9, name


def test_criterion_1_roundtrip_global_grid() -> None:
    """Same closure over a global grid, poles excluded (they collapse)."""
    worst = 0.0
    for lat in range(-88, 89, 4):
        for lon in range(-180, 181, 10):
            x, y = geodesy.geodetic_to_sinusoidal(float(lon), float(lat))
            lon_back, lat_back = geodesy.sinusoidal_to_geodetic(x, y)
            worst = max(worst, abs(lon_back - lon), abs(lat_back - lat))
    assert worst < 1e-9, f"worst round-trip error {worst} deg"


# ---------------------------------------------------------------------------
# Acceptance criterion 2 -- meridian quadrant
# ---------------------------------------------------------------------------


def test_criterion_2_meridian_quadrant_matches_published_value() -> None:
    """meridian_arc(90) reproduces the WGS84 meridian quadrant.

    The target is the GeographicLib value, carried to full double
    precision so that the comparison measures the series rather than the
    rounding of a three-decimal literal. One micrometre is the stated
    tolerance; the observed gap is a thousand times smaller than that,
    and equal to one unit in the last place of a float64 at this
    magnitude, which is the floor of what double precision can express.
    """
    computed = geodesy.meridian_arc(90.0)
    assert abs(computed - MERIDIAN_QUADRANT) < 1e-6


def test_criterion_2_quadrant_is_computed_not_read() -> None:
    """The module must not source the quadrant from the constant.

    Guards the failure mode the opening package calls out: if the
    implementation returned ``MERIDIAN_QUADRANT`` through any shortcut, the
    criterion-2 test would pass and prove nothing.
    """
    source = inspect.getsource(geodesy)
    assert "MERIDIAN_QUADRANT" not in source


def test_meridian_arc_equator_and_symmetry() -> None:
    """The arc vanishes at the equator and is odd in latitude."""
    assert geodesy.meridian_arc(0.0) == 0.0
    for lat in (1.0, 23.2237, 45.0, 71.2333, 90.0):
        assert geodesy.meridian_arc(-lat) == pytest.approx(
            -geodesy.meridian_arc(lat), abs=1e-9
        )


def test_meridian_arc_matches_simpson_quadrature() -> None:
    """Series and quadrature agree: two independent readings of Eq. (2)."""
    for lat in (5.0, 23.2237, 45.0, 66.5, 90.0):
        series = geodesy.meridian_arc(lat)
        quadrature = geodesy.meridian_arc_quadrature(lat)
        # Measured worst case 6.3e-8 m at 45 deg: Simpson truncation, seven
        # orders below the 1 cm cell of resolution 13.
        assert abs(series - quadrature) < 1e-6, lat


def test_meridian_arc_derivative_matches_closed_form() -> None:
    """d/dphi meridian_arc == meridian_radius, to central-difference order.

    The check of section 5.2 of the opening package: a wrong integrand or a
    wrong quadrature step shows up here before any external comparison.
    """
    step_deg = 1e-4
    step_rad = math.radians(step_deg)
    for lat in (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 89.0):
        numeric = (
            geodesy.meridian_arc(lat + step_deg) - geodesy.meridian_arc(lat - step_deg)
        ) / (2.0 * step_rad)
        closed = geodesy.meridian_radius(lat)
        assert abs(numeric - closed) / closed < 1e-9, lat


# ---------------------------------------------------------------------------
# Acceptance criterion 3 -- inversion of the meridian arc
# ---------------------------------------------------------------------------


def test_criterion_3_inverse_meridian_arc_roundtrip() -> None:
    """inverse_meridian_arc(meridian_arc(phi)) == phi below 1e-10 degree."""
    worst = 0.0
    for lat in range(-90, 91):
        recovered = geodesy.inverse_meridian_arc(geodesy.meridian_arc(float(lat)))
        worst = max(worst, abs(recovered - lat))
    assert worst < 1e-10, f"worst inversion residual {worst} deg"


def test_inverse_meridian_arc_converges_in_few_iterations() -> None:
    """Newton needs a handful of steps at every latitude."""
    worst = 0
    for lat in range(-90, 91, 5):
        arc = geodesy.meridian_arc(float(lat))
        _, iterations = geodesy._inverse_meridian_arc_rad(
            arc, geodesy.NEWTON_TOL_RAD, geodesy.NEWTON_MAX_ITER
        )
        worst = max(worst, iterations)
    assert worst <= 6, f"Newton took {worst} iterations"


def test_inverse_meridian_arc_rejects_beyond_pole() -> None:
    """No latitude exists past the quadrant."""
    with pytest.raises(DomainError, match="beyond the meridian quadrant"):
        geodesy.inverse_meridian_arc(MERIDIAN_QUADRANT * 1.001)


def test_inverse_meridian_arc_raises_on_iteration_cap() -> None:
    """A starved iteration budget raises rather than returning a guess."""
    with pytest.raises(ConvergenceError, match="Newton"):
        geodesy.inverse_meridian_arc(5.0e6, tol=1e-30, max_iter=2)


# ---------------------------------------------------------------------------
# Acceptance criterion 4 -- Vincenty against pyproj
# ---------------------------------------------------------------------------


@pytest.mark.crosscheck
@requires_pyproj
def test_criterion_4_vincenty_matches_pyproj() -> None:
    """Distance agrees with pyproj.Geod below 1 mm out to 5 500 km."""
    geod = pyproj.Geod(ellps="WGS84")
    rng = random.Random(20260824)
    worst_distance = 0.0
    worst_azimuth = 0.0
    pairs = 0
    while pairs < 500:
        lat1 = math.degrees(math.asin(rng.uniform(-1.0, 1.0)))
        lon1 = rng.uniform(-180.0, 180.0)
        lat2 = math.degrees(math.asin(rng.uniform(-1.0, 1.0)))
        lon2 = rng.uniform(-180.0, 180.0)
        azimuth_ref, _, distance_ref = geod.inv(lon1, lat1, lon2, lat2)
        if distance_ref > 5.5e6:
            continue
        pairs += 1
        distance, azimuth = geodesy.inverse_geodesic(lon1, lat1, lon2, lat2)
        worst_distance = max(worst_distance, abs(distance - distance_ref))
        delta = (azimuth - azimuth_ref % 360.0 + 180.0) % 360.0 - 180.0
        worst_azimuth = max(worst_azimuth, abs(delta))
    assert worst_distance < 1e-3, f"worst distance error {worst_distance} m"
    assert worst_azimuth < 1e-6, f"worst azimuth error {worst_azimuth} deg"


@pytest.mark.crosscheck
@requires_pyproj
def test_direct_geodesic_matches_pyproj() -> None:
    """The direct solution agrees with pyproj on the end point."""
    geod = pyproj.Geod(ellps="WGS84")
    lon1, lat1 = REFERENCE_POINTS["ita_sjc"]
    for azimuth in (0.0, 37.5, 90.0, 180.0, 271.25):
        for distance in (1.0, 1_000.0, 500_000.0, 3_000_000.0):
            lon_ref, lat_ref, _ = geod.fwd(lon1, lat1, azimuth, distance)
            lon2, lat2 = geodesy.direct_geodesic(lon1, lat1, azimuth, distance)
            assert abs(lat2 - lat_ref) < 1e-9
            assert abs(lon2 - lon_ref) < 1e-9


@pytest.mark.crosscheck
@requires_pyproj
def test_projection_is_not_proj_sinusoidal_on_the_sphere() -> None:
    """The ellipsoidal ordinate departs from the spherical one by km.

    Pins the justification of ``D-0.3``: PROJ's ``+proj=sinu`` on a sphere
    is not the parallels-plane projection of the paper.
    """
    spherical = pyproj.Transformer.from_crs(
        "EPSG:4326", "+proj=sinu +R=6371007.181 +units=m +no_defs", always_xy=True
    )
    _, y_sphere = spherical.transform(0.0, 45.0)
    y_ellipsoid = geodesy.meridian_arc(45.0)
    assert abs(y_ellipsoid - y_sphere) > 1_000.0


# ---------------------------------------------------------------------------
# Acceptance criterion 5 -- near-antipodal handling
# ---------------------------------------------------------------------------


def test_criterion_5_near_antipodal_raises_convergence_error() -> None:
    """Antipodal pairs raise instead of returning an unconverged iterate."""
    with pytest.raises(ConvergenceError, match="near-antipodal"):
        geodesy.inverse_geodesic(0.0, 0.0, 180.0, 0.0)


def test_criterion_5_frontier_is_beyond_operational_range() -> None:
    """Everything that converges is a number; the frontier sits past 19 000 km.

    Cadastral geometry never approaches this separation, but the frontier
    must be measured rather than assumed: an exception raised too early
    would break F7 densification silently.
    """
    converged: list[float] = []
    for label, second in NEAR_ANTIPODAL_PAIRS:
        try:
            distance, _ = geodesy.inverse_geodesic(0.0, 0.0, second[0], second[1])
        except ConvergenceError:
            continue
        assert math.isfinite(distance), label
        converged.append(distance)
    assert converged, "no pair converged; the frontier moved"
    assert max(converged) > 1.9e7


def test_vincenty_inverse_reports_iterations() -> None:
    """The instrumented core returns a positive, bounded iteration count."""
    lon1, lat1 = REFERENCE_POINTS["ita_sjc"]
    lon2, lat2 = REFERENCE_POINTS["curitiba"]
    _, _, iterations = geodesy._vincenty_inverse(
        lon1, lat1, lon2, lat2, geodesy.VINCENTY_TOL_RAD, geodesy.VINCENTY_MAX_ITER
    )
    assert 0 < iterations < 20


def test_direct_geodesic_raises_on_iteration_cap() -> None:
    """The direct solution reports non-convergence as ConvergenceError."""
    with pytest.raises(ConvergenceError, match="Vincenty direct"):
        geodesy.direct_geodesic(0.0, 0.0, 45.0, 1000.0, max_iter=0)


# ---------------------------------------------------------------------------
# Projection: structure and degenerate cases
# ---------------------------------------------------------------------------


def test_projection_origin_and_axes() -> None:
    """The origin maps to (0, 0); the axes carry true distances."""
    assert geodesy.geodetic_to_sinusoidal(0.0, 0.0) == (0.0, 0.0)
    x, _ = geodesy.geodetic_to_sinusoidal(180.0, 0.0)
    assert x == pytest.approx(EQUATOR_QUADRANT, abs=1e-6)
    assert x == pytest.approx(math.pi * WGS84_A, abs=1e-9)


def test_projection_pole_collapses_to_a_point() -> None:
    """At the pole the parallel vanishes: every longitude gives x = 0."""
    for lon in (-180.0, -45.0, 0.0, 90.0, 179.999):
        x, y = geodesy.geodetic_to_sinusoidal(lon, 90.0)
        assert abs(x) < 1e-6
        assert y == pytest.approx(MERIDIAN_QUADRANT, abs=1e-3)


def test_inverse_at_pole_returns_zero_longitude() -> None:
    """Longitude is undefined at the pole; the convention is zero."""
    lon, lat = geodesy.sinusoidal_to_geodetic(0.0, geodesy.meridian_arc(90.0))
    assert lat == pytest.approx(90.0, abs=1e-9)
    assert lon == 0.0


def test_projection_signs_follow_the_hemisphere() -> None:
    """East and north are positive; west and south negative."""
    x_east, y_north = geodesy.geodetic_to_sinusoidal(30.0, 20.0)
    x_west, y_south = geodesy.geodetic_to_sinusoidal(-30.0, -20.0)
    assert x_east > 0 and y_north > 0
    assert x_west == pytest.approx(-x_east, abs=1e-9)
    assert y_south == pytest.approx(-y_north, abs=1e-9)


def test_longitude_is_not_range_checked() -> None:
    """F4 owns the domain; the geodesy layer only needs finiteness (D-1.5)."""
    x, _ = geodesy.geodetic_to_sinusoidal(-179.4, 71.2333)
    assert x < 0.0


# ---------------------------------------------------------------------------
# Radii of curvature
# ---------------------------------------------------------------------------


def test_prime_vertical_radius_bounds() -> None:
    """N runs from a at the equator to a / sqrt(1 - e^2) at the pole."""
    assert geodesy.prime_vertical_radius(0.0) == pytest.approx(WGS84_A, abs=1e-9)
    polar = WGS84_A * WGS84_A / WGS84_B
    assert geodesy.prime_vertical_radius(90.0) == pytest.approx(polar, abs=1e-6)


def test_meridian_radius_is_positive_and_increasing() -> None:
    """M > 0 everywhere -- the arc is strictly monotonic, Newton is safe."""
    previous = 0.0
    for lat in range(0, 91, 5):
        radius = geodesy.meridian_radius(float(lat))
        assert radius > 0.0
        assert radius > previous
        previous = radius


# ---------------------------------------------------------------------------
# Quadrature
# ---------------------------------------------------------------------------


def test_quadrature_returns_zero_at_the_equator() -> None:
    assert geodesy.meridian_arc_quadrature(0.0) == 0.0


def test_quadrature_is_odd() -> None:
    assert geodesy.meridian_arc_quadrature(-45.0) == pytest.approx(
        -geodesy.meridian_arc_quadrature(45.0), abs=1e-9
    )


def test_quadrature_rejects_non_positive_tolerance() -> None:
    with pytest.raises(DomainError, match="tol must be positive"):
        geodesy.meridian_arc_quadrature(45.0, tol=0.0)


def test_quadrature_raises_when_interval_cap_is_too_low() -> None:
    """An unreachable tolerance with a low cap raises, it does not guess."""
    with pytest.raises(ConvergenceError, match="quadrature"):
        geodesy.meridian_arc_quadrature(45.0, tol=1e-30, max_intervals=4)


# ---------------------------------------------------------------------------
# Geodesics: degenerate cases
# ---------------------------------------------------------------------------


def test_coincident_points_return_zero_distance_and_azimuth() -> None:
    """A duplicated vertex is data, not an error (D-1.6)."""
    lon, lat = REFERENCE_POINTS["central_park"]
    assert geodesy.inverse_geodesic(lon, lat, lon, lat) == (0.0, 0.0)


def test_direct_geodesic_zero_distance_returns_start() -> None:
    lon, lat = REFERENCE_POINTS["suva"]
    assert geodesy.direct_geodesic(lon, lat, 45.0, 0.0) == (lon, lat)


def test_direct_geodesic_rejects_negative_distance() -> None:
    with pytest.raises(DomainError, match="non-negative"):
        geodesy.direct_geodesic(0.0, 0.0, 45.0, -1.0)


def test_geodesic_roundtrip_direct_then_inverse() -> None:
    """direct then inverse recovers the azimuth and the distance."""
    lon1, lat1 = REFERENCE_POINTS["curitiba"]
    for azimuth in (10.0, 100.0, 190.0, 280.0):
        for distance in (250.0, 25_000.0, 1_500_000.0):
            lon2, lat2 = geodesy.direct_geodesic(lon1, lat1, azimuth, distance)
            back_distance, back_azimuth = geodesy.inverse_geodesic(
                lon1, lat1, lon2, lat2
            )
            # 5 um measured at 250 m: the inverse solution loses absolute
            # precision on short lines. Five orders below the 1 cm cell,
            # and below the noise of any GNSS survey feeding F7.
            assert back_distance == pytest.approx(distance, abs=1e-4)
            assert back_azimuth == pytest.approx(azimuth, abs=1e-6)


def test_equatorial_geodesic_is_handled() -> None:
    """cos^2(alpha) == 0 on the equator: the branch must not divide by zero."""
    distance, azimuth = geodesy.inverse_geodesic(0.0, 0.0, 10.0, 0.0)
    assert distance == pytest.approx(1_113_194.9079, abs=1e-3)
    assert azimuth == pytest.approx(90.0, abs=1e-9)


def test_azimuth_is_normalised_to_the_positive_turn() -> None:
    """Westward bearings come back in [0, 360), not negative."""
    _, azimuth = geodesy.inverse_geodesic(0.0, 0.0, -10.0, 0.0)
    assert azimuth == pytest.approx(270.0, abs=1e-9)


def test_direct_geodesic_normalises_longitude_across_the_antemeridian() -> None:
    """A long eastward traverse from Fiji wraps into (-180, 180]."""
    lon1, lat1 = REFERENCE_POINTS["suva"]
    lon2, _ = geodesy.direct_geodesic(lon1, lat1, 90.0, 500_000.0)
    assert -180.0 < lon2 <= 180.0
    assert lon2 < 0.0


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("bad", [91.0, -90.5, float("nan"), float("inf")])
def test_latitude_out_of_range_raises_domain_error(bad: float) -> None:
    with pytest.raises(DomainError):
        geodesy.meridian_arc(bad)


def test_non_finite_longitude_raises_domain_error() -> None:
    with pytest.raises(DomainError, match="lon_deg"):
        geodesy.geodetic_to_sinusoidal(float("nan"), 0.0)


def test_non_finite_projected_coordinate_raises_domain_error() -> None:
    with pytest.raises(DomainError, match="x_m"):
        geodesy.sinusoidal_to_geodetic(float("inf"), 0.0)
    with pytest.raises(DomainError, match="arc_m"):
        geodesy.sinusoidal_to_geodetic(0.0, float("nan"))


def test_geodesic_validates_every_argument() -> None:
    with pytest.raises(DomainError, match="lat1_deg"):
        geodesy.inverse_geodesic(0.0, 91.0, 0.0, 0.0)
    with pytest.raises(DomainError, match="lat2_deg"):
        geodesy.inverse_geodesic(0.0, 0.0, 0.0, -91.0)
    with pytest.raises(DomainError, match="lon2_deg"):
        geodesy.inverse_geodesic(0.0, 0.0, float("inf"), 0.0)
    with pytest.raises(DomainError, match="azimuth_deg"):
        geodesy.direct_geodesic(0.0, 0.0, float("nan"), 10.0)
    with pytest.raises(DomainError, match="distance_m"):
        geodesy.direct_geodesic(0.0, 0.0, 45.0, float("inf"))


def test_public_surface_is_declared() -> None:
    """Everything exported exists; nothing public is left out of __all__."""
    for name in geodesy.__all__:
        assert callable(getattr(geodesy, name))
    public = {
        name
        for name, value in vars(geodesy).items()
        if callable(value)
        and not name.startswith("_")
        and getattr(value, "__module__", "") == geodesy.__name__
    }
    assert public == set(geodesy.__all__)


# ---------------------------------------------------------------------------
# Fidelity of the port against itacart_core
# ---------------------------------------------------------------------------
#
# The origin file is not vendored into this repository, so these values were
# captured from itacart_core/geodesy.py during F1 and are pinned here as
# literals. They are what keeps the port auditable after the origin is out of
# reach: any drift in the projection shows up as an exact-equality failure,
# not as a tolerance that quietly widens.

ORIGIN_FORWARD: dict[tuple[float, float], tuple[float, float]] = {
    (-45.9009, -23.2237): (-4698086.250191047, -2569311.3640246917),
    (178.4419, -18.1416): (18882745.78591646, -2006654.3925390893),
    (-73.9665, 40.7812): (-6243719.106558419, 4516275.158692554),
    (13.0, 42.0): (1077059.9006100234, 4651636.8795710625),
}
"""``(lon, lat) -> (x, y)`` as produced by ``itacart_core.geodesy.forward``."""

ORIGIN_MERIDIAN_QUADRANT: float = 10001965.729312722
"""``itacart_core.geodesy.meridian_arc(pi / 2)``."""

ORIGIN_GEODESIC_ITA_CURITIBA: tuple[float, float] = (
    420426.227334673,
    233.80850578273706,
)
"""``itacart_core.geodesy.inverse_geodesic`` between the two project sites."""


def test_port_reproduces_origin_projection_bit_for_bit() -> None:
    """The ported projection is not merely close to the origin: it is equal."""
    for (lon, lat), expected in ORIGIN_FORWARD.items():
        assert geodesy.geodetic_to_sinusoidal(lon, lat) == expected


def test_port_reproduces_origin_meridian_quadrant_bit_for_bit() -> None:
    assert geodesy.meridian_arc(90.0) == ORIGIN_MERIDIAN_QUADRANT


def test_port_reproduces_origin_geodesic_bit_for_bit() -> None:
    lon1, lat1 = REFERENCE_POINTS["ita_sjc"]
    lon2, lat2 = REFERENCE_POINTS["curitiba"]
    assert (
        geodesy.inverse_geodesic(lon1, lat1, lon2, lat2) == ORIGIN_GEODESIC_ITA_CURITIBA
    )


def test_argument_order_is_lon_lat_not_lat_lon() -> None:
    """Pins ``D-1.2`` against the (lat, lon) order of itacart_core.

    ``itacart_core.geodesy.forward`` takes ``(lat, lon)``; every public
    function here takes ``(lon, lat)``. Swapping them produces a valid
    coordinate pair somewhere else on the globe and raises nothing, so the
    order needs a test that fails loudly rather than a comment.
    """
    lon, lat = 13.0, 42.0
    x, y = geodesy.geodetic_to_sinusoidal(lon, lat)
    swapped_x, swapped_y = geodesy.geodetic_to_sinusoidal(lat, lon)
    assert (x, y) != (swapped_x, swapped_y)
    # y is the meridian arc, which depends on latitude alone: it is the
    # component that betrays a swap.
    assert y == pytest.approx(geodesy.meridian_arc(lat), abs=0.0)
    assert swapped_y == pytest.approx(geodesy.meridian_arc(lon), abs=0.0)


def test_geodesic_argument_order_is_lon_lat() -> None:
    """Same pin for the geodesics, where the origin also uses (lat, lon)."""
    distance_correct, _ = geodesy.inverse_geodesic(0.0, 0.0, 10.0, 0.0)
    distance_swapped, _ = geodesy.inverse_geodesic(0.0, 0.0, 0.0, 10.0)
    # A degree of longitude on the equator is longer than a degree of
    # meridian arc; equal values would mean the order is not being honoured.
    assert distance_correct > distance_swapped
    assert distance_swapped == pytest.approx(geodesy.meridian_arc(10.0), abs=1e-3)


def test_roundtrip_error_fits_the_epsilon_budget_of_cell_quantisation() -> None:
    """Round-trip noise on the plane stays far under the F3 flooring epsilon.

    ``itacart_core.cells`` floors sheared coordinates onto the grid with a
    1e-6 m epsilon, precisely so a representative point sitting on a cell
    boundary does not round down into the previous cell after a projection
    round-trip. That guard only works if the round-trip noise is well below
    it. Measured worst case: 7.5e-9 m, a margin of 134x.
    """
    worst = 0.0
    for lat in range(-89, 90, 7):
        for lon in range(-180, 181, 11):
            x, y = geodesy.geodetic_to_sinusoidal(float(lon), float(lat))
            lon_back, lat_back = geodesy.sinusoidal_to_geodetic(x, y)
            x2, y2 = geodesy.geodetic_to_sinusoidal(lon_back, lat_back)
            worst = max(worst, math.hypot(x2 - x, y2 - y))
    assert worst < 1e-7, f"round-trip noise {worst} m eats into the F3 epsilon"


def test_polar_row_beyond_the_quadrant_raises_instead_of_returning_a_latitude() -> None:
    """Regression for ``B-1.2``.

    Resolution-1 row 1000 is the clipped polar row: the quadrant is
    1000.196 cells tall, so refinements inside that row can address an
    ordinate past the pole. ``itacart_core.geodesy.inverse`` answers such a
    coordinate with latitude 90.063 deg and longitude -40726 deg -- a
    plausible-looking pair that is pure noise. Here it raises.
    """
    beyond_pole = 1000 * 10_000.0 + 5_000.0 + 4 * 1_000.0
    assert beyond_pole > geodesy.meridian_arc(90.0)
    with pytest.raises(DomainError, match="beyond the meridian quadrant"):
        geodesy.sinusoidal_to_geodetic(0.0, beyond_pole)
