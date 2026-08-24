"""Ellipsoidal geodesy: the parallels-plane (sinusoidal) projection and Vincenty.

The forward projection implements Eq. (1) and (2) of the paper directly on
the WGS84 ellipsoid, so PROJ is never called at runtime::

    x = lambda * a * cos(phi) / sqrt(1 - e^2 sin^2 phi)
    y = integral_0^phi a (1 - e^2) / (1 - e^2 sin^2 Phi)^{3/2} dPhi

The meridian arc integral is evaluated by composite Simpson quadrature to a
tolerance of 1e-10 m; the inverse solves for ``phi`` by Newton iteration on
the same integral.

Origem: itacart_core/geodesy.py (F1) + extensao Vincenty (F1+).
Concordancia medida contra ``pyproj.Geod``: 12 um em pares de 5 500 km.
"""

from __future__ import annotations

__all__ = [
    "geodetic_to_sinusoidal",
    "sinusoidal_to_geodetic",
    "meridian_arc",
    "inverse_meridian_arc",
    "inverse_geodesic",
    "direct_geodesic",
]


# --------------------------------------------------------------------------
# Parallels-plane (sinusoidal) projection on the ellipsoid
# --------------------------------------------------------------------------


def geodetic_to_sinusoidal(lon: float, lat: float) -> tuple[float, float]:
    """Project geodetic coordinates onto the ellipsoidal sinusoidal plane.

    Implements Eq. (1) and (2) of the paper. Axis order is ``(lon, lat)``
    in, ``(x, y)`` out, both consistent with ``always_xy=True``.

    Args:
        lon: Longitude in decimal degrees, in ``[-180, 180]``.
        lat: Latitude in decimal degrees, in ``[-90, 90]``.

    Returns:
        ``(x, y)`` in metres on the projection plane.

    Raises:
        DomainError: If either coordinate is outside its valid range.
    """
    raise NotImplementedError


def sinusoidal_to_geodetic(x: float, y: float) -> tuple[float, float]:
    """Invert the ellipsoidal sinusoidal projection.

    Args:
        x: Easting in metres on the projection plane.
        y: Northing in metres on the projection plane.

    Returns:
        ``(lon, lat)`` in decimal degrees.

    Raises:
        DomainError: If the point lies outside the projection envelope.
    """
    raise NotImplementedError


def meridian_arc(lat: float) -> float:
    """Meridian arc length from the equator to ``lat``.

    Evaluates Eq. (2) by composite Simpson quadrature, refining until the
    estimate is stable to 1e-10 m.

    Args:
        lat: Latitude in decimal degrees.

    Returns:
        Arc length in metres, signed with the latitude.
    """
    raise NotImplementedError


def inverse_meridian_arc(distance: float) -> float:
    """Recover the latitude whose meridian arc equals ``distance``.

    Newton iteration on :func:`meridian_arc`; the derivative is the
    meridian radius of curvature, available in closed form.

    Args:
        distance: Arc length in metres from the equator, signed.

    Returns:
        Latitude in decimal degrees.

    Raises:
        ConvergenceError: If the iteration does not converge.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------
# Vincenty (1975) geodesics
# --------------------------------------------------------------------------


def inverse_geodesic(
    lat1: float, lon1: float, lat2: float, lon2: float
) -> tuple[float, float]:
    """Solve the inverse geodesic problem on WGS84.

    Origem: itacart_core/geodesy.py (F1+), ja validado.

    Args:
        lat1: Latitude of the first point, decimal degrees.
        lon1: Longitude of the first point, decimal degrees.
        lat2: Latitude of the second point, decimal degrees.
        lon2: Longitude of the second point, decimal degrees.

    Returns:
        ``(distance_m, initial_azimuth_deg)``.

    Raises:
        ConvergenceError: On nearly antipodal points that fail to converge.
    """
    raise NotImplementedError


def direct_geodesic(
    lat1: float, lon1: float, azimuth_deg: float, distance_m: float
) -> tuple[float, float]:
    """Solve the direct geodesic problem on WGS84.

    Origem: itacart_core/geodesy.py (F1+), ja validado.

    Args:
        lat1: Latitude of the origin, decimal degrees.
        lon1: Longitude of the origin, decimal degrees.
        azimuth_deg: Initial azimuth, decimal degrees clockwise from north.
        distance_m: Distance to travel along the geodesic, metres.

    Returns:
        ``(lat2, lon2)`` in decimal degrees.
    """
    raise NotImplementedError
