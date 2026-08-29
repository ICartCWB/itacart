"""Named geographic points shared by the test suite and the notebooks.

Plain data, no ``pytest`` import: a notebook cannot consume a fixture, and
the F1 brief requires the notebook and the tests to exercise **the same**
points. A fixture that diverges from the notebook is how an error passes
the test and shows up in the plot, or the other way round.

``tests/conftest.py`` should re-export these through its geographic
fixtures rather than restating the coordinates.

Every entry is ``(lon_deg, lat_deg)`` -- the core argument order of the
package. The dictionary keys are stable; boundary fixtures for the
extension zones are expected to build on ``suva`` and ``wrangel``.
"""

from __future__ import annotations

from typing import Final

REFERENCE_POINTS: Final[dict[str, tuple[float, float]]] = {
    # Origin and axes of the projection
    "origin": (0.0, 0.0),
    "greenwich_observatory": (-0.0015, 51.4779),
    "antemeridian_equator": (180.0, 0.0),
    "north_pole": (0.0, 90.0),
    "south_pole": (0.0, -90.0),
    # Project sites
    "ita_sjc": (-45.9009, -23.2237),
    "curitiba": (-49.2733, -25.4284),
    # Figure 7 of the paper (Central Park, resolutions 6 and 7)
    "central_park": (-73.9665, 40.7812),
    # Figure 8: the three angular-distortion regimes
    "italian_peninsula": (13.0, 42.0),
    "east_china_sea": (125.0, 30.0),
    "kamchatka": (159.0, 56.0),
    # Section 3.2: the antemeridian extension zones
    "suva": (178.4419, -18.1416),
    "wrangel": (-179.4000, 71.2333),
    # High latitude, far from the origin of the projection
    "longyearbyen": (15.6469, 78.2232),
    # Section 3.2: positions that must NOT fall in an extension zone
    "central_pacific": (-150.0, -10.0),
    "bering_sea_west_of_the_limit": (-150.0, 68.0),
    # Section 3.2: the prime-meridian triangles
    "prime_meridian_north": (0.0, 42.0),
    "just_west_of_greenwich": (-0.02, 51.4779),
}
"""Named ``(lon_deg, lat_deg)`` positions, in degrees."""

ROUNDTRIP_POINTS: Final[tuple[str, ...]] = (
    "origin",
    "greenwich_observatory",
    "ita_sjc",
    "curitiba",
    "central_park",
    "italian_peninsula",
    "east_china_sea",
    "kamchatka",
    "suva",
    "wrangel",
    "longyearbyen",
)
"""Points suitable for projection round-trip checks.

Excludes the poles, where the parallel circle collapses and longitude is
not recoverable, and the antemeridian point, whose longitude round-trips
to the same value only up to the sign convention chosen in F4.
"""

NEAR_ANTIPODAL_PAIRS: Final[tuple[tuple[str, tuple[float, float]], ...]] = (
    ("separation 175 deg", (175.0, 0.0)),
    ("separation 179 deg", (179.0, 0.0)),
    ("separation 179.5 deg", (179.5, 0.0)),
    ("separation 179.9 deg", (179.9, 0.0)),
    ("exactly antipodal", (180.0, 0.0)),
)
"""Second endpoints of pairs starting at ``origin``, of growing separation.

Used to locate the frontier where the Vincenty inverse stops converging
(acceptance criterion 5). The frontier itself is measured, not asserted:
what the test pins is that beyond it the failure is a
``ConvergenceError`` and not a plausible number.
"""
