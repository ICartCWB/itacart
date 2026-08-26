"""Immutable constants of the ITACaRT specification.

All values transcribed from Silva, Dietzsch & Shiguemori (2025),
*Revista Brasileira de Cartografia*, v. 77. DOI: 10.14393/rbcv77n0a-79281

Origem: itacart_core/resolutions.py + Tabela 1 do artigo.

The module holds data, not behaviour. The single exception is
:func:`refinement_alphabet`, a lookup over ``REFINEMENT_ALPHABET`` that
exists so the even/odd predicate is written once in the whole package
(see the ``B-0.1`` warning on :data:`QUINARY_CODES`).

Nota de leitura: os comentarios de secao seguem a numeracao do artigo.
Toda constante carrega a referencia da figura, tabela ou equacao de
origem, para que a conferencia contra o PDF seja mecanica.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from types import MappingProxyType
from typing import Final, Literal, NamedTuple

from .exceptions import ResolutionError

# --------------------------------------------------------------------------
# WGS84 ellipsoid (artigo, secao 3)
# --------------------------------------------------------------------------
# Defining parameters are ``a`` and ``1/f``; everything else is derived at
# import time (decision D-0.9). tests/unit/test_constants.py audits each
# derived value against the published NGA figure.

WGS84_A: Final[float] = 6378137.0
"""Semi-major axis, in metres. Defining parameter of the datum."""

WGS84_INV_F: Final[float] = 298.257223563
"""Inverse flattening (1/f). Defining parameter of the datum."""

WGS84_F: Final[float] = 1.0 / WGS84_INV_F
"""Flattening, derived from :data:`WGS84_INV_F`."""

WGS84_B: Final[float] = WGS84_A * (1.0 - WGS84_F)
"""Semi-minor axis, in metres. Derived: ``a(1 - f)``. Published: 6356752.314245."""

WGS84_E2: Final[float] = WGS84_F * (2.0 - WGS84_F)
"""First eccentricity squared. Derived: ``f(2 - f)``. Published: 6.694379990141e-3."""

WGS84_E: Final[float] = math.sqrt(WGS84_E2)
"""First eccentricity. Published: 8.1819190842622e-2."""

WGS84_EP2: Final[float] = WGS84_E2 / (1.0 - WGS84_E2)
"""Second eccentricity squared, ``e^2 / (1 - e^2)``. Used by the meridian arc."""

# D-1.10: source is GeographicLib, read through
# pyproj.Geod(ellps="WGS84").inv(0, 0, 0, 90). Deliberately NOT the value our
# own series computes: the point of this constant is to be an independent
# target, and seeding it from our arithmetic would make F1 compare the code
# against itself.
MERIDIAN_QUADRANT: Final[float] = 10001965.729312724
"""Meridian quadrant of the WGS84 ellipsoid, in metres.

External literal, kept fixed on purpose so that F1 has an independent
target: ``meridian_arc(90 deg)`` must reproduce it from Eq. (2) rather
than read it back from here.
"""

EQUATOR_QUADRANT: Final[float] = math.pi * WGS84_A
"""Half-equator length on the sinusoidal plane, ``pi * a``, in metres.

At ``phi = 0`` Eq. (1) reduces to ``x = lambda * a``, so a quadrant spans
``pi * a = 20 037 508.34 m`` along x.
"""

SINUSOIDAL_PROJ: Final[str] = (
    "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)
"""PROJ string of the equivalent projection.

NOTE: provided for interoperability only. The package computes the
ellipsoidal sinusoidal ("parallels plane") projection directly from
Eq. (1) and (2) of the paper; it does not call PROJ at runtime (D-0.3).
"""

PRIME_MERIDIAN_LON: Final[float] = 0.0
"""Longitude of the prime meridian, in degrees. Boundary of the triangular cells."""

ANTEMERIDIAN_LON: Final[float] = 180.0
"""Longitude of the antemeridian, in degrees. Boundary of the trapezoidal cells."""

# --------------------------------------------------------------------------
# Resolutions (artigo, Tabela 1)
# --------------------------------------------------------------------------

MIN_RESOLUTION: Final[int] = 0
MAX_RESOLUTION: Final[int] = 13

RESOLUTION_COUNT: Final[int] = MAX_RESOLUTION - MIN_RESOLUTION + 1
"""Number of resolution levels: 14, from global quadrant to 1 cm."""

QUADRANTS: Final[tuple[str, ...]] = ("NE", "NW", "SE", "SW")

Quadrant = Literal["NE", "NW", "SE", "SW"]

CELL_SIZE_M: Final[tuple[float | None, ...]] = (
    None,  # res 0 - global quadrant, no metric size
    10_000.0,
    5_000.0,
    1_000.0,
    500.0,
    100.0,
    50.0,
    10.0,
    5.0,
    1.0,
    0.5,
    0.1,
    0.05,
    0.01,
)
"""Base and height of a cell, in metres, indexed by resolution.

The cell is a parallelogram whose base and height are equal, so the
nominal area is the square of this value at every resolution.
"""

CELL_AREA_M2: Final[tuple[float | None, ...]] = (
    None,
    1e8,  # 100 km2
    2.5e7,  # 25 km2
    1e6,  # 1 km2
    250_000.0,
    10_000.0,
    2_500.0,
    100.0,
    25.0,
    1.0,
    0.25,  # 2 500 cm2
    0.01,  # 100 cm2
    0.0025,  # 25 cm2
    0.0001,  # 1 cm2
)
"""Nominal cell area, in square metres, indexed by resolution.

Nominal, not effective: trapezoidal cells at the antemeridian depart from
this value by design (artigo, secao 3.2 e Quadro 5, Req 21).
"""

VISUALIZATION_SCALE: Final[tuple[int | None, ...]] = (
    None,
    100_000_000,
    50_000_000,
    10_000_000,
    5_000_000,
    1_000_000,
    500_000,
    100_000,
    50_000,
    10_000,
    5_000,
    1_000,
    500,
    100,
)
"""Denominator of the visualization scale (Jenny et al., 2008).

Minimum line visible on paper taken as 0.1 mm, hence
``scale = cell_size_m / 1e-4``.
"""

ANALYSIS_SCALE: Final[tuple[int | None, ...]] = (
    None,
    20_000_000,
    10_000_000,
    2_000_000,
    1_000_000,
    200_000,
    100_000,
    20_000,
    10_000,
    2_000,
    1_000,
    200,
    100,
    20,
)
"""Denominator of the analysis scale (Tobler, 1987).

Sampling theory, hence ``scale = cell_size_m * 2 * 1000``.
"""

TOKENIZABLE_RESOLUTIONS: Final[tuple[int, ...]] = (1, 3, 5, 7, 9, 11, 13)
"""Resolutions whose cell size is an exact power of ten metres.

Resolution 1 and every odd resolution have sides of 10^4, 10^3, 10^2,
10^1, 10^0, 10^-1 and 10^-2 metres, so their areas are exact powers of
ten in square metres. These are the levels where one token maps onto a
whole standard metric unit of area, which is the "decimal convergence"
design criterion (artigo, Quadro 1 e Quadro 3, Blockchain Integration).

Even resolutions carry sides of the form 5 x 10^k and areas of the form
25 x 10^k, which are regular but not decimal units.
"""


class ResolutionSpec(NamedTuple):
    """One row of Table 1, addressable by resolution."""

    resolution: int
    cell_size_m: float | None
    cell_area_m2: float | None
    refinement_ratio: int | None
    alphabet: tuple[str, ...] | None
    visualization_scale: int | None
    analysis_scale: int | None


RES1_CELLS_X: Final[int] = 2003
"""Number of FULL 10 km cells along the equator, per quadrant.

``EQUATOR_QUADRANT / 10 km = 2003.75``, so columns ``0000``..``2002`` are
full and column ``2003`` is 0.75 wide. This is what reconciles the two
statements the paper makes: "approximately 2,003 cells" (secao 3.1) counts
the full ones, while ``RES1_MAX_INDEX`` (Tabela 1) addresses the partial
one as well. Source of the antemeridian trapezoids handled in F4.
"""

RES1_CELLS_Y: Final[int] = 1000
"""Number of FULL 10 km cells along the central meridian, per quadrant.

``MERIDIAN_QUADRANT / 10 km = 1000.20``, so rows ``0000``..``0999`` are
full and row ``1000`` is 0.20 tall. Same reading as :data:`RES1_CELLS_X`.
"""

RES1_MIN_INDEX: Final[str] = "0000/0000"
"""Lowest resolution-1 address of a quadrant (artigo, Tabela 1)."""

RES1_MAX_INDEX: Final[str] = "2003/1000"
"""Highest resolution-1 address of a quadrant (artigo, Tabela 1)."""

# --------------------------------------------------------------------------
# Refinement alphabets (artigo, secao 3.1 e Figura 3)
# --------------------------------------------------------------------------

QUATERNARY_GRID_SIZE: Final[int] = 2
"""Side of the 1-to-4 refinement grid: 2 x 2."""

QUINARY_GRID_SIZE: Final[int] = 5
"""Side of the 1-to-25 refinement grid: 5 x 5."""

QUATERNARY_CODES: Final[tuple[str, ...]] = ("1", "2", "3", "4")
"""1-to-4 refinement codes, used to address EVEN resolutions (2, 4, ..., 12).

Row-major from the SOUTHERN row, as drawn in Figure 3(c): ``1`` and ``2``
sit on the bottom row (south), ``3`` and ``4`` on the top row (north).
That layout is what makes the F6 rule work: vertical adjacency adds or
subtracts 2, horizontal adjacency adds or subtracts 1.
"""

QUINARY_CODES: Final[tuple[str, ...]] = tuple(
    f"{chr(ord('A') + row)}{col + 1}"
    for row in range(QUINARY_GRID_SIZE)
    for col in range(QUINARY_GRID_SIZE)
)
"""1-to-25 refinement codes (A1..E5), addressing ODD resolutions (3, 5, ..., 13).

WARNING: the Colab prototype had this inverted. Even resolutions are
QUATERNARY, odd resolutions (>= 3) are QUINARY.

Row-major from the SOUTHERN row, as drawn in Figure 3(d): row ``A`` is
the bottom row and row ``E`` the top one, so generation starts at the
south. The vertical wrap-around ``E <-> A`` used in F6 depends on this
orientation; reversing it silently flips north and south.
"""

REFINEMENT_RATIO: Final[tuple[int | None, ...]] = (
    None,  # res 0 - quadrants, not produced by refinement
    None,  # res 1 - base cells, addressed by XXXX/YYYY
    4,
    25,
    4,
    25,
    4,
    25,
    4,
    25,
    4,
    25,
    4,
    25,
)
"""Number of children produced when descending INTO each resolution.

Even resolutions come from a 1-to-4 subdivision (side halved), odd ones
from a 1-to-25 subdivision (side divided by five). Alternating the two is
what yields an exact decimal decade every two levels.
"""

REFINEMENT_ALPHABET: Final[tuple[tuple[str, ...] | None, ...]] = tuple(
    (None if ratio is None else (QUATERNARY_CODES if ratio == 4 else QUINARY_CODES))
    for ratio in REFINEMENT_RATIO
)
"""Alphabet addressing each resolution, derived from :data:`REFINEMENT_RATIO`.

``None`` at resolutions 0 and 1, which are not addressed by a refinement
code. Deriving instead of restating keeps the even/odd rule in one place.
"""

RESOLUTION_TABLE: Final[tuple[ResolutionSpec, ...]] = tuple(
    ResolutionSpec(
        resolution=resolution,
        cell_size_m=CELL_SIZE_M[resolution],
        cell_area_m2=CELL_AREA_M2[resolution],
        refinement_ratio=REFINEMENT_RATIO[resolution],
        alphabet=(QUADRANTS if resolution == 0 else REFINEMENT_ALPHABET[resolution]),
        visualization_scale=VISUALIZATION_SCALE[resolution],
        analysis_scale=ANALYSIS_SCALE[resolution],
    )
    for resolution in range(RESOLUTION_COUNT)
)
"""Table 1 as a sequence indexed by resolution.

Assembled from the tuples above rather than restated, so there is exactly
one place to correct if a row is ever found wrong.
"""


def refinement_alphabet(resolution: int) -> tuple[str, ...]:
    """Return the refinement alphabet addressing ``resolution``.

    This is the only place in the package where the even/odd rule is
    written, which is the whole point: ``B-0.1`` was a second, inverted
    copy of it living in the parser.

    Args:
        resolution: Resolution level, 2 to 13.

    Returns:
        :data:`QUATERNARY_CODES` for even resolutions,
        :data:`QUINARY_CODES` for odd ones.

    Raises:
        ResolutionError: If ``resolution`` is outside 0..13, or is 0 or 1.
            Resolution 0 is addressed by a quadrant code and resolution 1
            by an ``XXXX/YYYY`` pair; neither is a refinement.
    """
    if not isinstance(resolution, int) or isinstance(resolution, bool):
        raise ResolutionError(
            f"resolution must be an int, got {type(resolution).__name__}"
        )
    if not MIN_RESOLUTION <= resolution <= MAX_RESOLUTION:
        raise ResolutionError(
            f"resolution {resolution} outside " f"{MIN_RESOLUTION}..{MAX_RESOLUTION}"
        )
    alphabet = REFINEMENT_ALPHABET[resolution]
    if alphabet is None:
        raise ResolutionError(
            f"resolution {resolution} is not addressed by a refinement "
            "code: 0 uses a quadrant code, 1 uses an XXXX/YYYY pair"
        )
    return alphabet


# --------------------------------------------------------------------------
# Index syntax (artigo, secao 3.1)
# --------------------------------------------------------------------------

DESCENT_OPEN: Final[str] = "("
DESCENT_CLOSE: Final[str] = ")"
SIBLING_SEPARATOR: Final[str] = ","
RES1_SEPARATOR: Final[str] = "/"
RES1_DIGITS: Final[int] = 4

QUADRANT_CODE_LENGTH: Final[int] = 2
"""Length of the leading quadrant code: ``NE``, ``NW``, ``SE`` or ``SW``."""

INDEX_EXAMPLE_ATOMIC: Final[str] = "SE(1400/0374(3(C2(3))))"
"""Canonical single-cell example from the paper (secao 3.1)."""

INDEX_EXAMPLE_SIBLINGS: Final[str] = "...4(C1,C2)"
"""Canonical sibling example from the paper (secao 4): two cells, one index."""

# No index component is ever negative: quadrant mirroring absorbs the sign
# (artigo, secao 3.1).

# --------------------------------------------------------------------------
# Antemeridian extension zones (artigo, secao 3.2 e Figura 5)
# --------------------------------------------------------------------------

ExtensionZone = Literal["FIJI", "CHUKOTKA"]


class ExtensionZoneSpec(NamedTuple):
    """Bounds of one eastern-quadrant extension across the antemeridian."""

    name: str
    quadrant: str
    lon_limit: float
    lat_min: float
    lat_max: float
    description: str


EXTENSION_ZONES: Final[Mapping[str, ExtensionZoneSpec]] = MappingProxyType(
    {
        "FIJI": ExtensionZoneSpec(
            name="FIJI",
            quadrant="SE",
            lon_limit=-178.0,  # extends east quadrant to 178 W
            lat_min=-21.5,
            lat_max=-15.5,
            description="Fiji Islands",
        ),
        "CHUKOTKA": ExtensionZoneSpec(
            name="CHUKOTKA",
            quadrant="NE",
            lon_limit=-169.5,  # extends east quadrant to 169.5 W
            lat_min=64.0,
            lat_max=72.0,
            description=(
                "Russian mainland (Chukotka AO), Wrangel Island " "and nearby islands"
            ),
        ),
    }
)
"""Eastern-quadrant extensions across the antemeridian.

Antarctica along the 180th meridian is deliberately excluded (limited
cadastral application). Cells in the oceans and in Antarctica along that
boundary therefore have unequal areas.

The keys ``"FIJI"`` and ``"CHUKOTKA"`` are load-bearing: F4 tests
``extension_zone_for_point`` against these exact strings.
"""

EXTENSION_LON_PRECISION: Final[float] = 0.5
"""Longitude precision adopted for extension limits, in degrees.

The coarsest step that clears every landmass except Antarctica, hence
the .5 in the Fiji and Chukotka limits.
"""

# --------------------------------------------------------------------------
# Cell shapes (artigo, secao 3.2 e Figuras 4 e 6)
# --------------------------------------------------------------------------

CellShape = Literal["parallelogram", "triangle", "trapezoid"]

CELL_SHAPES: Final[tuple[str, ...]] = ("parallelogram", "triangle", "trapezoid")
"""Every shape a cell may take, in order of decreasing frequency."""

PARALLELOGRAM_BASE_ANGLE_DEG: Final[float] = 45.0
"""Acute angle of the cell on the projection plane, before ellipsoidal distortion."""

TRIANGLE_BASE_TO_HEIGHT_RATIO: Final[float] = 2.0
"""Base-to-height ratio of the isosceles prime-meridian cell (Figura 4a).

The cell index sits at the midpoint of the base and the cell is mirrored
about the meridian, so base = 2 x height keeps its area equal to that of
a standard parallelogram of the same resolution.
"""
