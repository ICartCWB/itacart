"""Shared pytest fixtures.

Reference points and indices transcribed from the paper and from the
Colab prototype, so the suite pins behaviour against published values.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

import itacart
from tests.reference_points import REFERENCE_POINTS

# The suite must measure the working tree, never an installed copy that
# happens to be earlier on sys.path. Without this the whole session can
# pass green against code that is not the code under edit. docs/source/
# conf.py refuses the documentation build on the same grounds; the two
# guards exist so the property holds by construction rather than by
# remembering to check.

ROOT = Path(__file__).resolve().parents[1]
_loaded_from = Path(itacart.__file__).resolve().parent
if _loaded_from != (ROOT / "src" / "itacart").resolve():  # pragma: no cover
    raise RuntimeError(
        f"itacart was imported from {_loaded_from}, not from the working tree "
        f"at {ROOT / 'src' / 'itacart'}. An installed copy is shadowing it; "
        "uninstall it, reinstall in editable mode, or set PYTHONPATH to "
        "the working tree's src/."
    )

# Windows dev environment: a corrupt cert in the store can break aiohttp
# imports pulled in by optional geo dependencies.
if sys.platform == "win32":  # pragma: no cover
    import ssl

    ssl.SSLContext._load_windows_store_certs = (  # type: ignore[method-assign]
        lambda self, storename, purpose: []
    )


@pytest.fixture
def sydney_opera_house() -> tuple[float, float]:
    """(lon, lat) of the Sydney Opera House."""
    return (151.2150784, -33.8567529)


@pytest.fixture
def liberty_statue() -> tuple[float, float]:
    """(lon, lat) of the Statue of Liberty."""
    return (-74.0445142, 40.6892077)


@pytest.fixture
def praca_da_se() -> tuple[float, float]:
    """(lon, lat) of Praça da Sé, São Paulo."""
    return (-46.6328862, -23.5508962)


@pytest.fixture
def central_park_index() -> str:
    """Figure 7 of the paper: compositional fill at resolutions 6 and 7."""
    return (
        "NW(0625/0451(1(E1(3(B2(4(A2,B2,B3,B4,C2,C3,C4,C5,D1,D2,D3,D4,"
        "D5,E1,E2,E3,E4,E5)),B3(3(C1,D1,D2,E1,E2,E3,E4)),C2(1(A5,B5,C5,"
        "D4,D5,E4,E5),2,3(A4,A5,B3,B4,B5,C3,C4,C5,D2,D3,D4,D5,E2,E3,E4,"
        "E5),4),C3(1,2(A1,B1,B2,B3,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,"
        "E5),3,4),D2(1(A1,A2,A3,A4,A5,B2,B3,B4,B5,C3,C4,C5,D5),2,4(A3,"
        "A4,A5,B5)),D3(1,2,3(A1,A2,A3,A4,A5,B1,B2,B3,B4,B5,C2,C3,C4,C5,"
        "D4,D5),4),D4(1(D1,E1),3(A1,B1,C1,D1,E1)))))))"
    )


@pytest.fixture
def sydney_cell() -> str:
    """Atomic index for the Sydney Opera House, resolution 7.

    The string is what ``geo_to_cell`` returns for the opera house at
    resolution 7, and ``test_index.py`` has pinned it as a seven-level
    path since F2. The docstring said nine, which is the discrepancy
    ``test_fixture_invariants.py`` now forbids for every cell fixture.
    """
    return "SE(1400/0374(3(C2(3(C2(4(C1)))))))"


@pytest.fixture
def paper_example_index() -> str:
    """Example index from section 3.1 of the paper, resolution 4."""
    return "SE(1400/0374(3(C2(3))))"


@pytest.fixture
def suva() -> tuple[float, float]:
    """(lon, lat) of Suva, Fiji: east of 180 and inside the Fiji band."""
    return REFERENCE_POINTS["suva"]


@pytest.fixture
def wrangel() -> tuple[float, float]:
    """(lon, lat) of Wrangel Island: west of 180 and inside the Chukotka band."""
    return REFERENCE_POINTS["wrangel"]


@pytest.fixture
def greenwich_observatory() -> tuple[float, float]:
    """(lon, lat) of the Royal Observatory, metres west of the prime meridian."""
    return REFERENCE_POINTS["greenwich_observatory"]


@pytest.fixture
def central_pacific() -> tuple[float, float]:
    """(lon, lat) of open ocean in neither extension zone."""
    return REFERENCE_POINTS["central_pacific"]
