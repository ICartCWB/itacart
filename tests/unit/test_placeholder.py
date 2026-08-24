"""Sanity checks that hold before any implementation lands."""

from __future__ import annotations

import pytest

import itacart


def test_version_is_exposed() -> None:
    assert itacart.__version__


def test_every_exported_name_resolves() -> None:
    missing = [n for n in itacart.__all__ if not hasattr(itacart, n)]
    assert missing == []


def test_paper_doi_is_pinned() -> None:
    assert itacart.__paper_doi__ == "10.14393/rbcv77n0a-79281"


@pytest.mark.xfail(raises=NotImplementedError, reason="stub")
def test_geo_to_cell_stub(praca_da_se: tuple[float, float]) -> None:
    lon, lat = praca_da_se
    itacart.geo_to_cell(lon, lat, 13)
