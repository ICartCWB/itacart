"""Sanity checks that hold before any implementation lands."""

from __future__ import annotations

import itacart


def test_version_is_exposed() -> None:
    assert itacart.__version__


def test_every_exported_name_resolves() -> None:
    missing = [n for n in itacart.__all__ if not hasattr(itacart, n)]
    assert missing == []


def test_paper_doi_is_pinned() -> None:
    assert itacart.__paper_doi__ == "10.14393/rbcv77n0a-79281"


def test_get_parent_climbs_to_the_quadrant() -> None:
    """The former stub, as an assertion.

    It was written against a ``get_parent`` that raised, and stayed
    marked as expected-to-fail after the function landed. The marker was
    not strict, so the pass registered as a silent ``xpass`` that nobody
    read, which is the same failure mode ``test_geo_to_cell_stub`` had
    before it. Strict expected-failure is now on for the whole suite, so
    a marker that stops describing reality turns red instead of quiet.

    A resolution-1 cell has the quadrant as its parent, and the quadrant
    is where the climb stops: it is resolution 0 and has no coordinate
    pair to truncate.
    """
    cell = "SE(1400/0374)"
    assert itacart.get_resolution(cell) == 1
    assert itacart.get_parent(cell) == "SE"
    assert itacart.get_parent(cell, 0) == "SE"


def test_geo_to_cell_reaches_the_finest_resolution(
    praca_da_se: tuple[float, float],
) -> None:
    """The former stub, as an assertion.

    Resolution 13 is a one-centimetre cell, which is where the quantizer
    is most exposed to the accumulated rounding of thirteen refinement
    steps. Round-tripping the anchor back through ``geo_to_cell`` is what
    proves the descent landed in the cell it named rather than in a
    neighbour.
    """
    lon, lat = praca_da_se
    cell = itacart.geo_to_cell(lon, lat, 13)
    assert itacart.is_valid_index(cell)
    assert itacart.get_resolution(cell) == 13
    assert itacart.geo_to_cell(*itacart.cell_to_anchor(cell), 13) == cell
