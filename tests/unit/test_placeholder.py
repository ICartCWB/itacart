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
def test_hierarchy_is_still_a_stub() -> None:
    """The one entry point of this file that has not landed yet.

    Kept as the file's own reminder that ``hierarchy`` is unimplemented.
    When F5 lands, this turns red rather than passing quietly, which is
    the failure mode ``test_geo_to_cell_stub`` had: it was written in F0
    against a ``geo_to_cell`` that raised, F3 implemented the function,
    and the non-strict marker turned the pass into a silent ``xpass``
    that nobody read.
    """
    itacart.get_parent("SE(1400/0374)")


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
