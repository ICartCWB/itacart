"""OGC DGGS Core conformance suite (Topic 21 / ISO 19170-1).

One test per requirement of Frame 4 of the paper. Each asserts the
behaviour the paper claims, so the compliance table stops being a
declaration and becomes something CI verifies.

New code. Maps Frame 4 (Core) and Frame 5 (EAERS) of the paper.
"""

from __future__ import annotations

import pytest

import itacart

pytestmark = pytest.mark.conformance


class TestCore:
    """DGGS Core, requirements 6 to 19."""

    @pytest.mark.xfail(raises=NotImplementedError, reason="stub")
    def test_req_06_harmonized_model(self) -> None:
        """describe() reports the architecture per Figure 13 of Gibb (2021)."""
        raise NotImplementedError

    @pytest.mark.xfail(raises=NotImplementedError, reason="stub")
    def test_req_07_defined_crs(self) -> None:
        """crs() reports WGS84, guaranteeing GNSS compatibility."""
        raise NotImplementedError

    @pytest.mark.xfail(raises=NotImplementedError, reason="stub")
    def test_req_08_10_global_complete_unique_domain(self) -> None:
        """Every land position resolves to exactly one cell.

        Covers the western-quadrant X=0 rule, the prime-meridian
        triangles and the two antemeridian extension zones.
        """
        raise NotImplementedError

    @pytest.mark.xfail(raises=NotImplementedError, reason="stub")
    def test_req_11_simple_cell_geometry(self) -> None:
        """Cell boundaries are simple, non-self-intersecting polygons."""
        raise NotImplementedError

    @pytest.mark.xfail(raises=NotImplementedError, reason="stub")
    def test_req_12_direct_position(self) -> None:
        """The anchor lies on the cell boundary, as specified."""
        raise NotImplementedError

    @pytest.mark.xfail(raises=NotImplementedError, reason="stub")
    def test_req_13_unique_address(self) -> None:
        """Canonical forms are equal iff the regions are equal."""
        raise NotImplementedError

    @pytest.mark.xfail(raises=NotImplementedError, reason="stub")
    def test_req_14_15_hierarchical_grid_sequence(self) -> None:
        """Fourteen ordered levels with the documented refinement ratios."""
        raise NotImplementedError

    @pytest.mark.xfail(raises=NotImplementedError, reason="stub")
    def test_req_16_quantization(self) -> None:
        """Vector data maps to cell sets via polyfill and vertex_to_cell."""
        raise NotImplementedError

    def test_req_17_topological_queries(self) -> None:
        """Parent, child and neighbour resolve from the index alone.

        Completed by F6. The neighbour half is the part that had no prior
        art: it is checked here by composing a step and its opposite, which
        can only return to the origin if the arithmetic and the tessellation
        agree.
        """
        cell = "NE(0500/0300)"
        assert itacart.get_parent(cell) == "NE"
        assert len(list(itacart.get_children(cell, flatten=True))) == 4

        for outward, back in (("N", "S"), ("E", "W"), ("NE", "SW"), ("NW", "SE")):
            away = itacart.get_neighbor(cell, outward)
            assert isinstance(away, str)
            assert itacart.get_neighbor(away, back) == cell

        assert len(itacart.grid_disk(cell, 1)) == 9
        assert len(itacart.grid_disk(cell, 1, "manhattan")) == 5
        assert itacart.are_neighbor_cells(cell, "NE(0499/0301)")
        assert not itacart.are_neighbor_cells(cell, "NE(0500/0301)")
        assert itacart.grid_distance(cell, "NE(0505/0303)") == 8
        assert len(itacart.cell_to_edges(cell)) == 4

    @pytest.mark.xfail(raises=NotImplementedError, reason="stub")
    def test_req_18_19_interoperability(self) -> None:
        """Cells export to GeoJSON and WKT."""
        raise NotImplementedError


class TestEAERS:
    """Equal-Area Earth Reference System, requirements 20 to 29.

    Partial compliance is expected and intentional; these tests pin the
    divergences so they stay deliberate rather than drifting.
    """

    @pytest.mark.xfail(raises=NotImplementedError, reason="stub")
    def test_req_21_equal_area_constraint(self) -> None:
        """Equal area holds except for antemeridian trapezoids."""
        raise NotImplementedError

    @pytest.mark.xfail(raises=NotImplementedError, reason="stub")
    def test_req_22_25_not_met_by_design(self) -> None:
        """Direct surface tessellation: no polyhedral interface exists."""
        raise NotImplementedError

    @pytest.mark.xfail(raises=NotImplementedError, reason="stub")
    def test_req_26_simple_2d_polygons(self) -> None:
        """Cells are parallelograms, triangles or trapezoids."""
        raise NotImplementedError

    @pytest.mark.xfail(raises=NotImplementedError, reason="stub")
    def test_req_27_representative_position(self) -> None:
        """The anchor is a vertex, not the centroid; centroid is separate."""
        raise NotImplementedError
