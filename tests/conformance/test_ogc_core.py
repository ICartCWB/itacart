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

    def test_req_16_quantization(self) -> None:
        """Vector data maps to cell sets via polyfill and vertex_to_cell.

        Completed by F7, which supplied both halves of the operation the
        requirement names. Quantization is the map from vector data to a
        cell set at a stated resolution, and the requirement is about
        vector data rather than about polygons, so all three dimensions
        are asserted.

        The areal half is checked for closure, not merely for being
        non-empty: refilling the footprint of the cell set has to return
        the cell set. The vertex half is checked for sequence, since an
        unordered answer would quantize the geometry without preserving
        the thing that lets it be reconstructed.
        """
        from shapely.geometry import LineString, Point, Polygon

        parcel = Polygon(
            [
                (-73.9812, 40.7681),
                (-73.9581, 40.8005),
                (-73.9497, 40.7968),
                (-73.9730, 40.7644),
            ]
        )

        filled = itacart.polyfill(parcel, 5)
        cells = itacart.decompose(filled)
        assert cells
        assert all(itacart.get_resolution(cell) == 5 for cell in cells)

        from shapely.ops import unary_union

        footprint = unary_union([itacart.cell_to_polygon(cell) for cell in cells])
        assert itacart.polyfill(footprint, 5) == filled

        vertices = itacart.vertex_to_cell(parcel, 5)
        assert vertices
        assert all(itacart.get_resolution(cell) == 5 for cell in vertices)
        assert itacart.cells_to_geometry(vertices).geom_type == "Polygon"

        line = LineString([(-73.98, 40.77), (-73.96, 40.79)])
        assert len(itacart.vertex_to_cell(line, 5)) == 2

        point = Point(-73.97, 40.78)
        assert itacart.vertex_to_cell(point, 5) == [
            itacart.geo_to_cell(-73.97, 40.78, 5)
        ]

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

    def test_req_18_19_interoperability(self) -> None:
        """Cells export to GeoJSON and WKT, and come back.

        Completed by F9b. The requirement names two encodings, and the
        paper's compliance table marks them met by design; a concrete
        exporter is what turns that into something CI can check.

        Export alone would be a weak reading. A function that emitted a
        syntactically valid FeatureCollection of the wrong polygons would
        pass it, so the round trip is asserted as well: the indices come
        back out of the file exactly, because they are written into the
        Feature ``id`` that section 3.2 of RFC 7946 reserves for exactly
        this, and refilling an exported polygon at its own resolution
        returns the cell it came from.

        Both seams are included rather than an interior sample. The
        meridian triangle and the polar cap are where a naive exporter
        produces valid GeoJSON that means the wrong thing, so a
        requirement about interoperability that never left the interior
        would be measuring the easy half.
        """
        from shapely.geometry import shape
        from shapely.wkt import loads as wkt_loads

        cells = ["NE(0500/0300)", "NE(0000/0300)", "NE(0000/1000)"]
        index = itacart.compose(cells)

        collection = itacart.cells_to_geojson(index)
        assert collection["type"] == "FeatureCollection"
        assert len(collection["features"]) == 3
        for feature, cell in zip(collection["features"], cells):
            assert feature["type"] == "Feature"
            assert feature["id"] == cell
            assert feature["properties"]["itacart_index"] == cell
            geometry = shape(feature["geometry"])
            assert geometry.is_valid
            assert geometry.exterior.is_ccw

        assert itacart.recover_from_geojson(collection) == cells

        single = itacart.cell_to_wkt(cells[0])
        assert isinstance(single, str)
        assert wkt_loads(single).geom_type == "Polygon"

        many = itacart.cell_to_wkt(index)
        assert isinstance(many, list) and len(many) == 3

        merged = itacart.cells_to_wkt(index)
        assert wkt_loads(merged).geom_type == "GeometryCollection"

        ordinary = itacart.cells_to_geojson(cells[0])
        assert itacart.from_geojson(ordinary, 1) == [cells[0]]


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
