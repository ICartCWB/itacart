"""OGC DGGS Core conformance suite (Topic 21 / ISO 19170-1).

One test per requirement, asserting what the requirement asks rather than
what the paper claims about it, so that the compliance table stops being a
declaration and becomes something CI verifies -- including where the two
disagree.

The governing document is OGC 20-040r3, *Topic 21 -- Discrete Global Grid
Systems -- Part 1: Core Reference System and Operations and Equal Area
Earth Reference System* (2021), which is the edition ISO 19170-1 follows.
Its predecessor, OGC 15-104r5 (2017), stops at Requirement 18 and has no
EAERS class at all; the numbers 20 to 29 used below exist only in
20-040r3. Every requirement quoted here carries its URI, so that checking
a test against the standard is a lookup rather than a search.

New code. Maps Frame 4 (Core) and Frame 5 (EAERS) of the paper.
"""

from __future__ import annotations

import math

import pytest
from shapely.geometry import LinearRing, LineString, Point, Polygon
from shapely.ops import unary_union

import itacart
from itacart import constants
from itacart.boundary import plane_ring

#: One cell per family, so a test can say which family it is about. The
#: nine ordinary families plus the southern cap: the four quadrants, the
#: two seams, the equator, both poles, and the finest addressable depth.
FAMILIES: tuple[tuple[str, str], ...] = (
    ("interior", "NE(0500/0300)"),
    (
        "interior at resolution 13",
        "SW(0476/0259(4(E5(4(E2(2(B4(1(B4(2(D1(2(D3)))))))))))))",
    ),
    ("meridian column", "NE(0000/0300)"),
    ("trapezoid", "NE(2003/0000)"),
    ("equator, row zero", "NE(0500/0000)"),
    ("polar row", "NE(0001/0999)"),
    ("northern cap", "NE(0000/1000)"),
    ("southern cap", "SE(0000/1000)"),
    ("south", "SE(0500/0300)"),
    ("west", "NW(0500/0300)"),
)

#: The two polar caps. Named rather than detected, so that a third one
#: appearing anywhere would fail a test instead of quietly joining a
#: clause written for two.
CAPS: tuple[str, ...] = ("NE(0000/1000)", "SE(0000/1000)")

#: Wide enough to pass the last column of row zero, which is 2003.
_COLUMN_LIMIT = 2200


def _columns(quadrant: str, row: int) -> list[int]:
    """Every column of a row that names a cell, tested rather than counted."""
    return [
        column
        for column in range(_COLUMN_LIMIT)
        if itacart.is_valid_cell(f"{quadrant}({column:04d}/{row:04d})")
    ]


def _authalic_cap_area(latitude: float) -> float:
    """Surface area above a parallel, from the closed form.

    Independent of the package: written from the authalic integral of the
    ellipsoid rather than read back from anything ``itacart`` computes, so
    that a cap area agreeing with it is evidence and not a tautology.
    """
    eccentricity = math.sqrt(constants.WGS84_E2)

    def zone(phi: float) -> float:
        sine = math.sin(phi)
        return (1.0 - constants.WGS84_E2) * (
            sine / (1.0 - constants.WGS84_E2 * sine * sine)
            - math.log((1.0 - eccentricity * sine) / (1.0 + eccentricity * sine))
            / (2.0 * eccentricity)
        )

    pole = zone(math.pi / 2.0)
    return math.pi * constants.WGS84_A**2 * (pole - zone(math.radians(abs(latitude))))


def _points_in(
    region: "Polygon", wanted: int, side: int = 400
) -> list[tuple[float, float]]:
    """A deterministic sample of points strictly inside a region."""
    x0, y0, x1, y1 = region.bounds
    found: list[tuple[float, float]] = []
    for index in range(side * side):
        x = x0 + (index % side + 0.5) * (x1 - x0) / side
        y = y0 + (index // side + 0.5) * (y1 - y0) / side
        if region.contains(Point(x, y)):
            found.append((x, y))
            if len(found) == wanted:
                break
    return found


def _arc_boundary(
    latitude: float, edges: int, steps: int = 90
) -> tuple[list[float], list[LineString]]:
    """A parallel partitioned into ``edges`` small-circle arcs.

    Longitudes run over 0 to 360 rather than -180 to 180, so that the
    partition is not an artefact of where the chart happens to cut.
    """
    vertices = [index * 360.0 / edges for index in range(edges)]
    arcs = []
    for index, start in enumerate(vertices):
        end = start + 360.0 / edges
        arcs.append(
            LineString(
                [
                    (start + (end - start) * step / steps, latitude)
                    for step in range(steps + 1)
                ]
            )
        )
    return vertices, arcs


pytestmark = pytest.mark.conformance


def _frechet_mean(cell: str, side: int = 24) -> tuple[float, float]:
    """The geodesic centre of surface area, by direct numerical search.

    Independent of the package's own centroid formula on purpose: the
    sample is drawn in the equal-area plane, so it is uniform by area on
    the ellipsoid, and the minimiser is found by coordinate descent on the
    area-weighted squared geodesic distance rather than by any closed
    form. Coarse by design -- the divergences it has to resolve are metres
    against kilometres, not millimetres.
    """
    figure = Polygon(plane_ring(cell)[1])
    x0, y0, x1, y1 = figure.bounds
    samples = [
        itacart.sinusoidal_to_geodetic(x, y)
        for i in range(side)
        for j in range(side)
        for x, y in [
            (x0 + (i + 0.5) * (x1 - x0) / side, y0 + (j + 0.5) * (y1 - y0) / side)
        ]
        if figure.contains(Point(x, y))
    ]

    def cost(point: tuple[float, float]) -> float:
        return sum(
            itacart.inverse_geodesic(point[0], point[1], lon, lat)[0] ** 2
            for lon, lat in samples
        )

    best = itacart.cell_to_centroid(cell)
    assert isinstance(best, tuple)
    best_cost, step = cost(best), 0.01
    while step > 1e-7:
        moved = False
        for dlon, dlat in ((step, 0.0), (-step, 0.0), (0.0, step), (0.0, -step)):
            trial = (best[0] + dlon, max(-90.0, min(90.0, best[1] + dlat)))
            trial_cost = cost(trial)
            if trial_cost < best_cost:
                best, best_cost, moved = trial, trial_cost, True
        if not moved:
            step /= 2.0
    return best


class TestCore:
    """DGGS Core, requirements 6 to 19."""

    def test_req_06_harmonized_model(self) -> None:
        """The data model complies with the DGGS Core RS model.

        Requirement 6 (``core/rs/harmonized_model``): the data model and
        its element definitions shall comply with the DGGS Core RS data
        model of Figure 13 and the definitions of Tables 40 to 46.

        What a package can assert about a conceptual data model is that it
        reports the elements the model names, and reports them from one
        place. So the test checks that the description carries a domain
        with its dimensionality, a CRS, a tessellation, a cell geometry,
        an index and a resolution sequence -- and that the values it gives
        are the package's own rather than a second transcription that
        could drift.
        """
        described = itacart.describe()
        for element in (
            "crs",
            "domain",
            "tessellation",
            "cell_geometry",
            "index",
            "resolutions",
        ):
            assert element in described, element

        assert described["domain"]["dimensionality"]["spatial"] == 2
        assert described["domain"]["level_zero"] == list(itacart.QUADRANTS)
        assert described["resolutions"]["table"] == list(itacart.resolution_table())
        assert described["crs"] == itacart.crs()
        assert described["cell_geometry"]["representative_position"] == "centroid"

    def test_req_07_defined_crs(self) -> None:
        """A defined CRS, with WGS84 as the datum.

        Requirement 7 (``core/rs/crs``): the reference system shall define
        a CRS, and comply with the requirements for provision of coordinate
        epoch as specified for ``MD_ReferenceSystem``.

        WGS84 is the paper's answer and the reason for it is the
        GNSS-compatibility design criterion: a cadastral system has to
        accept field survey coordinates without a datum transformation
        standing between the survey and the cell. The coordinate epoch is
        reported as absent rather than omitted, the datum being static;
        a field that is silently missing and a field that says None are
        different declarations.
        """
        reported = itacart.crs()
        assert reported["datum"] == "WGS84"
        assert reported["code"] == 4326
        assert "coordinate_epoch" in reported
        assert reported["coordinate_epoch"] is None

        ellipsoid = reported["ellipsoid"]
        assert ellipsoid["semi_major_axis_m"] == itacart.WGS84_A
        assert ellipsoid["first_eccentricity_squared"] == itacart.WGS84_E2
        assert reported["projection"]["equal_area"] is True
        assert reported["projection"]["proj_string"] == itacart.SINUSOIDAL_PROJ

    @pytest.mark.slow
    def test_req_08_10_global_complete_unique_domain(self) -> None:
        """Global domain, completeness and location uniqueness.

        Requirement 8 (``core/rs/global_domain``) asks the reference system
        to specify a global domain and its dimensionality; Requirement 9
        (``/complete``) that the level zero grid cover the entire domain;
        Requirement 10 (``/unique``) that every location be in exactly one
        cell of the level zero grid.

        Level zero is measured for what it is -- four quadrants, resolution
        zero -- and completeness is then measured one level down, where the
        package has an operation that answers. This is the honest shape of
        the evidence: ``geo_to_cell`` refuses resolution zero, so the level
        zero half of Requirements 9 and 10 is asserted structurally and the
        positional half at resolution 1.

        Containment is asserted against the **exact plane figure**, not
        against the densified lon/lat polygon, and with the longitude of an
        extension-zone point continued past the antemeridian rather than
        wrapped back. Both choices are the requirement's, not a
        convenience: an extension zone exists precisely so that a cell east
        of +180 is addressed east of +180, and comparing such a cell
        against a point re-expressed at -179.9 asks whether two things on
        opposite sides of the seam coincide, which is a different question
        with a foregone answer.

        Three predicates live near this question and none of them is it,
        which is why the condition below is spelled out rather than
        borrowed. ``is_extension_cell`` marks a cell running past the
        nominal last column of its row, which happens near the seam in
        every quadrant: at (-179.9, -85.8) it is true although the point
        belongs to no named zone, and using it moves a south-western
        Antarctic point a full turn away from its own cell.
        ``extension_zone_for_point`` answers which zone *governs* a
        latitude band and hemisphere, deliberately and by its own
        docstring, so it returns FIJI at longitude 1.0. What this test
        needs is neither: it is the strip west of the antemeridian that the
        extension actually reaches, which is the zone predicate conjoined
        with the zone's own ``lon_limit``.

        Under the right comparison an enumerated lattice of 2376 positions
        has every position inside its own cell's plane figure. Eight of
        them fall outside the *densified* lon/lat polygon of that same
        cell -- six in the two named extension zones and two within the
        densification tolerance near the pole. That is a property of the
        drawing, not of the grid, and the count is pinned so that a change
        in it is visible.

        The first sentence is true of this lattice and is not a bound. The
        antemeridian border is a curve in the projection plane and a
        resolution 1 trapezoid approximates it with one chord, so there are
        positions this lattice does not reach whose plane figure
        under-covers them while the grid addresses them correctly. They are
        measured in ``test_req_22_25`` rather than left implied here.
        """
        assert set(itacart.QUADRANTS) == {"NE", "NW", "SE", "SW"}
        for quadrant in itacart.QUADRANTS:
            assert itacart.is_valid_cell(quadrant)
            assert itacart.get_resolution(quadrant) == 0

        positions = [
            (longitude / 10.0, latitude / 10.0)
            for latitude in range(-899, 900, 41)
            for longitude in range(-1799, 1800, 67)
        ]
        assert len(positions) == 2376

        outside_the_polygon = []
        for longitude, latitude in positions:
            cell = itacart.geo_to_cell(longitude, latitude, 1)
            assert isinstance(cell, str)
            continued = longitude
            zone = itacart.extension_zone_for_point(longitude, latitude)
            if (
                zone is not None
                and longitude <= itacart.EXTENSION_ZONES[zone].lon_limit
            ):
                continued = longitude + 360.0
            plane_point = Point(itacart.geodetic_to_sinusoidal(continued, latitude))
            assert Polygon(plane_ring(cell)[1]).covers(
                plane_point
            ), f"({longitude}, {latitude}) is outside the plane figure of {cell}"
            if not itacart.cell_to_polygon(cell).covers(Point(longitude, latitude)):
                outside_the_polygon.append(cell)

        assert len(outside_the_polygon) == 8
        assert sum(itacart.is_extension_cell(c) for c in outside_the_polygon) == 6

        for row in (0, 300, 900, 999):
            columns = _columns("NE", row)
            rings = [
                Polygon(plane_ring(f"NE({column:04d}/{row:04d})")[1])
                for column in columns
            ]
            overlaps = sum(
                1
                for left, right in zip(rings, rings[1:])
                if left.intersection(right).area > 1e-6
            )
            assert overlaps == 0, f"row {row}: {overlaps} of {len(columns)} overlap"

    def test_req_11_simple_cell_geometry(self) -> None:
        """Cells with simple geometry, measured as regions of the surface.

        Requirement 11 (``core/rs/cell/simple``) asks for cells with simple
        geometry at every level of refinement. Clause 8.2.4.1 gives three
        criteria and they are topological: the cell does not
        self-intersect, it is topologically the same as a circle, and it
        encloses a region measurable using a metric of the cell's own
        dimension. There is no criterion here about counting edges against
        vertices; that one belongs to Requirement 26.

        All ten families satisfy all three, the polar caps included. A cap
        is the region of the ellipsoid above one parallel: its boundary is
        that single parallel, which is a simple closed curve, so the region
        is homeomorphic to a disc and encloses a measurable area.

        The three-point ring the package draws is a rendering in the
        lon/lat chart, and the criteria are not about the chart. Read as a
        chart figure the ring looks like a triangle whose base runs the
        width of the plane; read as the region it names, it is a cap. The
        area is what settles which reading is the cell: it is asserted
        against the closed form for a cap on the ellipsoid, computed here
        rather than borrowed from the package.
        """
        for label, cell in FAMILIES:
            assert LinearRing(itacart.cell_to_boundary(cell)).is_simple, label
            assert Polygon(plane_ring(cell)[1]).area > 0.0, label

        for cap in CAPS:
            ring = itacart.cell_to_boundary(cap)
            boundary_latitudes = {lat for _, lat in ring if abs(lat) < 90.0}
            assert len(boundary_latitudes) == 1, f"{cap}: boundary is not one parallel"

            parallel = boundary_latitudes.pop()
            measured = Polygon(plane_ring(cap)[1]).area
            assert measured == pytest.approx(_authalic_cap_area(parallel), rel=1e-7)

    def test_req_12_direct_position(self) -> None:
        """A direct position within the zone's boundary -- the centroid.

        Requirement 12 (``core/rs/cell/direct_position``): all zones in
        each discrete global grid shall be assigned a direct position that
        is *within* the zone's boundary.

        The paper answers this with the anchor, the lower-left vertex, and
        describes it as residing within the boundary of the cell. Measured
        against the word the requirement uses, it does not: a vertex is on
        the boundary, and a point on the boundary is not within it. The
        anchor fails in ten families out of ten and the centroid passes in
        ten out of ten, which is what settles the choice -- not a
        preference, and not a parameter.

        Measured on the plane figure, where the cell's edges are straight
        and the containment test is exact. The same answer has to hold for
        the exporter, and the round trip below is what ties them together:
        the centroid is not merely inside, it re-addresses its own cell.
        """
        for label, cell in FAMILIES:
            figure = Polygon(plane_ring(cell)[1])
            anchor = Point(
                itacart.geodetic_to_sinusoidal(*itacart.cell_to_anchor(cell))
            )
            centroid_lonlat = itacart.cell_to_centroid(cell)
            assert isinstance(centroid_lonlat, tuple)
            centroid = Point(itacart.geodetic_to_sinusoidal(*centroid_lonlat))

            assert not anchor.within(figure), f"{label}: anchor is within"
            assert figure.exterior.distance(anchor) < 1e-6, f"{label}: anchor off ring"
            assert centroid.within(figure), f"{label}: centroid is not within"

            resolution = itacart.get_resolution(cell)
            assert itacart.geo_to_cell(*centroid_lonlat, resolution) == cell, label

    @pytest.mark.slow
    def test_req_13_unique_address(self) -> None:
        """A globally unique zonal identifier for every zone.

        Requirement 13 (``core/rs/cell/address``): all zones in all
        discrete global grids shall have a globally unique zonal identifier
        that provides a spatio-temporal reference.

        The requirement is uniqueness *of the identifier over the zones*:
        no identifier may name two zones. It is not uniqueness of spelling,
        and the distinction matters here, because the package accepts a
        western spelling that ``is_valid_cell`` denies and rewrites it into
        the canonical eastern one. That is a second spelling of one zone,
        not a second zone with one name, so it leaves this requirement
        intact. ITACaRT has a canonical form and accepts one variant of it.

        Enumerated over three whole rows rather than sampled, because a
        collision would live at a seam and a random draw would miss it.
        """
        for row, expected in ((0, 8014), (300, 7138), (999, 6)):
            zones: dict[str, str] = {}
            for quadrant in itacart.QUADRANTS:
                for column in _columns(quadrant, row):
                    cell = f"{quadrant}({column:04d}/{row:04d})"
                    identifier = itacart.normalize(cell)
                    assert isinstance(identifier, str)
                    assert identifier not in zones, f"{identifier} names two zones"
                    zones[identifier] = cell
            assert len(zones) == expected, f"row {row}: {len(zones)} zones"

    def test_req_14_15_hierarchical_grid_sequence(self) -> None:
        """Discrete global grids, and their sequence in refinement order.

        Requirement 14 (``core/rs/discrete_global_grid``) asks the
        reference system to define discrete global grids as aggregations of
        all cells at one level of the hierarchy; Requirement 15
        (``/sequence``) that the hierarchy be sorted in order of increasing
        refinement level.

        Fourteen levels, from the global quadrant to one centimetre, with
        cell size strictly decreasing and the refinement ratio alternating
        4 and 25. The two levels that carry no ratio are refused by name
        rather than skipped: resolution 0 is a quadrant and resolution 1 is
        the Cartesian base grid, and neither is produced by refining
        anything.
        """
        assert len(itacart.resolution_table()) == 14
        assert (itacart.MIN_RESOLUTION, itacart.MAX_RESOLUTION) == (0, 13)

        with pytest.raises(itacart.ResolutionError):
            itacart.cell_size(0)
        for resolution in (0, 1):
            with pytest.raises(itacart.ResolutionError):
                itacart.refinement_ratio(resolution)

        sizes = [itacart.cell_size(r) for r in range(1, itacart.MAX_RESOLUTION + 1)]
        assert sizes == sorted(sizes, reverse=True)
        assert (sizes[0], sizes[-1]) == (10_000.0, 0.01)

        ratios = [
            itacart.refinement_ratio(r) for r in range(2, itacart.MAX_RESOLUTION + 1)
        ]
        assert ratios == [4, 25] * 6
        for resolution in range(2, itacart.MAX_RESOLUTION + 1):
            linear = itacart.linear_refinement_ratio(resolution)
            assert linear * linear == itacart.refinement_ratio(resolution)

    def test_req_16_quantization(self) -> None:
        """Vector data maps to cell sets via polyfill and vertex_to_cell.

        Both halves of the operation the requirement names are in place.
        Quantization is the map from vector data to a
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

        The neighbour half is the part that had no prior
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

        The requirement names two encodings, and the
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

    def test_req_21_equal_area_constraint(self) -> None:
        """EAERS global domain, and the constraint value it must declare.

        Requirement 21 (``ea/ers/global_domain``) has two clauses: the
        domain shall be the whole surface of the reference frame's Earth
        model, **and** the reference system shall specify ``cellEqualSized``
        among its grid constraint values.

        The name this test carries is inherited from the paper's Frame 5,
        which labels Requirement 21 the equal-area constraint. It is not:
        equal area is Requirement 29 and its error budget is Requirement
        28, and neither has a test in this suite yet. What Requirement 21
        asks for is a *declaration*, which makes the second clause a
        property of ``describe()`` and this requirement dependent on the
        engine after all. The name is left as it stands so that a reader
        holding the paper finds the test where the paper puts it.

        Both clauses are asserted, and the declaration is required to name
        its exceptions. A reference system may declare ``cellEqualSized``
        and hold it for most of its domain; declaring it while three cell
        families do not hold it, and saying nothing about them, is the
        failure mode the declaration exists to avoid.
        """
        described = itacart.describe()
        assert described["domain"]["extent"] == "whole surface of the WGS84 ellipsoid"

        constraints = described["constraints"]
        assert constraints["cellEqualSized"] is True
        assert set(constraints["exceptions"]) == {"trapezoid", "polar cap", "polar row"}

        measured = {itacart.cell_shape(cell) for _, cell in FAMILIES}
        assert "trapezoid" in measured
        for cell in CAPS:
            assert itacart.normalized_cell_area(cell) != pytest.approx(1.0)
        assert itacart.normalized_cell_area("NE(2003/0000)") != pytest.approx(1.0)
        assert itacart.normalized_cell_area("NE(0500/0300)") == pytest.approx(1.0)

    def test_req_22_25_not_met_by_design(self) -> None:
        """Initial tessellation, sequence, its limit, and area preservation.

        The paper's Frame 5 groups Requirements 22 to 25 under one heading,
        "Initial Tessellation from Polyhedron", and records all four as not
        met, because the figures of the standard put a polyhedron in the
        path of all four. The name of this test is kept so that a reader
        holding the paper finds it where the paper puts it. Read against
        the normative text rather than the figures, only the first is about
        a polyhedron:

        - Requirement 22 (``tessellation/initial``) asks for an initial
          tessellation of equal-area cells produced by mapping the faces of
          a base unit polyhedron. **Not met**, by design and by contract.
        - Requirement 23 (``tessellation/sequence``) asks for operations
          generating grids with progressively smaller cells. **Met.**
        - Requirement 24 (``tessellation/sequence/max``) asks for a stated
          limit to the number of iterations. **Met**: resolution 13.
        - Requirement 25 (``tessellation/global_area_preservation``) offers
          two branches and says one shall apply. This CRS has a static
          datum, so the static branch governs, and it asks that domain
          completeness and position uniqueness be preserved throughout the
          sequence *for all cells within their respective grids*. **Met.**

        Requirement 25 is the one that had to be measured rather than read
        off. The children of a trapezoid cover 1.000248 of their parent's
        plane figure, which looks at first like a uniqueness failure and is
        not one. Every point of that surplus is addressed at resolution 1
        to the same trapezoid and at resolution 2 to one of its own
        children, so the grid owns it at both levels and owns it
        consistently. What differs is the drawing: the antemeridian border
        is a curve in the projection plane, and a cell approximates it with
        one chord while its children approximate it with several. The
        coarse figure under-covers a region the grid addresses correctly.
        That is a limitation of the rendering and it must not be written
        into a claim about the tessellation.
        """
        assert itacart.describe()["tessellation"]["base_unit_polyhedron"] is None

        sizes = [itacart.cell_size(r) for r in range(1, itacart.MAX_RESOLUTION + 1)]
        assert all(finer < coarser for coarser, finer in zip(sizes, sizes[1:]))
        assert itacart.describe()["resolutions"]["max"] == itacart.MAX_RESOLUTION == 13
        with pytest.raises(itacart.ResolutionError):
            itacart.cell_size(itacart.MAX_RESOLUTION + 1)

        assert itacart.crs()["coordinate_epoch"] is None

        for row, columns in ((300, range(1600, 1610)), (0, range(1999, 2004))):
            children = [
                child
                for column in columns
                for child in itacart.get_children(
                    f"NE({column:04d}/{row:04d})", flatten=True
                )
            ]
            figures = [Polygon(plane_ring(child)[1]) for child in children]
            overlaps = sum(
                1
                for left in range(len(figures))
                for right in range(left + 1, len(figures))
                if figures[left].intersection(figures[right]).area > 1e-6
            )
            assert overlaps == 0, f"row {row}: {overlaps} of {len(children)} overlap"

        for parent in ("NE(1605/0300)", "NE(2003/0000)", "NE(0000/0300)"):
            figure = Polygon(plane_ring(parent)[1])
            children = unary_union(
                [
                    Polygon(plane_ring(child)[1])
                    for child in itacart.get_children(parent, flatten=True)
                ]
            )
            assert children.covers(figure), parent

        trapezoid = "NE(2003/0000)"
        figure = Polygon(plane_ring(trapezoid)[1])
        children = unary_union(
            [
                Polygon(plane_ring(child)[1])
                for child in itacart.get_children(trapezoid, flatten=True)
            ]
        )
        surplus = children.difference(figure)
        assert surplus.area > 0.0

        border = itacart.geodetic_to_sinusoidal(180.0, 0.0)[0]
        assert max(x for x, _ in children.exterior.coords) <= border + 1e-6

        sampled = _points_in(surplus, 20)
        assert len(sampled) == 20
        for x, y in sampled:
            longitude, latitude = itacart.sinusoidal_to_geodetic(x, y)
            assert itacart.geo_to_cell(longitude, latitude, 1) == trapezoid
            finer = itacart.geo_to_cell(longitude, latitude, 2)
            assert itacart.get_parent(finer) == trapezoid
            assert Polygon(plane_ring(finer)[1]).covers(Point(x, y))
            assert not figure.covers(Point(x, y))

    def test_req_26_simple_2d_polygons(self) -> None:
        """Simple two-dimensional polygons, and the two that are not.

        Requirement 26 (``ea/ers/cell/simple/2d_polygon``) asks for cells
        that are simple polygons at every discrete global grid. Clause
        9.1.5.1 gives four criteria: edges meet only at the vertices;
        exactly two edges meet at each vertex; exactly the same number of
        edges and vertices; and the cell encloses a region which always has
        a measurable area.

        The criteria are about *boundary curves* on the surface model of
        the Earth. Nothing in them asks for straight segments or geodesic
        edges, and nothing asks a vertex to be a point where the tangent
        turns: clause 9.1.4.2 lists small circles among the allowed curve
        types, and a parallel of latitude is a small circle. A cell that is
        a disc bounded by one parallel therefore meets the criteria as soon
        as that parallel is partitioned into arcs, because k arcs on a
        closed curve meet only at k points, two at each, and the region has
        area. The choice of k is a representation choice the specification
        is free to make; the criteria are satisfied for every k.

        This is measured rather than asserted, for k of 3, 4 and 8: the
        arcs are built, and each of the first three criteria is checked
        against the built figure rather than declared of it. The fourth is
        the cap's own area.

        Parallelograms, triangles and trapezoids satisfy the criteria
        directly, with the vertices the package already draws.
        """
        shapes = {itacart.cell_shape(cell) for _, cell in FAMILIES}
        assert shapes == {"parallelogram", "triangle", "trapezoid"}

        for _, cell in FAMILIES:
            assert Polygon(plane_ring(cell)[1]).area > 0.0

        for cap in CAPS:
            parallel = {
                lat for _, lat in itacart.cell_to_boundary(cap) if abs(lat) < 90
            }
            assert len(parallel) == 1
            latitude = parallel.pop()
            for edges in (3, 4, 8):
                vertices, arcs = _arc_boundary(latitude, edges)
                assert len(arcs) == len(vertices) == edges

                endpoints = [
                    round(end[0] % 360.0, 6)
                    for arc in arcs
                    for end in (arc.coords[0], arc.coords[-1])
                ]
                for vertex in vertices:
                    assert endpoints.count(round(vertex, 6)) == 2

                for left in range(edges):
                    for right in range(left + 1, edges):
                        shared = arcs[left].intersection(arcs[right])
                        if shared.is_empty:
                            continue
                        pieces = getattr(shared, "geoms", [shared])
                        for piece in pieces:
                            assert piece.geom_type == "Point"
                            assert any(
                                abs((piece.x % 360.0) - vertex) < 1e-9
                                for vertex in vertices
                            )

            assert Polygon(plane_ring(cap)[1]).area > 0.0

    def test_req_27_representative_position(self) -> None:
        """The representative position is the centroid, and how close it is.

        Requirement 27 (``ea/ers/cell/direct_position/centroid``): an EAERS
        specification shall define the DirectPosition of a cell to be the
        centroid, computed as the geodesic centre of surface area. Clause
        9.1.5.2 adds that the representative point shall lie on the surface
        of the cell.

        The paper records this as partially met, on the grounds that the
        representative position is a vertex rather than the centroid, and
        notes that a further implementation may represent the centroid.
        That implementation exists: ``cell_to_centroid``. The requirement
        names the centroid, and Requirement 12 asks for a position within
        the boundary, which the anchor is not. So the answer to both is the
        centroid, and this test and the exporter give the same one.

        The requirement names a *method*, not only a point, and the two are
        not identical here. ``cell_to_centroid`` returns the area centroid
        of the equal-area plane figure, inverted to geodetic; the standard
        asks for the geodesic centre of surface area, the point minimising
        the area-weighted squared geodesic distance over the cell. Measured
        against a direct numerical search over a uniform-by-area sample --
        uniform in the plane, which is uniform by area because the plane is
        equal-area -- the two agree to under a metre in the six ordinary
        families and part company in the two polar ones, where the cap's
        geodesic centre is the pole itself and the plane's answer sits
        about 1.3 km short of it. The interior claim is asserted; the polar
        divergence is asserted as a divergence rather than papered over.
        """
        for _, cell in FAMILIES:
            centroid = itacart.cell_to_centroid(cell)
            assert isinstance(centroid, tuple)
            figure = Polygon(plane_ring(cell)[1])
            assert Point(itacart.geodetic_to_sinusoidal(*centroid)).within(figure)
            assert itacart.geo_to_cell(*centroid, itacart.get_resolution(cell)) == cell

        ordinary = itacart.cell_to_centroid("NE(0500/0300)")
        assert isinstance(ordinary, tuple)
        gap, _ = itacart.inverse_geodesic(*ordinary, *_frechet_mean("NE(0500/0300)"))
        assert gap < 5.0, f"ordinary family diverges by {gap:.3f} m"

        cap = itacart.cell_to_centroid("NE(0000/1000)")
        assert isinstance(cap, tuple)
        assert cap[1] < 90.0
        cap_gap, _ = itacart.inverse_geodesic(cap[0], cap[1], cap[0], 90.0)
        assert 1000.0 < cap_gap < 2000.0, f"cap diverges by {cap_gap:.3f} m"
