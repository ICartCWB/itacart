"""Cells as vertices: the inverse of ``vertex_to_cell``, declared and enumerated.

A vertex written as a cell is a position together with how well it is
known. Reading it back asks one question -- which point of the cell stands
for the vertex -- and this module holds the package to a single answer on
every route, the centroid, before it measures what that answer costs.

**The population is enumerated, not drawn.** Positions sit on a lattice of
latitudes and longitudes in the four quadrants, at five resolutions from
10 km to 1 cm, plus the two extension zones and a lattice over the polar
caps at every resolution. The families this grid breaks rules in -- the
prime-meridian column, the border-absorbing column, the southern equator
row, the polar row, the caps, the extensions -- are counted, and the counts
are asserted, so that a family cannot drop out of the population and leave
its tests passing over nothing.

**Distances are geodesic and the farthest point is a vertex.** A cell's
edges are straight on the sinusoidal plane. Densifying every edge six-fold
on the plane, inverting and measuring moved no reach in the enumerated
population, so the reach of a cell -- the distance from its centroid to
its farthest point -- is taken over its vertices.

**Longitude is the shear, not latitude.** The shape of an ordinary cell
on the ground is set by ``u = |lon * sin(lat)|``, lon in radians, so a
latitude sweep at one longitude covers a fraction of the range. The
lattice spans ``u`` from 0 to pi.
"""

from __future__ import annotations

import math
from functools import lru_cache

import pytest
from shapely.geometry import Point, Polygon

import itacart
from itacart.boundary import plane_ring
from itacart.geodesy import (
    geodetic_to_sinusoidal,
    inverse_geodesic,
    sinusoidal_to_geodetic,
)

LATITUDES = (1e-7, 10.0, 20.0, 30.0, 45.0, 60.0, 70.0, 80.0, 85.0, 89.0, 89.9)
LONGITUDES = (1e-7, 0.5, 10.0, 30.0, 60.0, 90.0, 120.0, 150.0, 170.0, 179.0, 179.999)
RESOLUTIONS = (1, 2, 5, 9, 13)
SIGNS = ((1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0))

#: Fiji and Chukotka, each written both ways the package reads a position
#: past the antemeridian: west of it, and in the eastern spelling.
EXTENSION_POSITIONS = (
    (-179.5, -16.0),
    (-178.5, -18.5),
    (180.5, -21.0),
    (-179.0, 65.0),
    (-175.0, 68.0),
    (186.0, 71.0),
)

#: The caps begin at about 89.9823 degrees; the lattice runs to the pole.
CAP_LATITUDES = (89.9824, 89.9868, 89.9912, 89.9956, 90.0)
CAP_LONGITUDES = tuple(-180.0 + 45.0 * k for k in range(9))

#: Latitude below which the shear law holds within ORDINARY_LAW_TOLERANCE.
LAW_LATITUDE_LIMIT = 85.0


def _res1_row(cell: str) -> int:
    return int(itacart.split_components(cell)[1].split("/")[1])


def _family(cell: str) -> str:
    row = _res1_row(cell)
    column = int(itacart.split_components(cell)[1].split("/")[0])
    if row == 1000:
        return "cap"
    if row == 999:
        return "polar row"
    if itacart.is_extension_cell(cell):
        return "extension"
    if itacart.absorbs_border(cell):
        return "absorbing"
    if itacart.cell_shape(cell) == "triangle":
        return "meridian triangle"
    if itacart.is_quadrant_boundary_cell(cell):
        return "southern equator row"
    return "ordinary" if column > 0 else "meridian column"


@lru_cache(maxsize=None)
def _population() -> tuple[str, ...]:
    positions = [
        (sx * lon, sy * lat, resolution)
        for resolution in RESOLUTIONS
        for lat in LATITUDES
        for lon in LONGITUDES
        for sx, sy in SIGNS
    ]
    positions += [
        (lon, lat, resolution)
        for resolution in RESOLUTIONS
        for lon, lat in EXTENSION_POSITIONS
    ]
    positions += [
        (lon, sign * lat, resolution)
        for resolution in range(1, 14)
        for lat in CAP_LATITUDES
        for lon in CAP_LONGITUDES
        for sign in (1.0, -1.0)
    ]
    cells = {itacart.geo_to_cell(lon, lat, res) for lon, lat, res in positions}
    return tuple(sorted(cells))


def _distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    if a == b:
        return 0.0
    return inverse_geodesic(a[0], a[1], b[0], b[1])[0]


def _reach(cell: str, point: tuple[float, float]) -> float:
    """Farthest vertex of the cell from ``point``, in metres."""
    ring = itacart.cell_to_boundary(cell)
    return max(_distance(point, (lon, lat)) for lon, lat in ring)


def _shear(cell: str) -> float:
    lon, lat = itacart.cell_to_centroid(cell)
    return abs(math.radians(lon) * math.sin(math.radians(lat)))


def shear_law(u: float) -> float:
    """Reach of an ordinary parallelogram over its side, to first order.

    On the ground the cell's plane vertices ``(0, 0)``, ``(l, 0)``,
    ``(0, l)`` and ``(-l, l)`` sit at ``(0, 0)``, ``(l, 0)``, ``(u l, l)``
    and ``((u - 1) l, l)``, because a plane displacement ``(dx, dy)`` is a
    ground displacement ``(dx + u dy, dy)``. The centroid is at
    ``(u l / 2, l / 2)``, and the farthest vertex gives the law.
    """
    return max(0.5 * math.sqrt(1.0 + u * u), math.sqrt((1.0 - u / 2.0) ** 2 + 0.25))


def anchor_law(u: float) -> float:
    """Reach of the same parallelogram measured from its anchor instead."""
    return max(math.sqrt(1.0 + (1.0 - u) ** 2), math.sqrt(1.0 + u * u))


# --------------------------------------------------------------------------
# The population is what it says it is
# --------------------------------------------------------------------------


def test_the_population_names_every_family_and_says_how_many() -> None:
    """Every family is present, in every quadrant it exists in, and counted."""
    counts: dict[str, int] = {}
    quadrants: dict[str, set[str]] = {}
    for cell in _population():
        family = _family(cell)
        counts[family] = counts.get(family, 0) + 1
        quadrants.setdefault(family, set()).add(itacart.quadrant_of(cell))
    assert counts == {
        "absorbing": 62,
        "cap": 520,
        "extension": 87,
        "meridian column": 56,
        "meridian triangle": 86,
        "ordinary": 1744,
        "polar row": 338,
        "southern equator row": 76,
    }
    for family in ("absorbing", "ordinary", "polar row"):
        assert quadrants[family] == {"NE", "NW", "SE", "SW"}, family
    assert quadrants["cap"] == {"NE", "SE"}
    assert all(itacart.is_valid_cell(cell) for cell in _population())


# --------------------------------------------------------------------------
# One declared point, whatever the route
# --------------------------------------------------------------------------


def test_every_single_point_route_answers_with_the_centroid() -> None:
    """The rebuilt vertex, the H3 alias, the engine and the declaration agree.

    The answer is compared for equality, not closeness: the routes are
    required to be one computation, and a second computation that agreed
    to a nanometre would still be a second answer.
    """
    engine = itacart.ITACaRT()
    assert itacart.describe()["cell_geometry"]["representative_position"] == (
        "centroid"
    )
    for cell in _population():
        centroid = itacart.cell_to_centroid(cell)
        rebuilt = itacart.cells_to_geometry([cell], "Point")
        assert (rebuilt.x, rebuilt.y) == centroid, cell
        assert itacart.cell_to_latlng(cell) == (centroid[1], centroid[0]), cell
        assert engine.cell_to_centroid(cell) == centroid, cell


def test_a_rebuilt_ring_is_the_centroids_of_its_cells_in_order() -> None:
    ring = [(10.0, 45.0), (10.001, 45.0), (10.001, 45.001), (10.0, 45.001)]
    cells = itacart.vertex_to_cell(Polygon(ring), 13)
    rebuilt = itacart.cells_to_geometry(cells, "Polygon")
    assert list(rebuilt.exterior.coords)[:-1] == [
        itacart.cell_to_centroid(cell) for cell in cells
    ]


def test_the_anchor_is_on_the_boundary_and_the_rebuilt_vertex_is_not() -> None:
    """Why the anchor is not offered: a vertex is not a point within the cell."""
    for cell in ("NE(0500/0300(3(C2)))", "SW(0476/0260(4(A2(2(E3)))))"):
        ring = plane_ring(cell)[1]
        polygon = Polygon([(x - ring[0][0], y - ring[0][1]) for x, y in ring])
        anchor = geodetic_to_sinusoidal(*itacart.cell_to_anchor(cell))
        rebuilt = itacart.cells_to_geometry([cell], "Point")
        centre = geodetic_to_sinusoidal(rebuilt.x, rebuilt.y)
        shift = (-ring[0][0], -ring[0][1])
        on_edge = Point(anchor[0] + shift[0], anchor[1] + shift[1])
        inside = Point(centre[0] + shift[0], centre[1] + shift[1])
        assert polygon.exterior.distance(on_edge) < 1e-6
        assert polygon.contains(inside)


# --------------------------------------------------------------------------
# The round trip, enumerated
# --------------------------------------------------------------------------


def test_every_rebuilt_vertex_quantizes_back_to_the_cell_that_encoded_it() -> None:
    """Over the whole population, extension zones and caps included.

    The two families this could not claim before are the two the
    enumeration found: a cap position quantized to a spelling that names
    no cell, and a rebuilt extension vertex, written past 180 degrees as
    the package writes it, was refused by the quantizer.
    """
    for cell in _population():
        rebuilt = itacart.cells_to_geometry([cell], "Point")
        resolution = itacart.get_resolution(cell)
        assert itacart.geo_to_cell(rebuilt.x, rebuilt.y, resolution) == cell, cell


def _probes(cell: str, fraction: float) -> list[tuple[float, float]]:
    """Positions just inside each vertex, moved toward the centroid on the plane."""
    centre = geodetic_to_sinusoidal(*itacart.cell_to_centroid(cell))
    probes = []
    for lon, lat in itacart.cell_to_boundary(cell):
        x, y = geodetic_to_sinusoidal(lon, lat)
        probes.append(
            sinusoidal_to_geodetic(
                x + fraction * (centre[0] - x), y + fraction * (centre[1] - y)
            )
        )
    return probes


def test_a_vertex_comes_back_within_its_cells_reach_and_near_it_at_a_corner() -> None:
    """The reach is a bound the round trip reaches, not a loose ceiling.

    Every cell of the ordinary and meridian families is probed just inside
    each of its vertices, one per cent of the way to the centroid. Each
    probe goes through ``vertex_to_cell`` and ``cells_to_geometry``; it must
    come back to the same cell, within the cell's reach, and the probe at
    the farthest corner must come back at least 98 per cent of the reach
    away, which is as close as a probe one per cent inside can get.
    """
    probed = 0
    for cell in _population():
        if _family(cell) not in ("ordinary", "meridian column", "meridian triangle"):
            continue
        resolution = itacart.get_resolution(cell)
        centroid = itacart.cell_to_centroid(cell)
        reach = _reach(cell, centroid)
        farthest = 0.0
        for lon, lat in _probes(cell, 0.01):
            cells = itacart.vertex_to_cell(Point(lon, lat), resolution)
            assert cells == [cell], (cell, lon, lat, cells)
            back = itacart.cells_to_geometry(cells, "Point")
            error = _distance((lon, lat), (back.x, back.y))
            assert error <= reach, cell
            farthest = max(farthest, error)
        assert farthest >= 0.98 * reach, cell
        probed += 1
    assert probed == 1886


def test_the_reach_of_an_ordinary_cell_follows_the_shear_law() -> None:
    """Measured reach over side against the closed form, below 85 degrees.

    The law is first order in the side. Closer to the pole a finite cell
    spans enough of the meridian's convergence to exceed it, and the
    population above the limit is asserted to do so, so that the limit is
    earned rather than conservative.
    """
    below = []
    above = []
    for cell in _population():
        if _family(cell) != "ordinary":
            continue
        ratio = _reach(cell, itacart.cell_to_centroid(cell)) / itacart.cell_size(
            itacart.get_resolution(cell)
        )
        law = shear_law(_shear(cell))
        lat = abs(itacart.cell_to_centroid(cell)[1])
        (below if lat < LAW_LATITUDE_LIMIT else above).append(ratio / law)
    assert (len(below), len(above)) == (1456, 288)
    assert min(below) > 1.0 - 1e-3
    assert max(below) < 1.0 + 3e-3
    assert 1.0 + 3e-3 < max(above) < 1.0 + 2.1e-2


def test_the_shear_law_spans_its_range_in_the_population() -> None:
    """Instrument check: the lattice reaches both ends of the law."""
    shears = [_shear(cell) for cell in _population() if _family(cell) == "ordinary"]
    assert min(shears) < 1e-3
    assert max(shears) > 3.0
    assert shear_law(0.0) == pytest.approx(math.sqrt(5.0) / 2.0)
    assert shear_law(1.0) == pytest.approx(math.sqrt(2.0) / 2.0)
    assert shear_law(math.pi) == pytest.approx(math.sqrt(1.0 + math.pi**2) / 2.0)


def test_a_meridian_triangle_reaches_root_ten_over_three_sides() -> None:
    """The triangle has base two sides and height one; its centroid is a third up.

    A finite triangle departs from the value for the reason the
    parallelogram departs from its law, its side not being small against
    the distance to the pole: within 0.3 per cent below 85 degrees, and up
    to 1.4 per cent above.
    """
    triangles = [c for c in _population() if _family(c) == "meridian triangle"]
    ratios = {"below": [], "above": []}
    for cell in triangles:
        side = itacart.cell_size(itacart.get_resolution(cell))
        reach = _reach(cell, itacart.cell_to_centroid(cell)) / side
        lat = abs(itacart.cell_to_centroid(cell)[1])
        band = "below" if lat < LAW_LATITUDE_LIMIT else "above"
        ratios[band].append(reach / (math.sqrt(10.0) / 3.0))
    assert (len(ratios["below"]), len(ratios["above"])) == (70, 16)
    assert all(abs(ratio - 1.0) < 3e-3 for ratio in ratios["below"])
    assert all(abs(ratio - 1.0) < 1.5e-2 for ratio in ratios["above"])
    assert max(abs(ratio - 1.0) for ratio in ratios["above"]) > 3e-3


@pytest.mark.slow
def test_the_border_families_are_not_bounded_by_the_law_and_are_measured() -> None:
    """The named exemption, earned and pinned.

    Population: the last existing column of every resolution-1 row in the
    four quadrants, and the border, polar-row, cap and extension cells of
    the enumerated population above. Some exceed the shear law, which is what makes the
    exemption necessary; the largest reach is pinned by value, so that a
    change in either direction is seen.
    """
    side = itacart.cell_size(1)
    border = set(c for c in _population() if _family(c) in _BORDER_FAMILIES)
    for quadrant in ("NE", "NW", "SE", "SW"):
        for row in range(1001):
            last = itacart.last_lattice_column(quadrant, row, side)
            cell = f"{quadrant}({max(last, 0):04d}/{row:04d})"
            if itacart.is_valid_cell(cell) and _family(cell) in _BORDER_FAMILIES:
                border.add(cell)
    exceeding = 0
    largest = 0.0
    for cell in border:
        ratio = _reach(cell, itacart.cell_to_centroid(cell)) / itacart.cell_size(
            itacart.get_resolution(cell)
        )
        largest = max(largest, ratio)
        if ratio > shear_law(_shear(cell)) * (1.0 + 5e-3):
            exceeding += 1
    assert (len(border), exceeding, round(largest, 4)) == (4967, 1969, 2.0532)
    assert largest < 2.06


_BORDER_FAMILIES = ("absorbing", "polar row", "cap", "extension")


# --------------------------------------------------------------------------
# The single point of the opening measurement is a position, not a law
# --------------------------------------------------------------------------


SAO_PAULO = (-46.633, -23.55)


def test_the_worst_case_favours_the_centroid_in_every_ordinary_cell() -> None:
    """Measured against the anchor, the centroid's reach is never larger.

    The ratio follows the two laws: 1.26 on the prime meridian, a minimum
    of ``sqrt(5) - 1`` between, and exactly 2 from ``u = 1``, where the
    anchor becomes an end of the long diagonal.
    """
    ratios = []
    for cell in _population():
        if _family(cell) != "ordinary":
            continue
        centroid = _reach(cell, itacart.cell_to_centroid(cell))
        anchor = _reach(cell, itacart.cell_to_anchor(cell))
        ratios.append(anchor / centroid)
    floor = min(anchor_law(k / 1000.0) / shear_law(k / 1000.0) for k in range(3200))
    assert floor == pytest.approx(math.sqrt(5.0) - 1.0, abs=1e-6)
    assert min(ratios) > floor * (1.0 - 1e-3)
    assert max(ratios) < 2.0 + 1e-3
    assert anchor_law(0.0) / shear_law(0.0) == pytest.approx(2.0 * math.sqrt(0.4))
    assert anchor_law(2.0) / shear_law(2.0) == pytest.approx(2.0)


def test_the_anchor_wins_a_single_point_only_inside_its_own_corner() -> None:
    """The opening table's reversal, reproduced and then taken apart.

    At resolutions 5 and 9 the São Paulo point is closer to the centroid;
    at 11 and 13 to the anchor. Moving the same point to the far corner of
    its resolution-13 cell, without changing the resolution, reverses it:
    what decided the table was where the point sits in its cell.
    """
    lon, lat = SAO_PAULO
    winners = {}
    for resolution in (5, 9, 11, 13):
        cell = itacart.geo_to_cell(lon, lat, resolution)
        to_anchor = _distance((lon, lat), itacart.cell_to_anchor(cell))
        to_centroid = _distance((lon, lat), itacart.cell_to_centroid(cell))
        winners[resolution] = "anchor" if to_anchor < to_centroid else "centroid"
    assert winners == {5: "centroid", 9: "centroid", 11: "anchor", 13: "anchor"}

    cell = itacart.geo_to_cell(lon, lat, 13)
    probes = _probes(cell, 0.01)
    anchor = itacart.cell_to_anchor(cell)
    far = max(probes, key=lambda probe: _distance(probe, anchor))
    assert itacart.geo_to_cell(far[0], far[1], 13) == cell
    assert _distance(far, itacart.cell_to_centroid(cell)) < _distance(far, anchor)


# --------------------------------------------------------------------------
# Precision is the index's to declare
# --------------------------------------------------------------------------


@pytest.mark.parametrize("lat", [-23.55, 0.5, 45.0, 70.0])
@pytest.mark.parametrize("sx, sy", SIGNS)
def test_resolution_13_neighbours_survive_only_at_full_precision(
    lat: float, sx: float, sy: float
) -> None:
    """Six decimal places merge neighbours; the rebuilt coordinates do not.

    Each disk is a resolution-13 cell and its ring of neighbours. Rebuilt,
    every vertex is distinct and quantizes back to its own cell. Rounded to
    six decimal places, two neighbours become one point, and requantizing
    the rounded coordinates no longer returns the cells.
    """
    centre = itacart.geo_to_cell(sx * 46.633, sy * abs(lat), 13)
    disk = sorted(itacart.grid_disk(centre, 1))
    assert len(disk) == 9
    rebuilt = itacart.cells_to_geometry(disk, "MultiPoint")
    points = [(point.x, point.y) for point in rebuilt.geoms]
    assert len(set(points)) == len(points)
    assert [itacart.geo_to_cell(x, y, 13) for x, y in points] == disk

    rounded = [(round(x, 6), round(y, 6)) for x, y in points]
    assert len(set(rounded)) < len(rounded)
    assert [itacart.geo_to_cell(x, y, 13) for x, y in rounded] != disk
