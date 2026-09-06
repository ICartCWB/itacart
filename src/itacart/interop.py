"""Export to standard geospatial formats.

Satisfies OGC DGGS Core requirements 18 and 19 (interoperability
functions). The paper marks these as met by design; concrete exporters are
what turn that into a demonstrable claim.

Provenance: new; the notebook did this ad hoc with GeoPandas.

Three edge conventions meet in this module and they are not the same. A cell
is defined by straight edges in the *sinusoidal plane*, and that is the
convention its area is computed under. Geodesic paths are a second
convention, used wherever a distance on the ellipsoid is meant. GeoJSON
imposes a third: RFC 7946 section 3.1.1 defines the line between two
positions as a straight segment in longitude and latitude, and warns that it
may differ markedly from the path along the curved surface. **The exported
polygon is therefore not the cell.** The difference is a property of the
format, not a defect to repair, and it is measured in the test suite rather
than assumed.

Coordinates are never rounded. Six decimal places are roughly ten
centimetres, which is smaller than a cell only down to resolution 10; from
resolution 11 the whole cell is narrower than one rounding step and any
rounding collapses it into a degenerate polygon that is still syntactically
valid GeoJSON. Emitting full precision costs bytes and never lies.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any, cast

from shapely.geometry import MultiPolygon, Polygon, box, mapping, shape
from shapely.geometry.base import BaseGeometry
from shapely.geometry.polygon import orient
from shapely.ops import unary_union

from .boundary import cell_shape, extension_zone, is_valid_cell
from .cells import cell_to_boundary
from .constants import WGS84_A, WGS84_E2
from .exceptions import GeometryError, NonExistentCellError
from .geometry import Containment, polyfill
from .index import decompose
from .resolutions import effective_cell_area, get_resolution, nominal_cell_area

if TYPE_CHECKING:
    import geopandas

__all__ = [
    "cells_to_geojson",
    "cell_to_wkt",
    "cells_to_wkt",
    "to_geodataframe",
    "from_geodataframe",
    "from_geojson",
    "recover_from_geojson",
]

#: Default target length of a densified segment along the base parallel of a
#: polar cap, in metres. Matches the package-wide densification default.
_DEFAULT_MAX_SEGMENT_M = 1000.0

#: Property names written by :func:`cells_to_geojson` when
#: ``include_metadata`` is true. These are public API: once emitted, renaming
#: one breaks every consumer that reads it.
INDEX_PROPERTY = "itacart_index"
RESOLUTION_PROPERTY = "itacart_resolution"
SHAPE_PROPERTY = "itacart_shape"
NOMINAL_AREA_PROPERTY = "itacart_nominal_area_m2"
EFFECTIVE_AREA_PROPERTY = "itacart_effective_area_m2"
EXTENSION_ZONE_PROPERTY = "itacart_extension_zone"


def _require_cell(cell: str) -> None:
    """Refuse an index that names no cell.

    Six public functions return geometry for a cell that does not exist. The
    exporters do not inherit that: a caller who asks for GeoJSON of an
    unoccupied address gets a refusal rather than a plausible polygon, and
    the refusal is a package exception so that one ``except ITACaRTError``
    still guards the pipeline.
    """
    if not is_valid_cell(cell):
        raise NonExistentCellError(f"index names no cell: {cell}")


def _holds_pole(ring: list[tuple[float, float]]) -> bool:
    """Whether a boundary ring reaches the pole.

    This is the operational test for a polar cap, and it is deliberately not
    a longitude-jump test. Minus 180 and plus 180 are the same physical
    meridian, and at the pole longitude is undefined, so a detector built on
    the size of a longitude step cannot tell a cap from a feature that
    genuinely occupies both sides of the seam.

    Nor is it a proximity test, and that distinction is not cosmetic. This
    was a comparison against 89.999 degrees, justified by the measurement
    that no cell which is not a cap reaches beyond 89.9824. That
    measurement was true and its scope was resolution one, which was never
    said. Refined, ordinary cells walk toward the pole without limit: three
    of the cap's grandchildren pass 89.999 at resolution four and three
    more at resolution five, none of them caps, and each was handed to the
    cap routine and exported as a band around the whole parallel. The
    inflation ran from ten to seventeen times the cell's own area, and
    three distinct cells came out as one identical footprint. No fixed
    threshold below ninety survives refinement, because the cells keep
    coming.

    The pole is not approximated anywhere in this package, so it does not
    have to be approximated here. ``itacart.boundary.to_geodetic`` answers
    ``copysign(90.0, y)`` as a literal whenever a plane vertex reaches the
    meridian quadrant, which is exactly the construction that produces a
    cap. Asking for that literal asks the same question the grid answers,
    rather than a question that resembles it.
    """
    return any(abs(lat) >= 90.0 for _, lat in ring)


def _parallel_circumference_m(latitude: float) -> float:
    """Length of the full parallel at ``latitude`` on the WGS84 ellipsoid."""
    phi = math.radians(latitude)
    prime_vertical = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(phi) ** 2)
    return 2.0 * math.pi * prime_vertical * math.cos(phi)


def _polar_cap_polygon(
    ring: list[tuple[float, float]],
    max_segment_m: float = _DEFAULT_MAX_SEGMENT_M,
) -> Polygon:
    """Build the interoperable polygon of a polar cap.

    The three-vertex ring that :func:`itacart.cell_to_boundary` returns for a
    cap is a coordinate-domain figure, not a triangle on the ellipsoid. Its
    two meridional sides at minus 180 and plus 180 are two copies of the same
    physical seam, its base is a complete parallel, and its apex is a point
    at which every longitude names the same place. A GeoJSON consumer has no
    way to know that, so the exporter converts the cap into a figure that
    means the same thing under RFC 7946 rules: the base parallel densified
    across the whole longitude range, the surface opened along the seam, and
    the apex carried as a degenerate edge at the pole.

    The result is one polygon with one part. A cap does not occupy two sides
    of the seam; it is produced by the seam, so cutting it would invent a
    division that the surface does not have.
    """
    latitudes = [lat for _, lat in ring]
    northern = max(latitudes, key=abs) > 0.0
    base = min(latitudes) if northern else max(latitudes)
    circumference = _parallel_circumference_m(base)
    steps = max(4, math.ceil(circumference / max_segment_m))
    lons = [-180.0 + 360.0 * step / steps for step in range(steps + 1)]
    if northern:
        parallel = [(lon, base) for lon in lons]
        polar_edge = [(180.0, 90.0), (-180.0, 90.0)]
    else:
        parallel = [(lon, -90.0) for lon in lons]
        polar_edge = [(180.0, base), (-180.0, base)]
    positions = parallel + polar_edge
    return Polygon(positions + [positions[0]])


def _unwrap_longitudes(
    ring: list[tuple[float, float]],
) -> list[tuple[float, float]]:
    """Make a ring continuous in longitude by lifting it off the seam."""
    unwrapped = [ring[0]]
    for lon, lat in ring[1:]:
        previous = unwrapped[-1][0]
        while lon - previous > 180.0:
            lon -= 360.0
        while previous - lon > 180.0:
            lon += 360.0
        unwrapped.append((lon, lat))
    return unwrapped


def _cut_at_antimeridian(geometry: BaseGeometry) -> BaseGeometry:
    """Cut an arbitrary geometry that occupies both sides of the seam.

    RFC 7946 section 3.1.9 says geometry crossing the antimeridian should be
    split so that no part crosses it. This is the general rule and it applies
    to arbitrary input geometry on the way out. It is not what handles a
    polar cap: see :func:`_polar_cap_polygon` for why those are a different
    problem with a different answer.

    No ITACaRT cell reaches this path. Measured across every column of the
    equatorial, mid-latitude and polar rows in all four quadrants, the widest
    ring spans exactly 180 degrees and the last column of every row stops on
    the seam without passing it. The routine exists for geometry that did not
    come from the grid.
    """
    if not isinstance(geometry, Polygon):
        raise GeometryError(
            "antimeridian cutting supports Polygon input only, "
            f"got {geometry.geom_type}"
        )
    exterior = list(geometry.exterior.coords)
    lifted = Polygon(_unwrap_longitudes(exterior))
    west = lifted.intersection(box(-540.0, -90.0, 180.0, 90.0))
    east = lifted.intersection(box(180.0, -90.0, 540.0, 90.0))
    parts: list[Polygon] = []
    for piece, offset in ((west, 0.0), (east, -360.0)):
        if piece.is_empty:
            continue
        for polygon in getattr(piece, "geoms", [piece]):
            shifted = [(lon + offset, lat) for lon, lat in polygon.exterior.coords]
            parts.append(Polygon(shifted))
    if len(parts) == 1:
        return parts[0]
    return MultiPolygon(parts)


def _cell_polygon(cell: str) -> Polygon:
    """The RFC 7946 polygon of one atomic cell, closed and right-handed.

    Rings come out explicitly closed, first position identical to last, as
    section 3.1.6 requires: a triangle in four positions, a parallelogram in
    five. Orientation is normalised rather than trusted. Measured,
    :func:`itacart.cell_to_boundary` already returns counterclockwise rings
    in all four quadrants, but section 3.1.6 makes right-handedness a
    requirement on production, so the exporter enforces it instead of
    relying on a measurement of today's behaviour. That is a deliberate
    divergence from the binary serialisation, which preserves whatever
    orientation it was given because it exists to reproduce identity, while
    this module exists to interoperate.
    """
    ring = cast("list[tuple[float, float]]", cell_to_boundary(cell, close=True))
    if _holds_pole(ring):
        return orient(_polar_cap_polygon(ring), sign=1.0)
    return orient(Polygon(ring), sign=1.0)


def _cell_metadata(cell: str) -> dict[str, Any]:
    """The declared metadata block of one cell."""
    return {
        INDEX_PROPERTY: cell,
        RESOLUTION_PROPERTY: get_resolution(cell),
        SHAPE_PROPERTY: str(cell_shape(cell)),
        NOMINAL_AREA_PROPERTY: nominal_cell_area(get_resolution(cell)),
        EFFECTIVE_AREA_PROPERTY: effective_cell_area(cell),
        EXTENSION_ZONE_PROPERTY: extension_zone(cell),
    }


def cells_to_geojson(
    index: str,
    properties: dict[str, Any] | None = None,
    include_metadata: bool = True,
) -> dict[str, Any]:
    """Export a compositional index as a GeoJSON FeatureCollection.

    One cell becomes one Feature, always, and a Feature carries a Polygon or
    a MultiPolygon rather than being split into two Features. The index goes
    into the Feature ``id`` member, which RFC 7946 section 3.2 reserves for
    an identifier already in common use, and is repeated in ``properties``
    when metadata is asked for, because some consumers drop ``id`` on
    import. Coordinates are WGS84 decimal degrees and nothing else: section 4
    fixes the coordinate reference system and records that alternatives were
    removed from the specification for interoperability reasons.

    Args:
        index: Compositional index; every cell it names becomes a Feature.
        properties: Extra members merged into every Feature's ``properties``.
        include_metadata: Whether to add the ITACaRT metadata block.

    Returns:
        A FeatureCollection dictionary.

    Raises:
        NonExistentCellError: If a named cell does not exist.
    """
    features: list[dict[str, Any]] = []
    for cell in decompose(index):
        _require_cell(cell)
        members: dict[str, Any] = {}
        if include_metadata:
            members.update(_cell_metadata(cell))
        if properties is not None:
            members.update(properties)
        features.append(
            {
                "type": "Feature",
                "id": cell,
                "geometry": mapping(_cell_polygon(cell)),
                "properties": members,
            }
        )
    return {"type": "FeatureCollection", "features": features}


def cell_to_wkt(cell: str) -> str | list[str]:
    """Export a cell as WKT, one string per cell.

    Follows the convention of the six functions that already answer a list
    when handed a compositional index naming several cells: the list is in
    the order :func:`itacart.decompose` returns.

    Args:
        cell: Atomic cell or compositional index.

    Returns:
        A WKT string, or a list of them in decompose order.

    Raises:
        NonExistentCellError: If a named cell does not exist.
    """
    cells = decompose(cell)
    for atom in cells:
        _require_cell(atom)
    if len(cells) == 1:
        return str(_cell_polygon(cells[0]).wkt)
    return [str(_cell_polygon(atom).wkt) for atom in cells]


def cells_to_wkt(index: str, dissolve: bool = False) -> str:
    """Export a whole index as a single WKT geometry.

    With ``dissolve`` false the result is a GEOMETRYCOLLECTION holding one
    polygon per cell. With ``dissolve`` true adjacent cells are merged, and
    the two seams behave differently: along the prime meridian the shared
    cell is a single cell claimed geometrically from both sides, so the union
    closes without a gap; across the antimeridian, cells that are neighbours
    on the ellipsoid are not neighbours in longitude and latitude, so the
    union leaves them apart and the result stays multi-part.

    Args:
        index: Compositional index.
        dissolve: Whether to merge adjacent cells into one geometry.

    Returns:
        A WKT string.

    Raises:
        NonExistentCellError: If a named cell does not exist.
    """
    cells = decompose(index)
    for atom in cells:
        _require_cell(atom)
    polygons = [_cell_polygon(atom) for atom in cells]
    if dissolve:
        return str(unary_union(polygons).wkt)
    return str(
        "GEOMETRYCOLLECTION (" + ", ".join(polygon.wkt for polygon in polygons) + ")"
    )


def geometry_to_geojson(geometry: BaseGeometry) -> dict[str, Any]:
    """Serialise arbitrary geometry, cutting it at the antimeridian.

    This is the export path for geometry that did not come from the grid.
    Anything that occupies both sides of the seam is split per RFC 7946
    section 3.1.9 so that no part crosses it. Polar caps are not handled
    here; they are a cell concern and have their own routine.

    Args:
        geometry: Any Shapely polygon in WGS84 degrees.

    Returns:
        A GeoJSON geometry dictionary.

    Raises:
        GeometryError: If the geometry type is not supported.
    """
    if not isinstance(geometry, Polygon):
        raise GeometryError(f"unsupported geometry type: {geometry.geom_type}")
    lons = [lon for lon, _ in geometry.exterior.coords]
    if max(lons) - min(lons) > 180.0:
        return dict(mapping(_cut_at_antimeridian(geometry)))
    return dict(mapping(geometry))


def recover_from_geojson(obj: dict[str, Any]) -> list[str]:
    """Read back the indices that :func:`cells_to_geojson` wrote.

    This is the recovery path and it is exact: the index is in the file, in
    the Feature ``id`` member, so nothing is inferred from geometry and
    nothing is lost. It is the inverse that makes the exporter checkable
    rather than merely plausible.

    Args:
        obj: A FeatureCollection produced by :func:`cells_to_geojson`.

    Returns:
        One index per Feature, in file order.

    Raises:
        GeometryError: If a Feature carries no recoverable index.
    """
    indices: list[str] = []
    for feature in obj.get("features", []):
        index = feature.get("id")
        if index is None:
            index = feature.get("properties", {}).get(INDEX_PROPERTY)
        if index is None:
            raise GeometryError("feature carries no ITACaRT index in id or properties")
        indices.append(str(index))
    return indices


def from_geojson(
    obj: dict[str, Any],
    resolution: int,
    containment: str = "center",
) -> list[str]:
    """Fill every geometry of a GeoJSON object into a compositional index.

    This is the filling path, and it is lossy by construction: the cells that
    cover a shape are not the shape. It is named apart from
    :func:`recover_from_geojson` on purpose, because the two answer different
    questions and sharing a name with a switching argument would hide which
    one a caller got.

    There is no ``crs`` argument. RFC 7946 section 4 fixes GeoJSON to WGS84
    and records that alternative coordinate reference systems were removed
    from the specification, so accepting one here would be accepting
    something the format does not have.

    Args:
        obj: A GeoJSON Feature, FeatureCollection or bare geometry.
        resolution: Target resolution level.
        containment: Predicate passed through to
            :func:`itacart.geometry.polyfill`.

    Returns:
        One compositional index per input geometry, in file order.
    """
    if obj.get("type") == "FeatureCollection":
        geometries = [feature["geometry"] for feature in obj.get("features", [])]
    elif obj.get("type") == "Feature":
        geometries = [obj["geometry"]]
    else:
        geometries = [obj]
    return [
        polyfill(
            shape(geometry),
            resolution,
            containment=cast("Containment", containment),
        )
        for geometry in geometries
    ]


def to_geodataframe(index: str, crs: str = "EPSG:4326") -> "geopandas.GeoDataFrame":
    """Export a compositional index as a GeoDataFrame.

    Requires the ``geo`` extra. Installing it brings ``pyproj`` back
    transitively, which does not undo the removal of ``pyproj`` from the
    package: the package itself still computes its geodesy directly, and this
    is an optional extra pulling in a reprojection library for the caller's
    benefit.

    Args:
        index: Compositional index.
        crs: Target coordinate reference system; reprojected from EPSG:4326.

    Returns:
        One row per cell, carrying the declared metadata columns.

    The index is checked before the extra is imported, so an address
    naming no cell is refused with the package's own exception whether or
    not GeoPandas is installed. The other order made the answer depend on
    the environment: with the extra present the caller was told the cell
    does not exist, and without it the caller was told to install
    GeoPandas, which was true and was not the problem.

    Raises:
        NonExistentCellError: If a named cell does not exist.
        ImportError: If GeoPandas is not installed.
    """
    collection = cells_to_geojson(index, include_metadata=True)
    geopandas_module = _import_geopandas()
    frame = geopandas_module.GeoDataFrame.from_features(
        collection["features"], crs="EPSG:4326"
    )
    if crs != "EPSG:4326":
        frame = frame.to_crs(crs)
    return frame


def from_geodataframe(
    gdf: "geopandas.GeoDataFrame",
    resolution: int,
    containment: str = "center",
) -> list[str]:
    """Fill every geometry of a GeoDataFrame into a compositional index.

    Requires the ``geo`` extra. This is a thin adapter: it reprojects to
    EPSG:4326 and delegates to :func:`from_geojson`, so the filling logic
    exists in exactly one place, and that place is the path that needs no
    optional extra and is therefore testable in any environment.

    Args:
        gdf: Input frame in any CRS; reprojected to EPSG:4326 as needed.
        resolution: Target resolution level, 1 to 13.
        containment: Predicate passed through to
            :func:`itacart.geometry.polyfill`.

    Returns:
        One compositional index per row, in row order.

    Raises:
        ImportError: If GeoPandas is not installed.
    """
    _import_geopandas()
    frame = gdf if gdf.crs is None else gdf.to_crs("EPSG:4326")
    return from_geojson(
        dict(frame.__geo_interface__), resolution, containment=containment
    )


def _import_geopandas() -> Any:
    """Import GeoPandas, or say plainly which extra is missing."""
    try:
        import geopandas as geopandas_module
    except ImportError as exc:  # pragma: no cover - exercised by injection
        raise ImportError(
            "to_geodataframe and from_geodataframe need the 'geo' extra: "
            "pip install itacart[geo]"
        ) from exc
    return geopandas_module
