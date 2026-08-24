"""Export to standard geospatial formats.

Satisfies OGC DGGS Core requirements 18 and 19 (interoperability
functions). The paper marks these as met by design; concrete exporters are
what turn that into a demonstrable claim.

Origem: novo (o notebook fazia isso ad hoc com GeoPandas).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import geopandas

__all__ = [
    "cells_to_geojson",
    "cell_to_wkt",
    "cells_to_wkt",
    "to_geodataframe",
    "from_geodataframe",
]


def cells_to_geojson(
    index: str,
    properties: dict[str, Any] | None = None,
    include_metadata: bool = True,
) -> dict[str, Any]:
    """Export cells as a GeoJSON FeatureCollection.

    Each cell becomes a Polygon feature in EPSG:4326. With
    ``include_metadata`` each feature carries its index, resolution,
    shape class, nominal and effective area, and its extension zone when
    it has one.

    Args:
        index: Compositional index string.
        properties: Extra properties copied onto every feature.
        include_metadata: Attach the per-cell ITACaRT properties.

    Returns:
        A GeoJSON FeatureCollection as a plain mapping.
    """
    raise NotImplementedError


def cell_to_wkt(cell: str) -> str | list[str]:
    """Cell boundary as a WKT polygon.

    Args:
        cell: Compositional index string.

    Returns:
        A WKT string for a single cell, or a positionally aligned list.
    """
    raise NotImplementedError


def cells_to_wkt(index: str, dissolve: bool = False) -> str:
    """Whole index as one WKT geometry.

    Args:
        index: Compositional index string.
        dissolve: Merge adjacent cells into a single polygon instead of
            emitting a MULTIPOLYGON of individual cells.

    Returns:
        A WKT string.
    """
    raise NotImplementedError


def to_geodataframe(index: str, crs: str = "EPSG:4326") -> "geopandas.GeoDataFrame":
    """Export cells as a GeoDataFrame.

    Requires the ``geo`` extra.

    Args:
        index: Compositional index string.
        crs: Output CRS. ``EPSG:4326`` for mapping, or the sinusoidal
            PROJ string to inspect cells undistorted on the plane.

    Returns:
        One row per cell, with geometry and ITACaRT metadata columns.

    Raises:
        ImportError: If GeoPandas is not installed.
    """
    raise NotImplementedError


def from_geodataframe(
    gdf: "geopandas.GeoDataFrame",
    resolution: int,
    containment: str = "center",
) -> list[str]:
    """Fill every geometry of a GeoDataFrame into a compositional index.

    Requires the ``geo`` extra.

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
    raise NotImplementedError
