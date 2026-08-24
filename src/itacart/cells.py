"""Quantization and inverse geometry: the core address <-> position mapping.

Satisfies OGC DGGS Core requirements 11 (simple cell geometry), 12 (direct
position) and 16 (quantization).

Origem: itacart_core/cells.py + sinusoidal_coordinates_to_dggs do notebook.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from shapely.geometry import Polygon

__all__ = [
    "geo_to_cell",
    "cell_to_anchor",
    "cell_to_centroid",
    "cell_to_boundary",
    "cell_to_polygon",
    "sinusoidal_to_cell",
    "cell_to_sinusoidal",
]


# --------------------------------------------------------------------------
# Quantization (OGC Req 16)
# --------------------------------------------------------------------------


def geo_to_cell(lon: float, lat: float, resolution: int) -> str:
    """Address the cell containing a geodetic position.

    Projects onto the ellipsoidal sinusoidal plane, then descends the
    hierarchy: the resolution-1 Cartesian pair, then alternating 1-to-4
    and 1-to-25 refinements. The parallelogram-to-square residual shear
    is applied at each step so the descent stays in integer arithmetic.

    Boundary cases are handled by :mod:`itacart.boundary`: positions on
    the prime meridian fall in triangular cells, and positions beyond the
    antemeridian resolve only inside a defined extension zone.

    Args:
        lon: Longitude in decimal degrees.
        lat: Latitude in decimal degrees.
        resolution: Target resolution level, 1 to 13.

    Returns:
        The atomic compositional index of the containing cell.

    Raises:
        ResolutionError: If ``resolution`` is out of range.
        DomainError: If the position is outside the addressable domain.
        AntemeridianError: If the position lies beyond the antemeridian
            and outside every extension zone.
    """
    raise NotImplementedError


def sinusoidal_to_cell(x: float, y: float, resolution: int) -> str:
    """Address a cell from projection-plane coordinates.

    Skips the forward projection; useful when coordinates are already on
    the plane, as in bulk cell filling.

    Args:
        x: Easting in metres on the sinusoidal plane.
        y: Northing in metres on the sinusoidal plane.
        resolution: Target resolution level, 1 to 13.

    Returns:
        The atomic compositional index of the containing cell.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------
# Inverse geometry (OGC Req 11, 12)
# --------------------------------------------------------------------------


def cell_to_anchor(cell: str) -> tuple[float, float] | list[tuple[float, float]]:
    """Representative position of a cell, as defined by ITACaRT.

    ITACaRT designates the lower-left vertex, not the centroid, so that
    addressing behaves like a Cartesian system for surveyors. This is the
    position that satisfies OGC Core requirement 12, and it is also why
    EAERS requirement 27 is only partially met.

    Use :func:`cell_to_centroid` when a centre point is wanted.

    Args:
        cell: Compositional index string.

    Returns:
        ``(lon, lat)`` for a single cell, or a positionally aligned list.
    """
    raise NotImplementedError


def cell_to_centroid(cell: str) -> tuple[float, float] | list[tuple[float, float]]:
    """Geodetic centroid of a cell.

    Computed from the actual vertices, so triangular and trapezoidal
    cells report their own centroid rather than a parallelogram
    approximation.

    Args:
        cell: Compositional index string.

    Returns:
        ``(lon, lat)`` for a single cell, or a positionally aligned list.
    """
    raise NotImplementedError


def cell_to_boundary(
    cell: str, close: bool = False
) -> list[tuple[float, float]] | list[list[tuple[float, float]]]:
    """Geodetic vertices bounding a cell, counter-clockwise.

    Vertex count depends on :func:`itacart.boundary.cell_shape`: four for
    a parallelogram, three for a prime-meridian triangle, four for a
    trapezoid with one base clipped at the boundary line.

    Args:
        cell: Compositional index string.
        close: Repeat the first vertex at the end, as ring conventions
            such as GeoJSON require.

    Returns:
        A vertex list for a single cell, or a positionally aligned list
        of vertex lists.
    """
    raise NotImplementedError


def cell_to_polygon(cell: str) -> "Polygon | list[Polygon]":
    """Cell boundary as a Shapely polygon in EPSG:4326.

    Args:
        cell: Compositional index string.

    Returns:
        A polygon for a single cell, or a positionally aligned list.
    """
    raise NotImplementedError


def cell_to_sinusoidal(
    cell: str,
) -> tuple[float, float] | list[tuple[float, float]]:
    """Anchor of a cell on the projection plane.

    Args:
        cell: Compositional index string.

    Returns:
        ``(x, y)`` in metres for a single cell, or a positionally aligned
        list.
    """
    raise NotImplementedError
