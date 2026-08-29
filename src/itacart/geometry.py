"""Vector geometry against the grid: filling, vertex mapping, canonical form.

Two ways to represent a vector feature, both described in section 4 of the
paper:

**Cell filling** - the feature is the set of cells it covers. Area follows
from a cell count with no projection distortion, but the index grows
verbose at fine resolutions.

**Vertex representation** - only the cells holding the polygon vertices are
kept, as in conventional vector data. Far more compact, and the form the
binary encodings in :mod:`itacart.serialization` are built on.

Provenance: ``itacart_core/cell_filling.py``, ``densification.py``,
``geometry_blob.py`` (``canonicalize_rings``) and
``cadastral_processor/vertex_extractor.py``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from shapely.geometry import LineString, Point, Polygon
    from shapely.geometry.base import BaseGeometry

__all__ = [
    "polyfill",
    "count_internal_cells",
    "vertex_to_cell",
    "cells_to_geometry",
    "densify_orthodromic",
    "densify_segment",
    "canonicalize_rings",
]

Containment = Literal["center", "intersects", "contains"]


# --------------------------------------------------------------------------
# Cell filling
# --------------------------------------------------------------------------


def polyfill(
    geometry: "BaseGeometry",
    resolution: int,
    containment: Containment = "center",
    compact: bool = False,
    n_jobs: int = 1,
) -> str:
    """Rasterise a geometry into a compositional index.

    Descends the hierarchy over the sinusoidal plane rather than testing
    cells one by one, so the cost tracks the boundary rather than the
    area.

    ``center`` keeps cells whose anchor falls inside, ``intersects``
    keeps every cell touching the geometry, and ``contains`` keeps only
    cells wholly inside. For a legally meaningful parcel area,
    ``contains`` paired with :func:`count_internal_cells` gives the strict
    interior.

    Args:
        geometry: A Shapely geometry in EPSG:4326. Densify it first with
            :func:`densify_orthodromic` when its edges are geodesics
            rather than straight lines on the plane.
        resolution: Target resolution level, 1 to 13.
        containment: Predicate deciding whether a cell is kept.
        compact: Return a mixed-resolution compacted index instead of a
            uniform one.
        n_jobs: Worker count; above 1 requires the ``parallel`` extra.

    Returns:
        A compositional index string covering the geometry.

    Raises:
        AntemeridianError: If the geometry crosses 180 degrees outside an
            extension zone.
        UnsupportedGeometryTypeError: On unsupported geometry types.
    """
    raise NotImplementedError


def count_internal_cells(polygon: "Polygon", resolution: int, n_jobs: int = 1) -> int:
    """Count cells strictly inside a polygon.

    Fast path over :func:`polyfill` with ``containment="contains"``: the
    count is accumulated during descent without materialising indices,
    which matters at resolution 13 where a modest parcel holds millions
    of cells.

    Since every cell carries the same area, the count times
    :func:`itacart.resolutions.nominal_cell_area` gives a distortion-free
    area, which is the property the paper builds tokenization on.

    Provenance: ``itacart_core/engine.py`` (``polygon_to_cells_count``).

    Args:
        polygon: A Shapely polygon in EPSG:4326.
        resolution: Target resolution level, 1 to 13.
        n_jobs: Worker count; above 1 requires the ``parallel`` extra.

    Returns:
        Number of cells wholly inside the polygon.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------
# Vertex representation
# --------------------------------------------------------------------------


def vertex_to_cell(
    geometry: "Point | LineString | Polygon",
    resolution: int,
    dedupe_consecutive: bool = True,
) -> list[str]:
    """Map geometry vertices to cells, preserving sequence.

    Order is topological, not sorted: it is what allows the original
    geometry to be reconstructed. Rings keep their winding, and holes
    follow the exterior.

    Consecutive vertices landing in the same cell are collapsed by
    default, since at resolution 13 that reflects survey precision rather
    than distinct corners. Non-consecutive repeats are kept, being a
    genuine self-touching pathology the caller should see.

    Provenance: ``cadastral_processor/vertex_extractor.py``.

    Args:
        geometry: A Shapely geometry in EPSG:4326.
        resolution: Target resolution level, 1 to 13.
        dedupe_consecutive: Collapse consecutive duplicates.

    Returns:
        Atomic index strings in traversal order.
    """
    raise NotImplementedError


def cells_to_geometry(
    cells: list[str], geometry_type: str = "Polygon"
) -> "BaseGeometry":
    """Rebuild a geometry from an ordered vertex cell list.

    Inverse of :func:`vertex_to_cell`. Reconstruction lands on cell
    anchors, so it is exact only to the resolution used: at resolution 13
    that is 1 cm.

    Args:
        cells: Atomic index strings in traversal order.
        geometry_type: OGC SFA type to build.

    Returns:
        A Shapely geometry in EPSG:4326.

    Raises:
        UnsupportedGeometryTypeError: On unsupported types.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------
# Densification
# --------------------------------------------------------------------------


def densify_orthodromic(polygon: "Polygon", max_segment_m: float = 1000.0) -> "Polygon":
    """Insert intermediate vertices along geodesics.

    Applied to the exterior ring and every hole, with
    ``n_segments = floor(d_geo / max_segment_m) + 1`` and points equally
    spaced in geodesic distance.

    Needed because a straight line on the sinusoidal plane is not a
    geodesic on the ellipsoid; without densification a long edge would
    fill the wrong cells in between.

    Provenance: ``itacart_core/densification.py``.

    Args:
        polygon: A Shapely polygon in EPSG:4326.
        max_segment_m: Longest segment to leave undensified, in metres.

    Returns:
        A densified polygon in EPSG:4326.

    Raises:
        DensificationError: If a segment cannot be densified.
    """
    raise NotImplementedError


def densify_segment(
    p1: tuple[float, float],
    p2: tuple[float, float],
    max_segment_m: float,
    edge_model: str = "WGS84_GEODESIC",
) -> list[tuple[float, float]]:
    """Densify a single segment.

    The building block of :func:`densify_orthodromic`, exposed for
    callers working segment by segment such as open LINESTRING handling.

    Provenance: ``itacart_core/geometry_blob.py`` (``densify_segment``).

    Args:
        p1: ``(lon, lat)`` of the start point.
        p2: ``(lon, lat)`` of the end point.
        max_segment_m: Longest segment to leave undensified, in metres.
        edge_model: Edge interpretation, currently ``"WGS84_GEODESIC"``.

    Returns:
        ``(lon, lat)`` pairs including both endpoints.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------
# Canonicalization
# --------------------------------------------------------------------------


def canonicalize_rings(rings: list[list[str]]) -> list[list[str]]:
    """Normalise rings to a single spelling per geometry.

    Rings are rotated to their minimum lexicographic cyclic rotation, so
    the same ring starting at a different vertex canonicalises to the
    same sequence. Rotation applies to closed rings only; LINESTRING is
    directional and is left as given.

    This is what makes a geometry content-addressable, and therefore what
    makes hashing it meaningful.

    Provenance: ``itacart_core/geometry_blob.py``
    (``MIN_LEX_CYCLIC_ROTATION``, v2).

    Args:
        rings: Rings as lists of atomic index strings, exterior first.

    Returns:
        The canonicalised rings, in input order.
    """
    raise NotImplementedError
