"""ITACaRT - ITA Cadastral Ellipsoidal Reference Tessellation.

An equal-area parallelogram Discrete Global Grid System for terrestrial
cadastral mapping, tessellated directly on the WGS84 ellipsoid.

Reference:
    Silva, I. N., Dietzsch, G., & Shiguemori, E. H. (2025). ITACaRT: An
    Equal-Area Parallelogram Discrete Global Grid System for Terrestrial
    Cadastral Mapping - Designed for Usability and Blockchain
    Integration. *Revista Brasileira de Cartografia*, 77.
    https://doi.org/10.14393/rbcv77n0a-79281

Quick start::

    import itacart

    cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=13)
    lon, lat = itacart.cell_to_centroid(cell)
    ring = itacart.cell_to_boundary(cell, close=True)

Every operation accepts a compositional index, which may address one cell
or a whole region. Operations returning one value per cell return a list
aligned with :func:`itacart.index.decompose` order.
"""

from __future__ import annotations

from typing import cast

__version__ = "0.1.0a4"
__paper_doi__ = "10.14393/rbcv77n0a-79281"

# -- Constants --------------------------------------------------------------

from .boundary import (
    absorbs_border,
    cell_shape,
    crosses_antemeridian,
    extension_bounds,
    extension_zone,
    extension_zone_for_point,
    is_boundary_cell,
    is_equal_area_cell,
    is_extension_cell,
    is_trapezoidal_cell,
    is_triangular_cell,
    is_valid_cell,
    last_lattice_column,
)
from .cells import (
    cell_to_anchor,
    cell_to_boundary,
    cell_to_centroid,
    cell_to_polygon,
    cell_to_sinusoidal,
    geo_to_cell,
    is_quadrant_boundary_cell,
    sinusoidal_to_cell,
)
from .constants import (
    ANALYSIS_SCALE,
    CELL_AREA_M2,
    CELL_SIZE_M,
    EXTENSION_ZONES,
    MAX_RESOLUTION,
    MIN_RESOLUTION,
    QUADRANTS,
    SINUSOIDAL_PROJ,
    VISUALIZATION_SCALE,
    WGS84_A,
    WGS84_E2,
    WGS84_F,
    refinement_alphabet,
)
from .engine import ITACaRT, conformance, crs, describe
from .exceptions import (
    AntemeridianError,
    ConvergenceError,
    DensificationError,
    DomainError,
    GeometryError,
    IncompatibleProfileError,
    InvalidIndexError,
    InvalidQuadrantError,
    InvalidRefinementCodeError,
    ITACaRTError,
    MalformedBlobError,
    MaxResolutionError,
    MinResolutionError,
    NonAtomicIndexError,
    NonExistentCellError,
    ResolutionError,
    SerializationError,
    UnsupportedGeometryTypeError,
)
from .geodesy import (
    direct_geodesic,
    geodetic_to_sinusoidal,
    inverse_geodesic,
    inverse_meridian_arc,
    meridian_arc,
    meridian_arc_quadrature,
    meridian_radius,
    prime_vertical_radius,
    sinusoidal_to_geodetic,
)
from .geometry import (
    canonicalize_rings,
    cells_to_geometry,
    count_internal_cells,
    densify_orthodromic,
    densify_segment,
    polyfill,
    vertex_to_cell,
)
from .hierarchy import (
    child_position,
    common_ancestor,
    compact_cells,
    contains,
    get_ancestors,
    get_children,
    get_descendants,
    get_parent,
    is_ancestor,
    uncompact_cells,
)
from .index import (
    base_cell_of,
    compose,
    count_cells,
    decompose,
    is_atomic,
    is_valid_index,
    iter_cells,
    normalize,
    parse,
    quadrant_of,
    split_components,
)
from .interop import (
    cell_to_wkt,
    cells_to_geojson,
    cells_to_wkt,
    from_geodataframe,
    to_geodataframe,
)
from .resolutions import (
    cell_size,
    effective_cell_area,
    get_resolution,
    is_tokenizable_resolution,
    linear_refinement_ratio,
    nominal_cell_area,
    refinement_ratio,
    resolution_for_scale,
    resolution_table,
    scale_for_resolution,
)
from .serialization import (
    count_vertices,
    decode_geometry,
    decode_tree,
    deserialize_from_blob,
    encode_geometry,
    encode_tree,
    geometry_hash,
    geometry_to_tree,
    is_ancestor_binary,
    iter_leaves,
    prefix_at_resolution_binary,
    resolution_of_binary,
    serialize_to_blob,
    validate_geometry,
    validate_tree,
)
from .topology import (
    are_neighbor_cells,
    cell_to_edges,
    cells_to_directed_edge,
    deflect,
    directed_edge_to_cells,
    get_neighbor,
    grid_disk,
    grid_distance,
    grid_ring,
)

# -- Exceptions -------------------------------------------------------------


# -- Core: projection and addressing ---------------------------------------


# -- Geodesy ----------------------------------------------------------------


# -- Index ------------------------------------------------------------------


# -- Resolutions ------------------------------------------------------------


# -- Hierarchy --------------------------------------------------------------


# -- Topology ---------------------------------------------------------------


# -- Boundary behaviour -----------------------------------------------------


# -- Vector geometry --------------------------------------------------------


# -- Serialization ----------------------------------------------------------


# -- Interoperability -------------------------------------------------------


# -- Engine -----------------------------------------------------------------


# -- H3-compatible aliases --------------------------------------------------
# Argument order follows H3 v4: (lat, lng) rather than (lon, lat).


def latlng_to_cell(lat: float, lng: float, resolution: int) -> str:
    """H3-style alias of :func:`geo_to_cell` taking ``(lat, lng)``."""
    return geo_to_cell(lng, lat, resolution)


def cell_to_latlng(cell: str) -> tuple[float, float]:
    """H3-style alias of :func:`cell_to_centroid` returning ``(lat, lng)``."""
    lon, lat = cast("tuple[float, float]", cell_to_centroid(cell))
    return (lat, lon)


cell_to_parent = get_parent
cell_to_children = get_children
cell_to_boundary_latlng = cell_to_boundary


__all__ = [
    # metadata
    "__version__",
    "__paper_doi__",
    # constants
    "WGS84_A",
    "WGS84_F",
    "WGS84_E2",
    "SINUSOIDAL_PROJ",
    "QUADRANTS",
    "MIN_RESOLUTION",
    "MAX_RESOLUTION",
    "CELL_SIZE_M",
    "CELL_AREA_M2",
    "VISUALIZATION_SCALE",
    "ANALYSIS_SCALE",
    "EXTENSION_ZONES",
    # exceptions
    "ITACaRTError",
    "InvalidIndexError",
    "NonAtomicIndexError",
    "ResolutionError",
    "DomainError",
    "NonExistentCellError",
    "AntemeridianError",
    "GeometryError",
    "SerializationError",
    "MalformedBlobError",
    # projection / addressing
    "geo_to_cell",
    "sinusoidal_to_cell",
    "cell_to_anchor",
    "cell_to_centroid",
    "cell_to_boundary",
    "cell_to_polygon",
    "cell_to_sinusoidal",
    # geodesy
    "geodetic_to_sinusoidal",
    "sinusoidal_to_geodetic",
    "inverse_geodesic",
    "direct_geodesic",
    # index
    "is_valid_index",
    "is_atomic",
    "decompose",
    "compose",
    "normalize",
    "count_cells",
    "iter_cells",
    "split_components",
    "quadrant_of",
    "base_cell_of",
    # resolutions
    "get_resolution",
    "refinement_ratio",
    "cell_size",
    "nominal_cell_area",
    "effective_cell_area",
    "scale_for_resolution",
    "resolution_for_scale",
    "resolution_table",
    "is_tokenizable_resolution",
    # hierarchy
    "get_parent",
    "get_children",
    "get_ancestors",
    "get_descendants",
    "is_ancestor",
    "common_ancestor",
    "contains",
    "compact_cells",
    "uncompact_cells",
    # topology
    "grid_disk",
    "grid_ring",
    "grid_distance",
    "are_neighbor_cells",
    "get_neighbor",
    "cells_to_directed_edge",
    "directed_edge_to_cells",
    "cell_to_edges",
    # boundary
    "cell_shape",
    "is_valid_cell",
    "is_boundary_cell",
    "is_triangular_cell",
    "is_trapezoidal_cell",
    "is_equal_area_cell",
    "is_extension_cell",
    "extension_zone",
    "extension_bounds",
    "crosses_antemeridian",
    "absorbs_border",
    "last_lattice_column",
    # geometry
    "polyfill",
    "count_internal_cells",
    "vertex_to_cell",
    "cells_to_geometry",
    "densify_orthodromic",
    "densify_segment",
    "canonicalize_rings",
    # serialization
    "encode_tree",
    "decode_tree",
    "serialize_to_blob",
    "deserialize_from_blob",
    "encode_geometry",
    "decode_geometry",
    "geometry_hash",
    "geometry_to_tree",
    # interop
    "cells_to_geojson",
    "cell_to_wkt",
    "cells_to_wkt",
    "to_geodataframe",
    "from_geodataframe",
    # engine
    "ITACaRT",
    "describe",
    "crs",
    "conformance",
    # H3-style aliases
    "latlng_to_cell",
    "cell_to_latlng",
    "cell_to_parent",
    "cell_to_children",
    "extension_zone_for_point",
    "is_quadrant_boundary_cell",
    "refinement_alphabet",
    "ConvergenceError",
    "DensificationError",
    "IncompatibleProfileError",
    "InvalidQuadrantError",
    "InvalidRefinementCodeError",
    "MaxResolutionError",
    "MinResolutionError",
    "UnsupportedGeometryTypeError",
    "inverse_meridian_arc",
    "meridian_arc",
    "meridian_arc_quadrature",
    "meridian_radius",
    "prime_vertical_radius",
    "child_position",
    "parse",
    "linear_refinement_ratio",
    "count_vertices",
    "is_ancestor_binary",
    "iter_leaves",
    "prefix_at_resolution_binary",
    "resolution_of_binary",
    "validate_geometry",
    "validate_tree",
    "deflect",
]
