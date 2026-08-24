"""Boundary behaviour: prime meridian, antemeridian, extension zones.

This module is what makes the domain global, complete and unique, and so
carries OGC DGGS Core requirements 8, 9 and 10 on its own.

Three departures from the uniform parallelogram grid are specified in
section 3.2 of the paper:

**Prime meridian.** Cells straddling 0 degrees are isosceles triangles
whose base is twice their height, mirrored about the meridian. This
applies at resolution 1; finer resolutions subdivide the adjacent eastern
cells rather than creating separate western ones. Consequently
resolution-1 cells with ``X = 0`` in the western quadrants do not exist.

**Antemeridian.** Uniform areas cannot be held everywhere at 180 degrees,
so the eastern quadrants are extended over the inhabited landmasses that
the meridian crosses: Fiji and the Chukotka/Wrangel group. Antarctica is
excluded. Cells in the oceans and in Antarctica along that boundary
therefore have unequal areas.

**Trapezoids.** Where a parallelogram vertex on the side opposite the
prime meridian would cross a boundary line, that side is clipped to the
line, extending the longer base into a trapezoid. The rule applies only
to the cells in that condition, not to their subdivisions.

Origem: novo (nao implementado no notebook nem no itacart_core).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .constants import CellShape, ExtensionZone

if TYPE_CHECKING:
    from shapely.geometry.base import BaseGeometry

__all__ = [
    "cell_shape",
    "is_valid_cell",
    "is_boundary_cell",
    "is_triangular_cell",
    "is_trapezoidal_cell",
    "is_extension_cell",
    "extension_zone",
    "extension_zone_for_point",
    "extension_bounds",
    "crosses_antemeridian",
    "is_equal_area_cell",
]


# --------------------------------------------------------------------------
# Shape classification
# --------------------------------------------------------------------------


def cell_shape(cell: str) -> CellShape | list[CellShape]:
    """Geometric class of a cell.

    Drives vertex count in :func:`itacart.cells.cell_to_boundary` and area
    computation in :func:`itacart.resolutions.effective_cell_area`.

    Args:
        cell: Compositional index string.

    Returns:
        One of ``"parallelogram"``, ``"triangle"`` or ``"trapezoid"`` for
        a single cell, or a positionally aligned list.
    """
    raise NotImplementedError


def is_triangular_cell(cell: str) -> bool | list[bool]:
    """Whether a cell is a prime-meridian isosceles triangle.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    raise NotImplementedError


def is_trapezoidal_cell(cell: str) -> bool | list[bool]:
    """Whether a cell has been clipped at a boundary line.

    Trapezoidal cells break the equal-area guarantee; this is the check
    that should gate any area-sensitive cadastral computation.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    raise NotImplementedError


def is_equal_area_cell(cell: str) -> bool | list[bool]:
    """Whether a cell honours the nominal area of its resolution.

    True for parallelograms and prime-meridian triangles, false for
    trapezoids. Complement of :func:`is_trapezoidal_cell`, named for the
    property callers actually care about.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------
# Domain validity (OGC Req 8-10)
# --------------------------------------------------------------------------


def is_valid_cell(cell: str) -> bool | list[bool]:
    """Whether a cell exists in the ITACaRT domain.

    Stricter than :func:`itacart.index.is_valid_index`, which only checks
    syntax. This additionally applies:

    - the western-quadrant ``X = 0`` non-existence rule,
    - the resolution-1 extent bounds per quadrant,
    - antemeridian coverage, so oceanic cells beyond 180 degrees outside
      an extension zone are rejected.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    raise NotImplementedError


def is_boundary_cell(cell: str) -> bool | list[bool]:
    """Whether a cell touches any grid discontinuity.

    True for prime-meridian triangles, antemeridian cells and
    extension-zone edge cells. A coarse screen; use :func:`cell_shape` or
    :func:`extension_zone` when the specific condition matters.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------
# Extension zones
# --------------------------------------------------------------------------


def is_extension_cell(cell: str) -> bool | list[bool]:
    """Whether a cell lies inside an eastern-quadrant extension.

    Args:
        cell: Compositional index string.

    Returns:
        A boolean for a single cell, or a positionally aligned list.
    """
    raise NotImplementedError


def extension_zone(cell: str) -> ExtensionZone | None | list[ExtensionZone | None]:
    """Which extension zone contains a cell, if any.

    Callers that need to branch on the specific zone should use this
    rather than :func:`is_extension_cell`, since the two zones sit in
    different quadrants and have different longitude limits.

    Args:
        cell: Compositional index string.

    Returns:
        ``"FIJI"``, ``"CHUKOTKA"`` or ``None`` for a single cell, or a
        positionally aligned list.
    """
    raise NotImplementedError


def extension_zone_for_point(lon: float, lat: float) -> ExtensionZone | None:
    """Which extension zone contains a geodetic position, if any.

    The point-level counterpart of :func:`extension_zone`, called during
    quantization before any cell exists.

    Args:
        lon: Longitude in decimal degrees.
        lat: Latitude in decimal degrees.

    Returns:
        The zone name, or ``None`` outside both zones.
    """
    raise NotImplementedError


def extension_bounds(zone: ExtensionZone) -> tuple[float, float, float, float]:
    """Geodetic bounds of an extension zone.

    Args:
        zone: ``"FIJI"`` or ``"CHUKOTKA"``.

    Returns:
        ``(min_lon, min_lat, max_lon, max_lat)`` in decimal degrees, with
        ``min_lon`` on the western side of the antemeridian.

    Raises:
        ValueError: If ``zone`` is not a defined zone.
    """
    raise NotImplementedError


def crosses_antemeridian(geometry: "BaseGeometry") -> bool:
    """Whether a geometry spans the 180th meridian.

    Called by :func:`itacart.geometry.polyfill` before descent, so the
    caller gets :class:`~itacart.exceptions.AntemeridianError` rather than
    a silently truncated cell set.

    Args:
        geometry: A Shapely geometry in EPSG:4326.

    Returns:
        ``True`` if the geometry crosses 180 degrees longitude.
    """
    raise NotImplementedError
