"""The OGC-facing engine object.

Carries OGC DGGS Core requirements 6 (harmonized model) and 7 (defined
CRS), which are properties of the system rather than of any one cell, and
so have nowhere else to live.

The functional API is complete on its own; this class exists to expose
ITACaRT as a described, introspectable DGGS to standards-aware consumers
such as an OGC API-DGGS server.

Origem: itacart_core/engine.py (IDGGSEngine), ampliada para cobrir
os requisitos declarativos.
"""

from __future__ import annotations

from typing import Any

__all__ = ["ITACaRT", "describe", "crs", "conformance"]


class ITACaRT:
    """Stateful facade over the ITACaRT reference system.

    Wraps the functional API so a caller can fix defaults once and reuse
    them, and provides the descriptive metadata the OGC model asks for.

    Args:
        default_resolution: Resolution assumed when a call omits one.
        edge_model: Default edge interpretation for geometry operations.
        max_segment_m: Default densification bound in metres.
    """

    def __init__(
        self,
        default_resolution: int = 13,
        edge_model: str = "WGS84_GEODESIC",
        max_segment_m: float = 1000.0,
    ) -> None:
        raise NotImplementedError

    # -- Descriptive (OGC Req 6, 7) -------------------------------------

    def describe(self) -> dict[str, Any]:
        """Machine-readable description of the reference system.

        Covers identity and DOI, the CRS and datum, the tessellation
        method, cell geometry, the full resolution table, refinement
        ratios, and the boundary treatments. Shaped to feed an OGC
        API-DGGS ``/dggs/{dggrsId}`` response.

        Returns:
            The description as a plain mapping.
        """
        raise NotImplementedError

    def crs(self) -> dict[str, Any]:
        """The coordinate reference system this grid is defined on.

        Reports WGS84 as the datum, satisfying requirement 7 and the
        GNSS-compatibility design criterion, alongside the ellipsoidal
        sinusoidal projection used internally.

        Returns:
            A mapping with datum, ellipsoid parameters, the PROJ string
            and the projection's defining equations.
        """
        raise NotImplementedError

    def conformance(self) -> dict[str, Any]:
        """Self-reported conformance against DGGS Core and EAERS.

        Mirrors Frames 4 and 5 of the paper: Core fully met, EAERS
        partially, with the divergences recorded as deliberate design
        trade-offs. Requirements 22 to 25 are not met because ITACaRT
        tessellates the ellipsoid directly instead of projecting a
        polyhedron.

        Returns:
            One record per requirement, each with an identifier, a
            status and a justification.
        """
        raise NotImplementedError

    # -- Delegating operations ------------------------------------------

    def geo_to_cell(self, lon: float, lat: float, resolution: int | None = None) -> str:
        """Address the cell containing a position.

        See :func:`itacart.cells.geo_to_cell`.
        """
        raise NotImplementedError

    def cell_to_centroid(self, cell: str) -> tuple[float, float]:
        """Geodetic centroid of a cell.

        See :func:`itacart.cells.cell_to_centroid`.
        """
        raise NotImplementedError

    def cell_to_boundary(self, cell: str) -> list[tuple[float, float]]:
        """Geodetic vertices bounding a cell.

        See :func:`itacart.cells.cell_to_boundary`.
        """
        raise NotImplementedError


def describe() -> dict[str, Any]:
    """Module-level shortcut to :meth:`ITACaRT.describe` with defaults."""
    raise NotImplementedError


def crs() -> dict[str, Any]:
    """Module-level shortcut to :meth:`ITACaRT.crs` with defaults."""
    raise NotImplementedError


def conformance() -> dict[str, Any]:
    """Module-level shortcut to :meth:`ITACaRT.conformance` with defaults."""
    raise NotImplementedError
