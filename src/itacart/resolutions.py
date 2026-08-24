"""Resolution table helpers.

Satisfies OGC DGGS Core requirements 14 and 15 (hierarchical grid sequence).

Origem: itacart_core/resolutions.py + Tabela 1 do artigo.
"""

from __future__ import annotations

from typing import Literal

__all__ = [
    "get_resolution",
    "refinement_ratio",
    "cell_size",
    "nominal_cell_area",
    "effective_cell_area",
    "scale_for_resolution",
    "resolution_for_scale",
    "resolution_table",
    "is_tokenizable_resolution",
]

ScaleKind = Literal["visualization", "analysis"]


def get_resolution(cell: str) -> int:
    """Return the resolution level of a cell index.

    For a compositional index holding several cells, every terminal cell
    must sit at the same level; otherwise the index is mixed-resolution
    and the caller should use :func:`itacart.index.decompose` first.

    Args:
        cell: Compositional index string.

    Returns:
        Resolution level, 0 to 13.

    Raises:
        InvalidIndexError: If the index is malformed.
        ResolutionError: If terminal cells sit at differing resolutions.
    """
    raise NotImplementedError


def refinement_ratio(resolution: int) -> int:
    """Number of children produced when refining *into* ``resolution``.

    Even resolutions come from a 1-to-4 subdivision, odd resolutions from
    a 1-to-25 subdivision. Resolution 1 is a special case: it is the base
    10 km grid addressed by Cartesian integers, not a refinement.

    Args:
        resolution: Target resolution level, 2 to 13.

    Returns:
        ``4`` for even resolutions, ``25`` for odd ones.

    Raises:
        ResolutionError: If ``resolution`` is below 2 or above 13.
    """
    raise NotImplementedError


def cell_size(resolution: int) -> float:
    """Base and height of a cell, in metres.

    Args:
        resolution: Resolution level, 1 to 13.

    Returns:
        Edge length in metres, from 10 000 m down to 0.01 m.

    Raises:
        ResolutionError: If ``resolution`` is 0 or out of range.
    """
    raise NotImplementedError


def nominal_cell_area(resolution: int) -> float:
    """Area of a standard parallelogram cell, in square metres.

    This is the value from Table 1 of the paper. It is exact for every
    parallelogram cell anywhere on the ellipsoid, and for the triangular
    prime-meridian cells, whose base is twice their height.

    It is NOT correct for trapezoidal cells at the antemeridian and at
    extension-zone boundaries; use :func:`effective_cell_area` when the
    cell may be one of those.

    Args:
        resolution: Resolution level, 1 to 13.

    Returns:
        Area in square metres.

    Raises:
        ResolutionError: If ``resolution`` is 0 or out of range.
    """
    raise NotImplementedError


def effective_cell_area(cell: str) -> float | list[float]:
    """True area of a specific cell, accounting for boundary clipping.

    Equals :func:`nominal_cell_area` for parallelogram and triangular
    cells. For trapezoidal cells the longer base is clipped at the
    boundary line, so the area is computed from the actual vertices.

    Cadastral use makes this distinction legally significant: returning
    the nominal area for a clipped cell would misstate a parcel.

    Args:
        cell: Compositional index string.

    Returns:
        Area in square metres for a single terminal cell, or a
        positionally aligned list when the index holds several.
    """
    raise NotImplementedError


def scale_for_resolution(resolution: int, kind: ScaleKind = "visualization") -> int:
    """Cartographic scale denominator matching a resolution.

    ``visualization`` follows the minimum visible line of 0.1 mm on paper
    maps (Jenny et al., 2008). ``analysis`` follows sampling theory,
    ``scale = resolution * 2 * 1000`` (Tobler, 1987).

    Args:
        resolution: Resolution level, 1 to 13.
        kind: Which of the two scale families to report.

    Returns:
        Scale denominator, e.g. ``10_000`` for 1:10 000.

    Raises:
        ResolutionError: If ``resolution`` is 0 or out of range.
    """
    raise NotImplementedError


def resolution_for_scale(denominator: int, kind: ScaleKind = "visualization") -> int:
    """Coarsest resolution adequate for a target cartographic scale.

    Inverse of :func:`scale_for_resolution`; the practical entry point for
    a surveyor who knows the map scale but not the grid level.

    Args:
        denominator: Scale denominator, e.g. ``1000`` for 1:1 000.
        kind: Which of the two scale families to interpret against.

    Returns:
        Resolution level, 1 to 13.
    """
    raise NotImplementedError


def resolution_table() -> list[dict[str, object]]:
    """Full Table 1 of the paper as structured records.

    Each record carries ``resolution``, ``cell_size_m``, ``cell_area_m2``,
    ``refinement``, ``index_alphabet``, ``visualization_scale`` and
    ``analysis_scale``. Supports the OGC ``describe`` operation and gives
    documentation a single source of truth.

    Returns:
        One record per resolution, ordered 0 to 13.
    """
    raise NotImplementedError


def is_tokenizable_resolution(resolution: int) -> bool:
    """Whether cells at this resolution carry a whole-number metric area.

    Decimal convergence is what makes one token equal one standard unit
    of area in a blockchain registry. Resolution 13 gives exactly 1 cm2.

    Origem: itacart_core/engine.py (IDGGSEngine).

    Args:
        resolution: Resolution level, 1 to 13.

    Returns:
        ``True`` when the nominal area is a whole number in its natural
        metric unit.
    """
    raise NotImplementedError
