"""Neighbourhood and adjacency on the parallelogram lattice.

Completes OGC DGGS Core requirement 17. Adjacency is resolved purely from
the index strings, with no coordinate arithmetic, following section 3.1 of
the paper:

- **Resolutions 0 and 1** - integer arithmetic on the ``XXXX/YYYY`` pair.
  Quadrant borders apply a deflection rule when ``X`` or ``Y`` reaches
  zero or a meridian boundary.
- **Even resolutions** - a 2x2 quad-tree pattern. Within one parent,
  vertical adjacency shifts the code by 2 and horizontal by 1. Across
  parents, ascend, step to the parent's neighbour, descend to the
  matching child.
- **Odd resolutions** - a 5x5 grid. Horizontal neighbours step the
  numeric component with wrap-around between columns 5 and 1, vertical
  neighbours step the alphabetic component with wrap-around between rows
  E and A.

Grid distance is measured on the lattice, not on the ellipsoid. On a
sheared parallelogram grid the two diverge, and they diverge more at high
latitudes and longitudes. Callers coming from H3 tend to assume otherwise.

Origem: novo (o notebook e o itacart_core nao implementam vizinhanca).
"""

from __future__ import annotations

from typing import Literal

__all__ = [
    "grid_disk",
    "grid_ring",
    "grid_distance",
    "are_neighbor_cells",
    "get_neighbor",
    "cells_to_directed_edge",
    "directed_edge_to_cells",
    "cell_to_edges",
    "deflect",
]

Metric = Literal["chebyshev", "manhattan"]
Direction = Literal["N", "S", "E", "W", "NE", "NW", "SE", "SW"]


def grid_disk(
    index: str,
    k_distance: int = 1,
    metric: Metric = "chebyshev",
    dedupe: bool = False,
    flatten: bool = False,
) -> list[str] | list[list[str]]:
    """Cells within ``k`` grid steps of the origin, the origin included.

    ``chebyshev`` admits diagonal steps, giving a filled rhombus on the
    lattice; ``manhattan`` admits only axis steps, giving a filled
    diamond.

    Args:
        index: Compositional index string.
        k_distance: Ring radius in grid steps.
        metric: Lattice metric to expand under.
        dedupe: When the input holds several cells, drop cells that
            appear in more than one disk. Breaks positional alignment,
            so it is off by default; the union is also obtainable by
            composing the aligned result and normalizing it.
        flatten: When the input holds several cells, ``False`` returns
            one list per input cell; ``True`` returns a single flat list.

    Returns:
        A cell list for a single origin, or a list of lists aligned with
        :func:`itacart.index.decompose`.

    Raises:
        ValueError: If ``k_distance`` is negative.
    """
    raise NotImplementedError


def grid_ring(
    index: str,
    k_distance: int = 1,
    metric: Metric = "chebyshev",
    dedupe: bool = False,
    flatten: bool = False,
) -> list[str] | list[list[str]]:
    """Cells at exactly ``k`` grid steps, the hollow shell of a disk.

    Args:
        index: Compositional index string.
        k_distance: Ring radius in grid steps.
        metric: Lattice metric to expand under.
        dedupe: Drop cells shared between rings of different origins.
        flatten: Return one flat list instead of one list per origin.

    Returns:
        A cell list for a single origin, or a positionally aligned list
        of lists.
    """
    raise NotImplementedError


def grid_distance(origin: str, destination: str, metric: Metric = "chebyshev") -> int:
    """Steps between two cells on the lattice.

    Not a geodesic distance. Use :func:`itacart.geodesy.inverse_geodesic`
    on the anchors or centroids for metric separation.

    Args:
        origin: Atomic index of the first cell.
        destination: Atomic index of the second cell.
        metric: Lattice metric to measure under.

    Returns:
        Step count, zero when the cells coincide.

    Raises:
        ResolutionError: If the cells sit at different resolutions.
        DomainError: If a path cannot be resolved across the boundary
            between the cells.
    """
    raise NotImplementedError


def get_neighbor(index: str, direction: Direction) -> str | None | list[str | None]:
    """Single adjacent cell in a given direction.

    Directions are lattice directions in the origin cell's own quadrant,
    already accounting for the mirroring that ITACaRT applies across the
    axes, so ``"N"`` means one step toward the pole in every quadrant.

    Args:
        index: Compositional index string.
        direction: One of the eight lattice directions.

    Returns:
        The neighbour index, or ``None`` where no cell exists on that
        side, for a single cell; or a positionally aligned list.
    """
    raise NotImplementedError


def are_neighbor_cells(origin: str, destination: str) -> bool:
    """Whether two cells share an edge.

    Edge adjacency only; cells meeting at a single vertex are not
    neighbours.

    Args:
        origin: Atomic index of the first cell.
        destination: Atomic index of the second cell.

    Returns:
        ``True`` when the cells are edge-adjacent.
    """
    raise NotImplementedError


def deflect(cell: str, direction: Direction) -> str | None:
    """Apply the boundary deflection rule for a step that leaves a quadrant.

    Isolated from :func:`get_neighbor` because it is the subtle part: a
    step across the prime meridian lands on a triangular cell, a step
    across a quadrant axis mirrors the coordinate, and a step across the
    antemeridian only resolves inside an extension zone.

    Args:
        cell: Atomic index string at the quadrant border.
        direction: Lattice direction of the step.

    Returns:
        The deflected neighbour index, or ``None`` when the step leaves
        the addressable domain.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------
# Directed edges
# --------------------------------------------------------------------------


def cells_to_directed_edge(origin: str, destination: str) -> str | list[str]:
    """Identifier of the directed edge between two adjacent cells.

    Validates edge adjacency first. Origins and destinations are paired
    by position, so N origins and N destinations give N edges; broadcast
    a single origin by repeating it.

    Args:
        origin: Compositional index of the origin cell or cells.
        destination: Compositional index of the destination cell or cells.

    Returns:
        An edge identifier for a single pair, or a positionally aligned
        list.

    Raises:
        ValueError: If the two indices hold different cell counts.
        DomainError: If any pair is not edge-adjacent.
    """
    raise NotImplementedError


def directed_edge_to_cells(edge: str) -> tuple[str, str] | list[tuple[str, str]]:
    """Recover the origin and destination of a directed edge.

    Args:
        edge: Directed edge identifier.

    Returns:
        ``(origin, destination)`` for a single edge, or a positionally
        aligned list.
    """
    raise NotImplementedError


def cell_to_edges(cell: str) -> list[str] | list[list[str]]:
    """Every directed edge leaving a cell.

    Four for a parallelogram or trapezoid, three for a prime-meridian
    triangle.

    Args:
        cell: Compositional index string.

    Returns:
        An edge list for a single cell, or a positionally aligned list of
        lists.
    """
    raise NotImplementedError
