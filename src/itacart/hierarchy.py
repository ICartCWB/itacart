"""Hierarchical navigation: ascent, descent, compaction.

Part of OGC DGGS Core requirement 17 (topological query functions), with
the neighbour half living in :mod:`itacart.topology`.

Every function here accepts a compositional index. Where the input holds
several cells, the return is positionally aligned with
:func:`itacart.index.decompose`.

Origem: itacart_core (analogos ASCII de binary_index) + notebook DGGS_Tree.
"""

from __future__ import annotations

from typing import Iterator

__all__ = [
    "get_parent",
    "get_children",
    "get_ancestors",
    "get_descendants",
    "is_ancestor",
    "common_ancestor",
    "contains",
    "compact_cells",
    "uncompact_cells",
    "child_position",
]


# --------------------------------------------------------------------------
# Ascent
# --------------------------------------------------------------------------


def get_parent(index: str, target_res: int | None = None) -> str | list[str]:
    """Ascend the hierarchy by lexical truncation.

    No floating-point work and no negative indices: the parent is a
    prefix of the child string, which is the property that makes the
    compositional index cheap to navigate.

    Args:
        index: Compositional index string.
        target_res: Resolution to ascend to. Defaults to one level up.

    Returns:
        The ancestor index for a single cell, or a positionally aligned
        list. Duplicates are preserved so alignment holds; call
        :func:`itacart.index.normalize` on the composed result to dedupe.

    Raises:
        MinResolutionError: If ``target_res`` is above resolution 0.
        ResolutionError: If ``target_res`` is finer than the input.
    """
    raise NotImplementedError


def get_ancestors(index: str) -> list[str] | list[list[str]]:
    """Every ancestor of a cell, coarse to fine.

    Args:
        index: Compositional index string.

    Returns:
        An ancestor chain for a single cell, or a positionally aligned
        list of chains.
    """
    raise NotImplementedError


def common_ancestor(index: str) -> str:
    """Deepest cell containing every terminal cell of the index.

    Underlies :func:`compact_cells` and :func:`itacart.index.compose`, and
    is independently useful for finding the grouping node of a vertex set.

    Args:
        index: Compositional index string, or a composed list of cells.

    Returns:
        The deepest common ancestor index. Returns the quadrant when the
        cells share only that, and raises when they do not even share a
        quadrant.

    Raises:
        DomainError: If the cells span more than one quadrant.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------
# Descent
# --------------------------------------------------------------------------


def get_children(
    index: str, target_res: int | None = None, flatten: bool = False
) -> Iterator[str] | Iterator[list[str]]:
    """Project direct descendants of an index.

    Yields 4 children when descending into an even resolution and 25 when
    descending into an odd one.

    Args:
        index: Compositional index string.
        target_res: Resolution to descend to. Defaults to one level down.
        flatten: When the input holds several cells, ``False`` yields one
            list per input cell, preserving positional alignment;
            ``True`` yields a single flat stream.

    Yields:
        Child index strings, or lists of them when ``flatten`` is
        ``False`` and the input is compositional.

    Raises:
        MaxResolutionError: If ``target_res`` exceeds resolution 13.
    """
    raise NotImplementedError


def get_descendants(index: str, target_res: int) -> Iterator[str]:
    """Stream every descendant at a target resolution.

    Cardinality grows fast: from resolution 1 to 13 a single base cell
    holds 10^12 cells. Always a generator, never a list.

    Args:
        index: Compositional index string.
        target_res: Resolution to expand to.

    Yields:
        Atomic index strings at ``target_res``.
    """
    raise NotImplementedError


def child_position(cell: str) -> int | list[int]:
    """Ordinal of a cell among its siblings.

    Zero-based against the refinement alphabet of the cell's resolution:
    ``1``-``4`` for even levels, ``A1``-``E5`` in row-major order for odd.

    Args:
        cell: Compositional index string.

    Returns:
        An ordinal for a single cell, or a positionally aligned list.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------
# Containment
# --------------------------------------------------------------------------


def is_ancestor(parent_index: str, child_index: str) -> bool:
    """Whether one cell contains another in the hierarchy.

    Origem: itacart_core/engine.py.

    Args:
        parent_index: Candidate ancestor, a single atomic index.
        child_index: Candidate descendant, a single atomic index.

    Returns:
        ``True`` if ``parent_index`` is a strict prefix of
        ``child_index``.
    """
    raise NotImplementedError


def contains(region: str, cell: str) -> bool | list[bool]:
    """Whether a compositional region covers a cell.

    Complements :func:`is_ancestor`, which only tests the vertical
    relation between two single cells. A region covers a cell when any of
    its terminal cells is that cell or an ancestor of it.

    Args:
        region: Compositional index string denoting the region.
        cell: Compositional index string to test.

    Returns:
        A boolean for a single test cell, or a positionally aligned list.
    """
    raise NotImplementedError


# --------------------------------------------------------------------------
# Compaction
# --------------------------------------------------------------------------


def compact_cells(index: str) -> str:
    """Replace exhaustive sibling sets by their parent, recursively.

    Runs to a fixed point, so a fully covered base cell collapses all the
    way to resolution 1.

    Distinct from :func:`itacart.index.normalize`: compaction changes
    which resolutions appear, while normalization only fixes spelling.
    Normalization already applies the completeness collapse, so this is
    the same reduction exposed under the name callers coming from H3
    expect.

    Args:
        index: Compositional index string.

    Returns:
        The compacted index, mixed-resolution in general.
    """
    raise NotImplementedError


def uncompact_cells(index: str, target_res: int) -> Iterator[str]:
    """Expand a compacted index back to uniform resolution.

    Args:
        index: Compositional index string, possibly mixed-resolution.
        target_res: Resolution to expand every cell to.

    Yields:
        Atomic index strings at ``target_res``.

    Raises:
        ResolutionError: If any terminal cell is finer than ``target_res``.
    """
    raise NotImplementedError
