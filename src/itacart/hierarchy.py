"""Hierarchical navigation: ascent, descent, compaction.

Part of OGC DGGS Core requirement 17 (topological query functions), with
the neighbour half living in :mod:`itacart.topology`.

Every function here accepts a compositional index. Where the input holds
several cells, the return is positionally aligned with
:func:`itacart.index.decompose`.

Ascent and descent are not symmetric, and the asymmetry is the whole
difficulty of this module. Ascent is a string operation: the parent is a
prefix of the child, so it costs a slice and never consults the domain.
Descent is not, because a cell whose outer side has been carried onto the
domain border reaches east past its own nominal column, and some of its
children are therefore spelled under the next column's resolution-1
prefix. Those cells exist, they are addressed uniquely, and there are
neither four nor twenty-five of them. Ascent from such a child lands on a
string that is a well-formed prefix but names no cell.

The module answers this by keeping the two claims apart. Functions that
ascend return prefixes and say so; :func:`itacart.boundary.is_valid_cell`
is the only arbiter of whether a prefix names a cell. Functions that
descend consult :func:`itacart.boundary.absorbs_border` first and take a
purely lexical path when it answers false, which it does for every cell
of the lattice except the last column of each row.
"""

from __future__ import annotations

from typing import Iterable, Iterator, Sequence, cast

from .constants import (
    MAX_RESOLUTION,
    MIN_RESOLUTION,
    QUADRANTS,
    RES1_DIGITS,
    RES1_SEPARATOR,
    refinement_alphabet,
)
from .exceptions import (
    DomainError,
    GeometryError,
    MaxResolutionError,
    MinResolutionError,
    NonExistentCellError,
    ResolutionError,
)
from .index import (
    BASE_CELL_RESOLUTION,
    QUADRANT_RESOLUTION,
    compose,
    decompose,
    split_components,
)

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


_DESCENT_OPEN = "("
_DESCENT_CLOSE = ")"

_OVERLAP_EPSILON_RATIO = 1e-6
"""Fraction of a child's nominal area below which overlap is only contact.

Border children are selected by intersecting plane rings, and two cells
that merely share an edge intersect in a sliver of rounding noise. The
threshold has to sit above that sliver and below the smallest genuine
overlap. Both of those scale with the cell, so the threshold scales with
it too, and the constant is a ratio rather than an area.

Measured, not assumed. The sliver between lexical neighbours is exactly
zero up to level 11 and at most 2.33e-9 of the nominal area at levels 12
and 13. The smallest genuine overlap, enumerated over all 4000
border-absorbing cells of resolution 1 and over 60 chains carried down
to resolution 13, is 2.05e-2 of the nominal area. This ratio sits 430
times above the first and 20500 times below the second, close to the
geometric mean of the two.

A fixed area cannot do this job. One square metre, which is what stood
here, is the entire nominal area of a resolution-9 cell and is larger
than every cell below it, so a child wholly inside its parent was
rejected for being the size it is supposed to be and the eastern border
of the grid had no refinement at all below resolution 8.
"""


# --------------------------------------------------------------------------
# Component paths
# --------------------------------------------------------------------------


def _render(components: Sequence[str]) -> str:
    """Render a component path back into a fully qualified index string."""
    quadrant, *rest = components
    if not rest:
        return quadrant
    body = _DESCENT_OPEN.join(rest)
    return f"{quadrant}{_DESCENT_OPEN}{body}{_DESCENT_CLOSE * len(rest)}"


def _resolution_of(components: Sequence[str]) -> int:
    """Resolution addressed by a component path.

    One component is a bare quadrant at resolution zero, two reach the
    resolution-1 lattice, and each further component descends one level.
    """
    return len(components) - 1


def _sort_key(components: Sequence[str]) -> tuple[object, ...]:
    """Deterministic ordering key for a component path.

    Quadrants order as they are declared, resolution-1 cells by the
    numeric pair, and refinements by their rank in the level alphabet.
    """
    key: list[object] = [QUADRANTS.index(components[0])]
    for position, component in enumerate(components[1:], start=BASE_CELL_RESOLUTION):
        if position == BASE_CELL_RESOLUTION:
            column, _, row = component.partition(RES1_SEPARATOR)
            key.append((int(column), int(row)))
        else:
            key.append((refinement_alphabet(position).index(component), 0))
    return tuple(key)


def _paths(index: str) -> list[list[str]]:
    """Component path of every cell addressed by an index, in write order."""
    return [split_components(cell) for cell in decompose(index)]


def _is_single(index: str) -> bool:
    """Whether an index addresses exactly one cell."""
    return len(decompose(index)) == 1


def _scalar_or_list(index: str, values: list[str]) -> str | list[str]:
    """Answer a scalar for a single-cell index and a list otherwise.

    Follows the convention of :func:`itacart.index.quadrant_of`: whoever
    asks about one cell wants one answer, not a list holding one.
    """
    if len(values) == 1 and _is_single(index):
        return values[0]
    return values


# --------------------------------------------------------------------------
# Ascent
# --------------------------------------------------------------------------


def get_parent(index: str, target_res: int | None = None) -> str | list[str]:
    """Ascend the hierarchy by lexical truncation.

    No floating-point work and no negative indices: the parent is a
    prefix of the child string, which is the property that makes the
    compositional index cheap to navigate.

    The return is a **prefix**, not necessarily a cell. For most of the
    grid the two coincide, but a cell that absorbs the domain border
    reaches past its own column and fathers children spelled under the
    next column's prefix; truncating such a child yields a well-formed
    string that names no cell. Callers that need the stronger claim pass
    the result to :func:`itacart.boundary.is_valid_cell`.

    Args:
        index: Compositional index string.
        target_res: Resolution to ascend to. Defaults to one level up.

    Returns:
        The ancestor prefix for a single cell, or a positionally aligned
        list. Duplicates are preserved so alignment holds; call
        :func:`itacart.index.normalize` on the composed result to dedupe.

    Raises:
        MinResolutionError: If ``target_res`` is above resolution 0.
        ResolutionError: If ``target_res`` is finer than the input.
    """
    if target_res is not None and target_res < MIN_RESOLUTION:
        raise MinResolutionError(
            f"resolution {target_res} is above the quadrant level {MIN_RESOLUTION}"
        )
    values: list[str] = []
    for components in _paths(index):
        current = _resolution_of(components)
        wanted = current - 1 if target_res is None else target_res
        if wanted < MIN_RESOLUTION:
            raise MinResolutionError(
                f"{_render(components)!r} is at resolution {current} "
                "and has no coarser prefix"
            )
        if wanted > current:
            raise ResolutionError(
                f"resolution {wanted} is finer than {_render(components)!r} "
                f"at resolution {current}"
            )
        values.append(_render(components[: wanted + 1]))
    return _scalar_or_list(index, values)


def get_ancestors(index: str) -> list[str] | list[list[str]]:
    """Every ancestor of a cell, coarse to fine.

    Strict: a cell is not among its own ancestors. The reflexive relation
    is :func:`contains`.

    Args:
        index: Compositional index string.

    Returns:
        An ancestor chain for a single cell, or a positionally aligned
        list of chains. Entries are prefixes under the same caveat as
        :func:`get_parent`.
    """
    chains: list[list[str]] = []
    for components in _paths(index):
        chains.append(
            [_render(components[: depth + 1]) for depth in range(len(components) - 1)]
        )
    if len(chains) == 1 and _is_single(index):
        return chains[0]
    return cast("list[str] | list[list[str]]", chains)


def common_ancestor(index: str) -> str:
    """Deepest cell containing every terminal cell of the index.

    Underlies :func:`compact_cells` and :func:`itacart.index.compose`, and
    is independently useful for finding the grouping node of a vertex set.

    Args:
        index: Compositional index string, or a composed list of cells.

    Returns:
        The deepest common ancestor index. Returns the quadrant when the
        cells share only that.

    Raises:
        DomainError: If the cells span more than one quadrant.
    """
    paths = _paths(index)
    quadrants = {components[0] for components in paths}
    if len(quadrants) > 1:
        raise DomainError(
            f"{index!r} spans quadrants {sorted(quadrants)} "
            "and has no common ancestor"
        )
    shared = list(paths[0])
    for components in paths[1:]:
        depth = 0
        limit = min(len(shared), len(components))
        while depth < limit and shared[depth] == components[depth]:
            depth += 1
        shared = shared[:depth]
    return _render(shared)


# --------------------------------------------------------------------------
# Descent
# --------------------------------------------------------------------------


def _descend(cell: str, code: str) -> str:
    """Append one refinement code to an atomic index string."""
    depth = cell.count(_DESCENT_OPEN)
    head = cell[: len(cell) - depth]
    return head + _DESCENT_OPEN + code + _DESCENT_CLOSE * (depth + 1)


def _shift_column(cell: str, step: int) -> str:
    """The same index with its resolution-1 column moved east by ``step``."""
    quadrant, _, rest = cell.partition(_DESCENT_OPEN)
    column, _, tail = rest.partition(RES1_SEPARATOR)
    shifted = f"{int(column) + step:0{RES1_DIGITS}d}"
    return f"{quadrant}{_DESCENT_OPEN}{shifted}{RES1_SEPARATOR}{tail}"


def _children_of(cell: str) -> list[str]:
    """Every existing cell one resolution below ``cell``, in canonical order.

    Two regimes. A cell that does not absorb the domain border refines
    into exactly the level alphabet, and the children are read off the
    string with no geometry at all. A cell that does absorb it reaches
    east past its own column, so both its own stem and that of the column
    immediately east are descended, and a candidate is kept when its
    polygon actually lies inside the parent. Descending only the own stem
    loses children.

    The order is the enumeration order: own stem before eastern stem, and
    within a stem the order of the level alphabet. This agrees with
    sorting the children by component path, because the eastern stem
    carries the larger column.
    """
    from . import boundary

    components = split_components(cell)
    current = _resolution_of(components)
    if current < BASE_CELL_RESOLUTION:
        raise ResolutionError(
            f"{cell!r} is a quadrant; its children are the resolution-1 "
            "lattice, addressed by a coordinate pair rather than by a "
            "refinement alphabet, and are not enumerated here"
        )
    if current >= MAX_RESOLUTION:
        raise MaxResolutionError(
            f"{cell!r} is at resolution {MAX_RESOLUTION} and has no children"
        )
    level = current + 1
    if not boundary.absorbs_border(cell):
        return [_descend(cell, code) for code in refinement_alphabet(level)]
    return _border_children_of(cell, level)


def _border_children_of(cell: str, level: int) -> list[str]:
    """Children of a border-absorbing cell, selected geometrically.

    A candidate whose ring is self-intersecting is refused loudly rather
    than repaired. Repairing it would answer with a child set derived
    from a polygon nobody meant to draw, and the caller would have no way
    to tell. The one place this fires is the polar triangle, whose
    refinement rings are malformed upstream.
    """
    from shapely.geometry import Polygon

    from . import boundary
    from .resolutions import cell_size

    body = Polygon(boundary.plane_ring(cell)[1])
    epsilon = _OVERLAP_EPSILON_RATIO * cell_size(level) ** 2
    children: list[str] = []
    for step in (0, 1):
        stem = _shift_column(cell, step)
        for code in refinement_alphabet(level):
            candidate = _descend(stem, code)
            if not boundary.is_valid_cell(candidate):
                continue
            outline = Polygon(boundary.plane_ring(candidate)[1])
            if not outline.is_valid:
                raise GeometryError(
                    f"the refinement ring of {candidate!r} is self-intersecting, "
                    f"so the children of {cell!r} cannot be selected by overlap"
                )
            if outline.intersection(body).area > epsilon:
                children.append(candidate)
    return children


def get_children(
    index: str, target_res: int | None = None, flatten: bool = False
) -> Iterator[str] | Iterator[list[str]]:
    """Project descendants of an index down to a target resolution.

    Yields 4 children when descending into an even resolution and 25 when
    descending into an odd one — away from the domain border. A cell
    whose outer side has been carried onto the border was measured to
    hold between two and six children instead, so the count comes from
    enumeration and never from
    :func:`itacart.resolutions.refinement_ratio`.

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
        ResolutionError: If ``target_res`` is not finer than the input.
    """
    if target_res is not None and target_res > MAX_RESOLUTION:
        raise MaxResolutionError(
            f"resolution {target_res} is finer than the maximum {MAX_RESOLUTION}"
        )
    plans: list[tuple[str, int]] = []
    for cell in decompose(index):
        current = _resolution_of(split_components(cell))
        wanted = current + 1 if target_res is None else target_res
        if wanted <= current:
            raise ResolutionError(
                f"resolution {wanted} is not finer than {cell!r} "
                f"at resolution {current}"
            )
        plans.append((cell, wanted))
    if flatten:
        return _flat_children(plans)
    return _grouped_children(plans)


def _flat_children(plans: Iterable[tuple[str, int]]) -> Iterator[str]:
    for cell, wanted in plans:
        yield from _descendants_of(cell, wanted)


def _grouped_children(plans: Iterable[tuple[str, int]]) -> Iterator[list[str]]:
    for cell, wanted in plans:
        yield list(_descendants_of(cell, wanted))


def _descendants_of(cell: str, target_res: int) -> Iterator[str]:
    """Stream descendants of one atomic cell, depth first."""
    for child in _children_of(cell):
        if _resolution_of(split_components(child)) == target_res:
            yield child
        else:
            yield from _descendants_of(child, target_res)


def get_descendants(index: str, target_res: int) -> Iterator[str]:
    """Stream every descendant at a target resolution.

    Cardinality grows fast: from resolution 1 to 13 a base cell away from
    the border undergoes twelve refinements, six of ratio 4 and six of
    ratio 25, hence ``(4 * 25) ** 6`` — one million million cells. Always
    a generator, never a list.

    Args:
        index: Compositional index string.
        target_res: Resolution to expand to.

    Yields:
        Atomic index strings at ``target_res``.

    Raises:
        MaxResolutionError: If ``target_res`` exceeds resolution 13.
        ResolutionError: If ``target_res`` is not finer than any input cell.
    """
    return cast("Iterator[str]", get_children(index, target_res, flatten=True))


def _parent_cell(cell: str) -> str:
    """The cell that actually fathers ``cell``, border cases included.

    The lexical prefix is the answer everywhere except under a
    border-absorbing parent, where the child may be spelled under the
    column immediately east. That column never holds a cell of its own in
    that row — it is exactly one past the last one — so a prefix that
    fails :func:`itacart.boundary.is_valid_cell` is the signal to step one
    column west.
    """
    from . import boundary

    components = split_components(cell)
    if _resolution_of(components) <= QUADRANT_RESOLUTION:
        raise MinResolutionError(f"{cell!r} is a quadrant and has no parent cell")
    prefix = _render(components[:-1])
    if boundary.is_valid_cell(prefix):
        return prefix
    western = _shift_column(prefix, -1)
    if boundary.is_valid_cell(western) and cell in _children_of(western):
        return western
    raise NonExistentCellError(f"{cell!r} has no parent cell in the domain")


def child_position(cell: str) -> int | list[int]:
    """Ordinal of a cell among its siblings.

    Zero-based against the refinement alphabet of the cell's resolution:
    ``1``-``4`` for even levels, ``A1``-``E5`` in row-major order for odd.

    Under a border-absorbing parent the alphabet rank is not an ordinal,
    because two siblings can carry the same code under different column
    prefixes. There the ordinal is the position in the enumeration of
    :func:`get_children`: own stem before eastern stem, alphabet order
    within each. Away from the border the two definitions agree, and the
    cheap one is used.

    Args:
        cell: Compositional index string.

    Returns:
        An ordinal for a single cell, or a positionally aligned list.

    Raises:
        MinResolutionError: If any addressed cell is a whole quadrant.
        ResolutionError: If any addressed cell is a resolution-1 cell.
        NonExistentCellError: If any addressed cell has no parent cell.
    """
    from . import boundary

    positions: list[int] = []
    for atom in decompose(cell):
        components = split_components(atom)
        current = _resolution_of(components)
        if current <= QUADRANT_RESOLUTION:
            raise MinResolutionError(f"{atom!r} is a quadrant and has no siblings")
        if current == BASE_CELL_RESOLUTION:
            raise ResolutionError(
                f"{atom!r} is a resolution-1 cell; its siblings are the "
                "lattice of its quadrant, which carries no refinement alphabet"
            )
        prefix = _render(components[:-1])
        if boundary.is_valid_cell(prefix) and not boundary.absorbs_border(prefix):
            positions.append(refinement_alphabet(current).index(components[-1]))
            continue
        positions.append(_children_of(_parent_cell(atom)).index(atom))
    if len(positions) == 1 and _is_single(cell):
        return positions[0]
    return positions


# --------------------------------------------------------------------------
# Containment
# --------------------------------------------------------------------------


def is_ancestor(parent_index: str, child_index: str) -> bool:
    """Whether one cell contains another in the hierarchy.

    Strict, so a cell is not its own ancestor; the reflexive relation is
    :func:`contains`. A predicate over index strings that never consults
    geometry: the union of a cell's children was measured to overrun the
    cell itself by about one part in a thousand, because the domain
    border is curved and each cell approximates it by a chord, so a
    geometric test would disagree at the border with the index that
    addresses the cell.

    Args:
        parent_index: Candidate ancestor, a single atomic index.
        child_index: Candidate descendant, a single atomic index.

    Returns:
        ``True`` if ``parent_index`` is a strict prefix of
        ``child_index``.
    """
    ancestor = split_components(parent_index)
    descendant = split_components(child_index)
    if len(ancestor) >= len(descendant):
        return False
    return descendant[: len(ancestor)] == ancestor


def contains(region: str, cell: str) -> bool | list[bool]:
    """Whether a compositional region covers a cell.

    Complements :func:`is_ancestor`, which only tests the vertical
    relation between two single cells. A region covers a cell when any of
    its terminal cells is that cell or an ancestor of it.

    Like :func:`is_ancestor`, a predicate over index strings only.

    Args:
        region: Compositional index string denoting the region.
        cell: Compositional index string to test.

    Returns:
        A boolean for a single test cell, or a positionally aligned list.
    """
    terminals = [tuple(components) for components in _paths(region)]
    answers: list[bool] = []
    for candidate in _paths(cell):
        answers.append(
            any(
                len(terminal) <= len(candidate)
                and tuple(candidate[: len(terminal)]) == terminal
                for terminal in terminals
            )
        )
    if len(answers) == 1 and _is_single(cell):
        return answers[0]
    return answers


# --------------------------------------------------------------------------
# Compaction
# --------------------------------------------------------------------------


def _is_complete(parent: str, present: set[tuple[str, ...]]) -> bool:
    """Whether ``present`` holds every child of ``parent``.

    No guard against a childless or over-deep parent: candidates come
    from :func:`_candidate_parents`, which only proposes the prefix of a
    cell, so the parent sits between resolution 1 and 12 and always has
    children. A defective ring still propagates, deliberately.
    """
    children = _children_of(parent)
    return all(tuple(split_components(child)) in present for child in children)


def _candidate_parents(present: set[tuple[str, ...]]) -> set[str]:
    """Cells that could absorb one of the paths in ``present``."""
    from . import boundary

    candidates: set[str] = set()
    for components in present:
        if _resolution_of(components) <= BASE_CELL_RESOLUTION:
            continue
        prefix = _render(components[:-1])
        if boundary.is_valid_cell(prefix):
            candidates.add(prefix)
            continue
        western = _shift_column(prefix, -1)
        if boundary.is_valid_cell(western):
            candidates.add(western)
    return candidates


def compact_cells(index: str) -> str:
    """Replace exhaustive sibling sets by their parent, recursively.

    Runs to a fixed point, so a fully covered base cell collapses all the
    way to resolution 1.

    Distinct from :func:`itacart.index.normalize`, and deliberately
    stricter. Normalization collapses a parent as soon as the level
    alphabet is complete, which overstates coverage under a
    border-absorbing parent: four quaternary children spell the whole
    alphabet, but such a parent may have five. Compaction counts the
    children the parent actually has, so a partition that is incomplete
    at the border is left uncompacted rather than silently claimed whole.
    In short, normalize operates on index syntax; compact_cells operates
    on the spatial semantics of the DGGS.

    Args:
        index: Compositional index string.

    Returns:
        The compacted index, mixed-resolution in general, with cells in
        deterministic order.
    """
    present = {tuple(components) for components in _paths(index)}
    changed = True
    while changed:
        changed = False
        for parent in sorted(_candidate_parents(present), key=len, reverse=True):
            path = tuple(split_components(parent))
            if not _is_complete(parent, present):
                continue
            for child in _children_of(parent):
                present.discard(tuple(split_components(child)))
            present.add(path)
            changed = True
    ordered = sorted(present, key=_sort_key)
    return compose(_render(components) for components in ordered)


def uncompact_cells(index: str, target_res: int) -> Iterator[str]:
    """Expand a compacted index back to uniform resolution.

    Args:
        index: Compositional index string, possibly mixed-resolution.
        target_res: Resolution to expand every cell to.

    Yields:
        Atomic index strings at ``target_res``.

    Raises:
        MaxResolutionError: If ``target_res`` exceeds resolution 13.
        ResolutionError: If any terminal cell is finer than ``target_res``.
    """
    if target_res > MAX_RESOLUTION:
        raise MaxResolutionError(
            f"resolution {target_res} is finer than the maximum {MAX_RESOLUTION}"
        )
    plans: list[tuple[str, bool]] = []
    for cell in decompose(index):
        current = _resolution_of(split_components(cell))
        if current > target_res:
            raise ResolutionError(
                f"{cell!r} is at resolution {current}, finer than the "
                f"requested {target_res}"
            )
        plans.append((cell, current == target_res))
    return _expand(plans, target_res)


def _expand(plans: Iterable[tuple[str, bool]], target_res: int) -> Iterator[str]:
    for cell, already_there in plans:
        if already_there:
            yield cell
        else:
            yield from _descendants_of(cell, target_res)
