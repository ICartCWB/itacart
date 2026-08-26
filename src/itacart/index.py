"""The compositional hierarchical index: parsing, composition, canonical form.

A single index string may denote one terminal cell or a whole region::

    SE(1400/0374(3(C2(3))))          one cell
    NW(0625/0451(1(E1(3(B2(4(A2,B2))))))) two sibling cells under one parent

Parentheses descend a level, commas separate siblings at the same level.
Every public operation in the package accepts this form; operations that
return one value per cell return a list aligned with
:func:`decompose` order.

Satisfies OGC DGGS Core requirement 13 (unique address): uniqueness holds
only against the canonical form produced by :func:`normalize`, since the
same region admits several spellings (``4(1,2,3,4)`` and ``4`` denote the
same space).

Origem: itacart_core/compositional_index.py + parser do notebook DGGS_Tree.

Grammar accepted by :func:`parse` (``D-2.2`` widens the first production to
a list, so that :func:`compose` over cells of several quadrants produces a
string this module can read back)::

    index      := root (',' root)*
    root       := QUADRANT subtree?
    subtree    := '(' nodelist ')'
    nodelist   := node (',' node)*
    node       := component subtree?
    component  := base | refinement
    base       := DIGITS '/' DIGITS          (resolution 1 only)
    refinement := a code of refinement_alphabet(resolution)

The even/odd rule is deliberately absent from that last production: it
lives in :func:`itacart.constants.refinement_alphabet` and nowhere else
(``D-0.11``). Writing it a second time here is exactly how ``B-0.1``
happened in the prototype.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final, Iterable, Iterator

from .constants import (
    DESCENT_CLOSE,
    DESCENT_OPEN,
    MAX_RESOLUTION,
    MIN_RESOLUTION,
    QUADRANT_CODE_LENGTH,
    QUADRANTS,
    RES1_DIGITS,
    RES1_MAX_INDEX,
    RES1_SEPARATOR,
    SIBLING_SEPARATOR,
    refinement_alphabet,
)
from .exceptions import (
    InvalidIndexError,
    InvalidQuadrantError,
    InvalidRefinementCodeError,
    NonAtomicIndexError,
)

__all__ = [
    "parse",
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
]

# --------------------------------------------------------------------------
# Levels
# --------------------------------------------------------------------------

QUADRANT_RESOLUTION: Final[int] = MIN_RESOLUTION
"""Resolution addressed by a quadrant code, namely zero."""

BASE_CELL_RESOLUTION: Final[int] = MIN_RESOLUTION + 1
"""Resolution addressed by an ``X/Y`` pair, namely one."""

GLOBE_RESOLUTION: Final[int] = MIN_RESOLUTION - 1
"""Resolution of the synthetic root of a parsed tree, one level up.

Its value is -1. The four quadrants sit at resolution 0, so their common
parent is the globe, one level coarser. It is not addressable — no code names it — but
having it lets a walker treat every real node uniformly, and lets one tree
carry the several quadrant roots that :func:`compose` may produce.
"""

# --------------------------------------------------------------------------
# Alphabets, derived from constants and never restated here
# --------------------------------------------------------------------------

_REFINEMENT_LEVELS: Final[range] = range(BASE_CELL_RESOLUTION + 1, MAX_RESOLUTION + 1)

_ALPHABET_SET: Final[dict[int, frozenset[str]]] = {
    resolution: frozenset(refinement_alphabet(resolution))
    for resolution in _REFINEMENT_LEVELS
}
"""Membership test per resolution, derived from ``REFINEMENT_ALPHABET``."""

_ALPHABET_RANK: Final[dict[int, dict[str, int]]] = {
    resolution: {
        code: position for position, code in enumerate(refinement_alphabet(resolution))
    }
    for resolution in _REFINEMENT_LEVELS
}
"""Position of each code inside its alphabet, used to order siblings.

Ranking by position rather than by the string itself keeps the ordering
stable if an alphabet ever grows past one digit, where ``"A10" < "A2"``
lexically.
"""

_RES1_MAX_X: Final[int] = int(RES1_MAX_INDEX.split(RES1_SEPARATOR)[0])
_RES1_MAX_Y: Final[int] = int(RES1_MAX_INDEX.split(RES1_SEPARATOR)[1])

_DELIMITERS: Final[frozenset[str]] = frozenset(
    {DESCENT_OPEN, DESCENT_CLOSE, SIBLING_SEPARATOR}
)


# --------------------------------------------------------------------------
# Internal tree
# --------------------------------------------------------------------------


@dataclass(slots=True)
class _Node:
    """One node of a parsed index tree.

    ``code`` is the empty string on the synthetic globe root and a real
    component everywhere else. A node with no children denotes the whole
    cell at its own resolution.
    """

    code: str
    resolution: int
    children: list[_Node] = field(default_factory=list)


# --------------------------------------------------------------------------
# Component validation
# --------------------------------------------------------------------------


def _validate_base(component: str) -> tuple[int, int]:
    """Validate a resolution-1 ``X/Y`` component and return ``(x, y)``."""
    if RES1_SEPARATOR not in component:
        raise InvalidIndexError(
            f"resolution-1 component must be 'X{RES1_SEPARATOR}Y', "
            f"got {component!r}"
        )
    x_str, _, y_str = component.partition(RES1_SEPARATOR)
    for part in (x_str, y_str):
        # str.isdigit() alone accepts superscripts and other Unicode digits
        # that int() then refuses; the ASCII guard keeps the two agreeing.
        if not part or not part.isascii() or not part.isdigit():
            raise InvalidIndexError(
                f"resolution-1 indices must be decimal integers, got {component!r}"
            )
    x, y = int(x_str), int(y_str)
    if not 0 <= x <= _RES1_MAX_X:
        raise InvalidIndexError(f"X index {x} outside [0, {_RES1_MAX_X}]")
    if not 0 <= y <= _RES1_MAX_Y:
        raise InvalidIndexError(f"Y index {y} outside [0, {_RES1_MAX_Y}]")
    return x, y


def _validate_refinement(component: str, resolution: int) -> None:
    """Validate a refinement code against the alphabet of its resolution."""
    if component not in _ALPHABET_SET[resolution]:
        alphabet = refinement_alphabet(resolution)
        raise InvalidRefinementCodeError(
            f"{component!r} is not a resolution-{resolution} refinement code: "
            f"expected one of {alphabet[0]}..{alphabet[-1]}"
        )


def _validate_component(component: str, resolution: int) -> None:
    if resolution == BASE_CELL_RESOLUTION:
        _validate_base(component)
    else:
        _validate_refinement(component, resolution)


def _canonical_base(component: str) -> str:
    """Zero-pad an ``X/Y`` component to the canonical width."""
    x, y = _validate_base(component)
    return f"{x:0{RES1_DIGITS}d}{RES1_SEPARATOR}{y:0{RES1_DIGITS}d}"


# --------------------------------------------------------------------------
# Recursive-descent parser
# --------------------------------------------------------------------------


class _Parser:
    """Recursive-descent parser over the grammar in the module docstring."""

    __slots__ = ("_i", "_n", "_s")

    def __init__(self, text: str) -> None:
        self._s = text
        self._i = 0
        self._n = len(text)

    def _peek(self) -> str | None:
        return self._s[self._i] if self._i < self._n else None

    def _expect(self, char: str) -> None:
        if self._peek() != char:
            raise InvalidIndexError(
                f"expected {char!r} at position {self._i} in {self._s!r}"
            )
        self._i += 1

    def _read_component(self) -> str:
        start = self._i
        while self._i < self._n and self._s[self._i] not in _DELIMITERS:
            self._i += 1
        end = self._i
        if end == start:
            raise InvalidIndexError(
                f"empty component at position {start} in {self._s!r}"
            )
        return self._s[start:end]

    def _parse_nodelist(self, parent: _Node, resolution: int) -> None:
        if resolution > MAX_RESOLUTION:
            raise InvalidIndexError(
                f"index descends past resolution {MAX_RESOLUTION} "
                f"at position {self._i} in {self._s!r}"
            )
        while True:
            component = self._read_component()
            _validate_component(component, resolution)
            node = _Node(component, resolution)
            parent.children.append(node)
            if self._peek() == DESCENT_OPEN:
                self._i += 1
                self._parse_nodelist(node, resolution + 1)
                self._expect(DESCENT_CLOSE)
            if self._peek() == SIBLING_SEPARATOR:
                self._i += 1
                continue
            break

    def _parse_root(self) -> _Node:
        start = self._i
        end = start + QUADRANT_CODE_LENGTH
        quadrant = self._s[start:end]
        if quadrant not in QUADRANTS:
            raise InvalidQuadrantError(
                f"{quadrant!r} is not a quadrant code: expected one of "
                f"{', '.join(QUADRANTS)}"
            )
        self._i += QUADRANT_CODE_LENGTH
        node = _Node(quadrant, QUADRANT_RESOLUTION)
        if self._peek() == DESCENT_OPEN:
            self._i += 1
            self._parse_nodelist(node, BASE_CELL_RESOLUTION)
            self._expect(DESCENT_CLOSE)
        return node

    def parse(self) -> _Node:
        if self._n == 0:
            raise InvalidIndexError("index string is empty")
        globe = _Node("", GLOBE_RESOLUTION)
        while True:
            globe.children.append(self._parse_root())
            if self._peek() == SIBLING_SEPARATOR:
                self._i += 1
                continue
            break
        if self._i != self._n:
            raise InvalidIndexError(
                f"trailing characters from position {self._i} in {self._s!r}"
            )
        return globe


def _parse_tree(index: str) -> _Node:
    """Parse into the internal tree, the substrate every operation walks."""
    if not isinstance(index, str):
        raise InvalidIndexError(f"index must be a string, got {type(index).__name__}")
    return _Parser(index.strip()).parse()


# --------------------------------------------------------------------------
# Tree walking
# --------------------------------------------------------------------------


def _leaf_paths(node: _Node, prefix: tuple[str, ...] = ()) -> Iterator[tuple[str, ...]]:
    """Yield the component path of every terminal cell, depth first."""
    if node.resolution == GLOBE_RESOLUTION:
        for child in node.children:
            yield from _leaf_paths(child, prefix)
        return
    path = prefix + (node.code,)
    if not node.children:
        yield path
        return
    for child in node.children:
        yield from _leaf_paths(child, path)


def _count_leaves(node: _Node) -> int:
    if node.resolution != GLOBE_RESOLUTION and not node.children:
        return 1
    return sum(_count_leaves(child) for child in node.children)


def _render_path(path: tuple[str, ...]) -> str:
    """Render one component path as a fully qualified atomic index."""
    quadrant, *rest = path
    if not rest:
        return quadrant
    body = DESCENT_OPEN.join(rest)
    return f"{quadrant}{DESCENT_OPEN}{body}{DESCENT_CLOSE * len(rest)}"


def _render_node(node: _Node) -> str:
    if not node.children:
        return node.code
    inner = SIBLING_SEPARATOR.join(_render_node(child) for child in node.children)
    return f"{node.code}{DESCENT_OPEN}{inner}{DESCENT_CLOSE}"


def _render_tree(globe: _Node) -> str:
    return SIBLING_SEPARATOR.join(_render_node(child) for child in globe.children)


# --------------------------------------------------------------------------
# Merging: how shared ancestry is factored out
# --------------------------------------------------------------------------


def _merged_copy(node: _Node) -> _Node:
    fresh = _Node(node.code, node.resolution)
    for child in node.children:
        _merge_into(fresh.children, child)
    return fresh


def _merge_into(siblings: list[_Node], node: _Node) -> None:
    """Fold ``node`` into ``siblings``, merging with an equal code if any.

    A childless occurrence wins: a node written without a subtree denotes
    the whole cell, which already contains whatever its other occurrence
    descends into.
    """
    for existing in siblings:
        if existing.code != node.code:
            continue
        if not existing.children or not node.children:
            existing.children = []
            return
        for child in node.children:
            _merge_into(existing.children, child)
        return
    siblings.append(_merged_copy(node))


def _merge_trees(trees: Iterable[_Node]) -> _Node:
    globe = _Node("", GLOBE_RESOLUTION)
    for tree in trees:
        for quadrant_node in tree.children:
            _merge_into(globe.children, quadrant_node)
    return globe


# --------------------------------------------------------------------------
# Canonicalisation
# --------------------------------------------------------------------------


def _pad_base_codes(node: _Node) -> None:
    """Rewrite resolution-1 codes to the canonical zero-padded width."""
    if node.resolution == BASE_CELL_RESOLUTION:
        node.code = _canonical_base(node.code)
    for child in node.children:
        _pad_base_codes(child)


def _collapse(node: _Node) -> None:
    """Replace a fully covered node by itself, deepest level first.

    Only refinement levels collapse. A quadrant holding every base cell is
    left alone: resolution 1 is addressed by the 2 004 x 1 001 product of
    Table 1, not by a refinement alphabet, and materialising it to test
    completeness would cost two million comparisons per node.
    """
    for child in node.children:
        _collapse(child)
    if not node.children or node.resolution < BASE_CELL_RESOLUTION:
        return
    # The node has children, so their resolution is addressed by a
    # refinement alphabet: the parser refuses to descend past
    # MAX_RESOLUTION, which is why this lookup needs no guard.
    whole_children = {child.code for child in node.children if not child.children}
    if whole_children == _ALPHABET_SET[node.resolution + 1]:
        node.children = []


def _sibling_key(node: _Node) -> tuple[int, int]:
    """Ordering key for siblings, by position in the level's alphabet."""
    if node.resolution == QUADRANT_RESOLUTION:
        return (QUADRANTS.index(node.code), 0)
    if node.resolution == BASE_CELL_RESOLUTION:
        x_str, _, y_str = node.code.partition(RES1_SEPARATOR)
        return (int(x_str), int(y_str))
    return (_ALPHABET_RANK[node.resolution][node.code], 0)


def _sort_siblings(node: _Node) -> None:
    for child in node.children:
        _sort_siblings(child)
    node.children.sort(key=_sibling_key)


# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------


def parse(index: str) -> dict[str, object]:
    """Parse an index string into a nested structural tree.

    The tree is the shared substrate for hierarchy, topology and
    serialization; most callers want :func:`decompose` instead.

    Every node carries ``code``, ``resolution`` and ``children``. The root
    is the globe — ``resolution`` -1 and ``code`` ``None`` — and carries
    two further keys: ``quadrant`` and ``base_cell``, each holding the one
    code when the index has exactly one and ``None`` when it addresses
    several. ``base_cell`` is the bare ``X/Y`` component, not the qualified
    form that :func:`base_cell_of` returns.

    Args:
        index: Compositional index string.

    Returns:
        A nested mapping with ``quadrant``, ``base_cell`` and ``children``.

    Raises:
        InvalidIndexError: On malformed syntax, unbalanced parentheses, or
            a refinement code outside the alphabet of its resolution.
    """
    return _to_mapping(_parse_tree(index))


def _to_mapping(node: _Node) -> dict[str, object]:
    is_globe = node.resolution == GLOBE_RESOLUTION
    mapping: dict[str, object] = {
        "code": None if is_globe else node.code,
        "resolution": node.resolution,
        "children": [_to_mapping(child) for child in node.children],
    }
    if is_globe:
        roots = node.children
        quadrant = roots[0] if len(roots) == 1 else None
        mapping["quadrant"] = None if quadrant is None else quadrant.code
        base_cell = None
        if quadrant is not None and len(quadrant.children) == 1:
            base_cell = quadrant.children[0].code
        mapping["base_cell"] = base_cell
    return mapping


def is_valid_index(index: str) -> bool:
    """Whether the string parses and every code is legal for its level.

    Structural validity only. It does not check whether the addressed
    cells actually exist in the ITACaRT domain; for that see
    :func:`itacart.boundary.is_valid_cell`, which also applies the
    western-quadrant ``X = 0`` rule.

    Args:
        index: Candidate index string.

    Returns:
        ``True`` if the string is syntactically well formed.
    """
    try:
        _parse_tree(index)
    except InvalidIndexError:
        return False
    return True


def is_atomic(index: str) -> bool:
    """Whether the index denotes exactly one terminal cell.

    Args:
        index: Compositional index string.

    Returns:
        ``True`` when the index has no sibling separators.

    Raises:
        InvalidIndexError: If the index is malformed. A malformed string
            denotes no cell at all, which is a different claim from "more
            than one"; :func:`is_valid_index` is the predicate for that.
    """
    return _count_leaves(_parse_tree(index)) == 1


def decompose(index: str) -> list[str]:
    """Expand a compositional index into a flat list of atomic indices.

    Order is deterministic: depth-first, left to right, matching the order
    siblings appear in the string. Every vectorised operation in the
    package aligns its output to this order.

    Args:
        index: Compositional index string.

    Returns:
        One fully qualified atomic index per terminal cell.

    Example:
        ``"NE(0001/0002(1(A1,A2)))"`` yields
        ``["NE(0001/0002(1(A1)))", "NE(0001/0002(1(A2)))"]``.
    """
    return list(iter_cells(index))


def compose(cells: Iterable[str]) -> str:
    """Fold a list of atomic indices back into one compositional index.

    Inverse of :func:`decompose`. Factors out shared ancestry so siblings
    collapse under their common parent. Cells may sit at different
    resolutions and in different quadrants; the result then carries
    several roots separated by commas.

    Siblings keep the order in which they first appear, and no completeness
    collapse is applied: rewriting spelling is :func:`normalize`'s job, and
    keeping the two apart is what makes ``compose(decompose(x)) == x`` an
    exact identity rather than an approximate one.

    Args:
        cells: Index strings. Atomic in the usual case, but a compositional
            entry is accepted and merged like any other.

    Returns:
        A single compositional index string.

    Raises:
        InvalidIndexError: If any entry is malformed, or the iterable is
            empty — no index denotes the empty region.
    """
    trees = [_parse_tree(cell) for cell in cells]
    if not trees:
        raise InvalidIndexError("cannot compose an empty collection of cells")
    return _render_tree(_merge_trees(trees))


def normalize(index: str) -> str:
    """Reduce an index to its canonical form.

    Two reductions are applied to a fixed point:

    1. **Completeness collapse** - a node whose children are all present
       is replaced by the node itself (``4(1,2,3,4)`` becomes ``4``).
    2. **Sibling ordering** - siblings are sorted by their refinement
       alphabet so a region has one spelling regardless of input order.

    Two further rewrites are needed for the canonical form to be unique in
    practice: repeated siblings are merged, since ``4(1),4(2)`` and
    ``4(1,2)`` are one region, and resolution-1 components are zero-padded
    to four digits, since ``1400/374`` and ``1400/0374`` are one cell.

    This is what makes OGC requirement 13 hold in practice: two indices
    denote the same region if and only if their canonical forms are equal.

    Origem: analogo ASCII de binary_index.recompose_to_prefix_form.

    Args:
        index: Compositional index string.

    Returns:
        The canonical index string.
    """
    tree = _parse_tree(index)
    _pad_base_codes(tree)
    merged = _merge_trees([tree])
    _collapse(merged)
    _sort_siblings(merged)
    return _render_tree(merged)


def count_cells(index: str) -> int:
    """Number of terminal cells addressed, without materialising them.

    Args:
        index: Compositional index string.

    Returns:
        Count of terminal cells.
    """
    return _count_leaves(_parse_tree(index))


def iter_cells(index: str) -> Iterator[str]:
    """Stream atomic indices in :func:`decompose` order.

    Preferable to :func:`decompose` for dense regions, where the flat list
    can be very large.

    Args:
        index: Compositional index string.

    Yields:
        Fully qualified atomic index strings.
    """
    for path in _leaf_paths(_parse_tree(index)):
        yield _render_path(path)


def split_components(cell: str) -> list[str]:
    """Break an atomic index into its per-level components.

    Args:
        cell: Atomic index string.

    Returns:
        Components ordered coarse to fine, starting with the quadrant and
        the base cell, e.g. ``["SE", "1400/0374", "3", "C2", "3"]``.

    Raises:
        NonAtomicIndexError: If the index holds more than one cell.
    """
    paths = list(_leaf_paths(_parse_tree(cell)))
    if len(paths) != 1:
        raise NonAtomicIndexError(f"{cell!r} addresses {len(paths)} cells, not one")
    return list(paths[0])


def quadrant_of(index: str) -> str | list[str]:
    """Global quadrant of the addressed cells.

    Args:
        index: Compositional index string.

    Returns:
        A quadrant code for a single cell, or a positionally aligned list.
    """
    return _ancestor_per_cell(index, 0)


def base_cell_of(index: str) -> str | list[str]:
    """Resolution-1 ancestor of the addressed cells.

    Args:
        index: Compositional index string.

    Returns:
        A base cell index for a single cell, or a positionally aligned
        list. Each entry is fully qualified, e.g. ``"SE(1400/0374)"``.

    Raises:
        InvalidIndexError: If any addressed cell is a whole quadrant, which
            has no resolution-1 ancestor.
    """
    return _ancestor_per_cell(index, 1)


def _ancestor_per_cell(index: str, depth: int) -> str | list[str]:
    """Ancestor at ``depth`` components: scalar for one cell, list for many.

    ``depth`` 0 is the quadrant, 1 the base cell. The scalar return follows
    the contract of the two public callers rather than ``D-0.2``: whoever
    asks the quadrant of one cell wants a code, not a list of one.
    """
    values: list[str] = []
    for path in _leaf_paths(_parse_tree(index)):
        if len(path) <= depth:
            raise InvalidIndexError(
                f"{_render_path(path)!r} has no component at depth {depth}"
            )
        values.append(_render_path(path[: depth + 1]))
    if len(values) == 1:
        return values[0]
    return values
