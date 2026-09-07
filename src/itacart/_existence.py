"""The existence contract: a spelling naming no cell neither enters nor leaves.

``is_valid_cell`` is the arbiter of existence, and until this module it was
an arbiter nobody consulted. Forty-six public callables accepted a spelling
it denies; thirteen of them answered about the *neighbouring* cell, silently,
because the western zero-column spelling folds onto its eastern twin, and one
returned that denied spelling to the caller inside a cell list.

The contract this module enforces has two halves and one deliberate hole.

Entering, a callable that answers *about a cell* -- its geometry, its shape,
its position, its area, its neighbourhood, its ancestry, its serialized form
-- requires the cell to exist, and refuses with :class:`NonExistentCellError`
when it does not. In a cadastral system, returning the area of a cell that is
not the one asked for is a grave error, which is the reasoning of ``D-0.5``
applied to existence rather than to clipping.

Leaving, no such callable may emit a spelling the predicate denies. That half
needs no runtime check and has none: measured over a corpus of valid cells
spanning four quadrants, nothing in the surface emits a denied spelling, and
the one emission that existed came from a denied spelling going *in*. Closing
the entrance closes the exit by construction. What the exit half gets is a
pin, in ``tests/unit/test_existence.py``, so the property fails loudly if it
ever stops holding.

The hole is intentional and is the older decision this contract does not
reopen. Callables that answer about the *string* rather than about the cell
-- its syntax, its parts, its canonical form, its resolution, its validity --
keep accepting the denied spelling, because analysing a spelling is not the
same act as answering for a cell. Two of them fold it onto the eastern twin
and two hand it straight back. That reach is written down in the exemption
table of the test module rather than left to be rediscovered.

Nothing here polices the interior. A topological computation is free to walk
through a non-canonical lattice coordinate on its way to an answer; the
contract is about the boundary of the public surface, not about the route
taken inside it.
"""

from __future__ import annotations

from typing import Iterable

from .exceptions import NonExistentCellError


def require_existing_cells(index: str | Iterable[str]) -> None:
    """Refuse an index, or any cell of it, that names no cell.

    Syntax is judged first and by the existing machinery. ``decompose``
    raises this package's syntactic exceptions for a malformed string, so a
    caller who passes rubbish still gets the error naming what is actually
    wrong rather than being told the cell does not exist. Only a string that
    parses reaches the existence question.

    Every atom is judged separately, and that is the whole point rather than
    a detail. ``is_valid_cell`` answers ``bool`` for one cell and
    ``list[bool]`` for a compositional index, so the obvious guard --
    ``if not is_valid_cell(index)`` -- measures the truth of a *list* when it
    is handed several cells, and a non-empty list is always true. Written
    that way, a compositional index carrying one denied cell among valid ones
    passes the guard untouched.

    Args:
        index: A compositional index string, or an iterable of them.

    Raises:
        InvalidIndexError: If the string is not a well-formed index. Raised
            by ``decompose``, and its subclasses carry the specific fault.
        NonExistentCellError: If any cell of the index names no cell under
            the specification.
    """
    # Deferred: ``boundary`` imports ``resolutions``, and ``resolutions`` is
    # one of the modules that has to call this guard. Importing at module
    # level would close that ring.
    from .boundary import is_valid_cell
    from .index import decompose

    entries = [index] if isinstance(index, str) else list(index)
    for entry in entries:
        for atom in decompose(entry):
            if not is_valid_cell(atom):
                raise NonExistentCellError(f"index names no cell: {atom}")
