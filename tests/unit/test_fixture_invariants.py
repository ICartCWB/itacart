"""Every cell fixture states the resolution it actually has.

A fixture that declares one resolution and carries another is not an
inert typo: a test reading the declaration and asserting against the
string proves nothing, and the discrepancy survives every green run.
``sydney_cell`` carried one from F1 to F7, was reported twice, and was
corrected pointwise once without the class of error being closed.

So the check is an invariant over the whole conftest rather than an
assertion about one fixture. It reads the source with ``ast``, which
means a fixture added tomorrow is covered the day it is added, with no
list to remember to extend. This is the lesson of the public-surface
invariant: pin the property, not the instance.

The declaration lives in the docstring because that is where a reader
looks for it, so that is what the invariant reads.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest

import itacart

CONFTEST = Path(__file__).resolve().parents[1] / "conftest.py"

_DECLARED = re.compile(r"\bresolutions?\s+(\d+)(?:\s+and\s+(\d+))?", re.IGNORECASE)


def _is_fixture(node: ast.FunctionDef) -> bool:
    """Whether a function carries a ``pytest.fixture`` decorator."""
    for decorator in node.decorator_list:
        target = decorator.func if isinstance(decorator, ast.Call) else decorator
        name = (
            target.attr
            if isinstance(target, ast.Attribute)
            else getattr(target, "id", "")
        )
        if name == "fixture":
            return True
    return False


def _returned_literal(node: ast.FunctionDef) -> str | None:
    """The string a fixture returns, when it returns a literal one."""
    for statement in ast.walk(node):
        if isinstance(statement, ast.Return) and statement.value is not None:
            try:
                value = ast.literal_eval(statement.value)
            except ValueError:
                return None
            return value if isinstance(value, str) else None
    return None


def _cell_fixtures() -> list[tuple[str, str, str]]:
    """``(fixture name, docstring, index string)`` for every cell fixture."""
    tree = ast.parse(CONFTEST.read_text(encoding="utf-8"))
    out: list[tuple[str, str, str]] = []
    for node in tree.body:
        if not isinstance(node, ast.FunctionDef) or not _is_fixture(node):
            continue
        value = _returned_literal(node)
        if value is None or not itacart.is_valid_index(value):
            continue
        out.append((node.name, ast.get_docstring(node) or "", value))
    return out


def test_the_conftest_holds_cell_fixtures_to_find() -> None:
    """The invariant is worthless if it silently matches nothing."""
    names = [name for name, _, _ in _cell_fixtures()]
    assert len(names) >= 3
    assert "sydney_cell" in names
    assert "central_park_index" in names


@pytest.mark.parametrize(
    "name, docstring, index",
    _cell_fixtures(),
    ids=[name for name, _, _ in _cell_fixtures()],
)
def test_declared_resolution_matches_the_parsed_one(
    name: str, docstring: str, index: str
) -> None:
    """The docstring's resolution is the one the index actually has."""
    match = _DECLARED.search(docstring)
    assert match is not None, (
        f"fixture {name!r} returns an index but its docstring declares no "
        "resolution; every cell fixture has to state one so that this "
        "invariant can check it"
    )
    declared = {int(group) for group in match.groups() if group is not None}
    present = {itacart.get_resolution(cell) for cell in itacart.decompose(index)}
    assert declared == present, (
        f"fixture {name!r} declares resolution(s) {sorted(declared)} but "
        f"parses to {sorted(present)}"
    )


def test_every_cell_fixture_is_a_cell_that_exists() -> None:
    """A fixture naming a structurally impossible cell would pin nonsense."""
    for name, _, index in _cell_fixtures():
        for cell in itacart.decompose(index):
            assert itacart.is_valid_cell(cell), (name, cell)
