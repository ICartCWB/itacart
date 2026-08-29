"""Tests for :mod:`itacart.exceptions`.

Two things are being protected here: that a caller can guard the whole
package with one ``except ITACaRTError``, and that the three names cited
verbatim in later acceptance criteria keep both their spelling and their
place in the hierarchy.
"""

from __future__ import annotations

import inspect

import pytest

from itacart import exceptions as exc

ALL_ERRORS = [
    obj
    for _, obj in inspect.getmembers(exc, inspect.isclass)
    if issubclass(obj, BaseException) and obj.__module__ == exc.__name__
]


def test_every_error_derives_from_the_base() -> None:
    assert ALL_ERRORS, "no exception classes found"
    for error in ALL_ERRORS:
        assert issubclass(error, exc.ITACaRTError)


def test_base_derives_from_exception_not_baseexception() -> None:
    """Never a BaseException: it must not slip past a bare ``except``."""
    assert issubclass(exc.ITACaRTError, Exception)
    assert exc.ITACaRTError.__bases__ == (Exception,)


def test_all_matches_the_module_contents() -> None:
    declared = set(exc.__all__)
    defined = {error.__name__ for error in ALL_ERRORS}
    assert declared == defined


def test_every_error_has_a_docstring() -> None:
    for error in ALL_ERRORS:
        assert error.__doc__, f"{error.__name__} has no docstring"


@pytest.mark.parametrize(
    ("name", "parent"),
    [
        ("InvalidIndexError", "ITACaRTError"),
        ("InvalidQuadrantError", "InvalidIndexError"),
        ("InvalidRefinementCodeError", "InvalidIndexError"),
        ("NonAtomicIndexError", "ITACaRTError"),
        ("ResolutionError", "ITACaRTError"),
        ("MaxResolutionError", "ResolutionError"),
        ("MinResolutionError", "ResolutionError"),
        ("DomainError", "ITACaRTError"),
        ("NonExistentCellError", "DomainError"),
        ("AntemeridianError", "DomainError"),
        ("GeometryError", "ITACaRTError"),
        ("UnsupportedGeometryTypeError", "GeometryError"),
        ("DensificationError", "GeometryError"),
        ("ConvergenceError", "ITACaRTError"),
        ("SerializationError", "ITACaRTError"),
        ("MalformedBlobError", "SerializationError"),
        ("IncompatibleProfileError", "SerializationError"),
    ],
)
def test_hierarchy_shape(name: str, parent: str) -> None:
    error = getattr(exc, name)
    assert error.__bases__ == (getattr(exc, parent),)


@pytest.mark.parametrize(
    "name",
    ["ConvergenceError", "InvalidRefinementCodeError", "AntemeridianError"],
)
def test_load_bearing_names_exist(name: str) -> None:
    """Cited verbatim by F1 crit. 5, F2 crit. 5 and F7 crit. 9."""
    assert hasattr(exc, name)
    assert issubclass(getattr(exc, name), exc.ITACaRTError)


def test_antemeridian_error_is_a_domain_error() -> None:
    """F7 raises it; F4 callers catching DomainError must see it."""
    with pytest.raises(exc.DomainError):
        raise exc.AntemeridianError("crosses 180 outside an extension zone")


def test_non_existent_cell_is_a_domain_error() -> None:
    """The western X = 0 rule of F4 is a domain condition, not a syntax one."""
    assert issubclass(exc.NonExistentCellError, exc.DomainError)
    assert not issubclass(exc.NonExistentCellError, exc.InvalidIndexError)


def test_invalid_refinement_code_is_an_index_error() -> None:
    """The alphabet inversion surfaces as a syntax error, caught by
    ``except InvalidIndexError``.
    """
    with pytest.raises(exc.InvalidIndexError):
        raise exc.InvalidRefinementCodeError("C2 at an even resolution")


def test_convergence_error_is_not_a_geometry_error() -> None:
    """Vincenty failing is a numerical event, not bad input geometry."""
    assert not issubclass(exc.ConvergenceError, exc.GeometryError)
    assert exc.ConvergenceError.__bases__ == (exc.ITACaRTError,)


def test_errors_carry_their_message() -> None:
    error = exc.ResolutionError("resolution 14 outside 0..13")
    assert str(error) == "resolution 14 outside 0..13"
