"""Exception hierarchy for the itacart package.

Provenance: adapted from ``itacart_core/exceptions.py``.

Every error raised by the package derives from :class:`ITACaRTError`, so a
caller can guard an entire pipeline with one ``except`` and never catch a
bare ``ValueError`` that leaked from somewhere else::

    ITACaRTError
    |
    +-- InvalidIndexError            index string is malformed
    |     +-- InvalidQuadrantError
    |     +-- InvalidRefinementCodeError
    +-- NonAtomicIndexError          several cells where one was required
    +-- ResolutionError              resolution out of range or inapplicable
    |     +-- MaxResolutionError
    |     +-- MinResolutionError
    +-- DomainError                  position or cell outside the domain
    |     +-- NonExistentCellError
    |     +-- AntemeridianError
    +-- GeometryError                input geometry invalid or unsupported
    |     +-- UnsupportedGeometryTypeError
    |     +-- DensificationError
    +-- ConvergenceError             iterative geodesic did not converge
    +-- SerializationError           binary encoding or decoding failed
          +-- MalformedBlobError
          +-- IncompatibleProfileError

Three names are load-bearing, cited verbatim in the acceptance criteria of
later phases: :class:`ConvergenceError` (F1),
:class:`InvalidRefinementCodeError` (F2) and :class:`AntemeridianError`
(F7). Renaming any of them breaks those phases by name, not by logic.
"""

from __future__ import annotations


class ITACaRTError(Exception):
    """Base class for every error raised by this package."""


# -- Index / addressing -----------------------------------------------------


class InvalidIndexError(ITACaRTError):
    """The compositional index string is syntactically malformed."""


class InvalidQuadrantError(InvalidIndexError):
    """Quadrant code is not one of NE, NW, SE, SW."""


class InvalidRefinementCodeError(InvalidIndexError):
    """A refinement code does not belong to the alphabet of its resolution.

    Even resolutions accept ``1``-``4``; odd resolutions accept ``A1``-``E5``.
    """


class NonAtomicIndexError(ITACaRTError):
    """Operation requires a single terminal cell but the index holds several.

    Raised only by the few operations with no meaningful vectorised
    semantics. Most of the public API accepts compositional indices and
    returns positionally aligned results instead.
    """


# -- Resolution -------------------------------------------------------------


class ResolutionError(ITACaRTError):
    """Resolution is outside the valid range or invalid for the operation."""


class MaxResolutionError(ResolutionError):
    """Attempt to refine beyond resolution 13."""


class MinResolutionError(ResolutionError):
    """Attempt to ascend above resolution 0."""


# -- Domain / boundary ------------------------------------------------------


class DomainError(ITACaRTError):
    """Coordinates or cell fall outside the addressable ITACaRT domain."""


class NonExistentCellError(DomainError):
    """Cell is structurally non-existent under the specification.

    Chiefly the resolution-1 cells with X index equal to 0 in the western
    quadrants, which are absorbed by the triangular prime-meridian cells.
    """


class AntemeridianError(DomainError):
    """Geometry crosses the antemeridian outside a defined extension zone."""


# -- Geometry ---------------------------------------------------------------


class GeometryError(ITACaRTError):
    """Input geometry is invalid or unsupported."""


class UnsupportedGeometryTypeError(GeometryError):
    """Geometry type is not one of the supported OGC SFA types."""


class DensificationError(GeometryError):
    """Densification could not be carried out for the given segment."""


# -- Geodesy ----------------------------------------------------------------


class ConvergenceError(ITACaRTError):
    """An iterative geodesic computation failed to converge (Vincenty).

    Raised chiefly for nearly antipodal pairs, where the inverse solution
    oscillates. Raising is deliberate: returning an unconverged value
    would look like a distance and silently poison whatever consumes it.
    """


# -- Serialization ----------------------------------------------------------


class SerializationError(ITACaRTError):
    """Base class for binary encoding and decoding failures."""


class MalformedBlobError(SerializationError):
    """Binary blob failed structural validation."""


class IncompatibleProfileError(SerializationError):
    """Blob declares a profile combination this build cannot honour."""


__all__ = [
    "ITACaRTError",
    "InvalidIndexError",
    "InvalidQuadrantError",
    "InvalidRefinementCodeError",
    "NonAtomicIndexError",
    "ResolutionError",
    "MaxResolutionError",
    "MinResolutionError",
    "DomainError",
    "NonExistentCellError",
    "AntemeridianError",
    "GeometryError",
    "UnsupportedGeometryTypeError",
    "DensificationError",
    "ConvergenceError",
    "SerializationError",
    "MalformedBlobError",
    "IncompatibleProfileError",
]
