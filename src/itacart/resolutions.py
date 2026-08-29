"""Resolution table helpers.

Satisfies OGC DGGS Core requirements 14 and 15 (hierarchical grid sequence).

Provenance: ``itacart_core/resolutions.py`` and Table 1 of the paper.

Correspondence with the origin, which stores sides as exact integer
centimetres and derives everything from them:

===============================  ==============================
itacart_core/resolutions.py      this module
===============================  ==============================
``base_length_m``                :func:`cell_size`
``cell_area_m2``                 :func:`nominal_cell_area`
``children_per_cell``            :func:`refinement_ratio`
``subdivisions_per_axis``        :func:`linear_refinement_ratio`
``is_tokenizable_resolution``    :func:`is_tokenizable_resolution`
``is_even_resolution``           -- (see ``constants.refinement_alphabet``)
===============================  ==============================

Values agree with the origin exactly, to zero ulps, on all 13 metric rows
for side, area and both ratios. The one divergence is resolution 1 under
:func:`is_tokenizable_resolution`, and it is F0's ``D-0.12``, not this
phase's.

The module is a typed, validating facade over the tuples in
:mod:`itacart.constants`. It restates no value from Table 1: every figure
returned here is read from a constant that F0 audited against the paper,
so a wrong row can only be wrong in one place.

Two distinctions in this module are load-bearing and easy to get wrong.

**Area ratio against linear ratio.** :func:`refinement_ratio` reports how
many children a refinement produces -- 4 or 25. The *side* of the child is
half or a fifth of the parent's, never a quarter or a twenty-fifth. The
descent in :mod:`itacart.cells` needs the linear factor and child
enumeration needs the area factor, and substituting one for the other
yields a wrong cell with no exception raised.
:func:`linear_refinement_ratio` exists so neither has to be spelled
inline (``D-3.1``).

**Nominal area against effective area.** :func:`nominal_cell_area` takes a
resolution and :func:`effective_cell_area` takes a cell. The differing
signatures are deliberate (``D-0.5``): a caller cannot reach for the
cheaper one without noticing that it never saw the cell.
"""

from __future__ import annotations

import math
from typing import Literal

from .constants import (
    ANALYSIS_SCALE,
    CELL_AREA_M2,
    CELL_SIZE_M,
    MAX_RESOLUTION,
    MIN_RESOLUTION,
    QUADRANTS,
    REFINEMENT_ALPHABET,
    REFINEMENT_RATIO,
    RESOLUTION_COUNT,
    TOKENIZABLE_RESOLUTIONS,
    VISUALIZATION_SCALE,
)
from .exceptions import ResolutionError
from .index import decompose, is_atomic, iter_cells, split_components

__all__ = [
    "get_resolution",
    "refinement_ratio",
    "linear_refinement_ratio",
    "cell_size",
    "nominal_cell_area",
    "effective_cell_area",
    "scale_for_resolution",
    "resolution_for_scale",
    "resolution_table",
    "is_tokenizable_resolution",
]

ScaleKind = Literal["visualization", "analysis"]

_SCALE_COLUMN: dict[str, tuple[int | None, ...]] = {
    "visualization": VISUALIZATION_SCALE,
    "analysis": ANALYSIS_SCALE,
}


# --------------------------------------------------------------------------
# Validation helpers
# --------------------------------------------------------------------------


def _check_resolution(resolution: int, minimum: int, reason: str) -> None:
    """Reject a resolution outside ``minimum..MAX_RESOLUTION``.

    ``bool`` is excluded explicitly: ``True`` is an ``int`` and would
    silently address resolution 1, which is a real level and so would not
    fail anywhere downstream.
    """
    if not isinstance(resolution, int) or isinstance(resolution, bool):
        raise ResolutionError(
            f"resolution must be an int, got {type(resolution).__name__}"
        )
    if not MIN_RESOLUTION <= resolution <= MAX_RESOLUTION:
        raise ResolutionError(
            f"resolution {resolution} outside {MIN_RESOLUTION}..{MAX_RESOLUTION}"
        )
    if resolution < minimum:
        raise ResolutionError(f"resolution {resolution} {reason}")


def _check_scale_kind(kind: str) -> tuple[int | None, ...]:
    """Resolve a scale family name to its column, rejecting unknown names."""
    try:
        return _SCALE_COLUMN[kind]
    except KeyError:
        raise ResolutionError(
            f"unknown scale kind {kind!r}; expected 'visualization' or 'analysis'"
        ) from None


def _scale_at(column: tuple[int | None, ...], resolution: int) -> int:
    """Read a scale column at a resolution known to carry a value."""
    denominator = column[resolution]
    assert denominator is not None  # resolutions 1..13 all carry both scales
    return denominator


# --------------------------------------------------------------------------
# Resolution of a cell
# --------------------------------------------------------------------------


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
    levels = {len(split_components(atom)) - 1 for atom in decompose(cell)}
    if len(levels) > 1:
        raise ResolutionError(
            "index is mixed-resolution, holding terminal cells at levels "
            f"{sorted(levels)}; decompose() it and ask per cell"
        )
    return levels.pop()


# --------------------------------------------------------------------------
# Refinement
# --------------------------------------------------------------------------


def refinement_ratio(resolution: int) -> int:
    """Number of children produced when refining *into* ``resolution``.

    Even resolutions come from a 1-to-4 subdivision, odd resolutions from
    a 1-to-25 subdivision. Resolution 1 is a special case: it is the base
    10 km grid addressed by Cartesian integers, not a refinement.

    This is an **area** ratio. The side of the child shrinks by
    :func:`linear_refinement_ratio`, which is its square root.

    Args:
        resolution: Target resolution level, 2 to 13.

    Returns:
        ``4`` for even resolutions, ``25`` for odd ones.

    Raises:
        ResolutionError: If ``resolution`` is below 2 or above 13.
    """
    _check_resolution(
        resolution,
        2,
        "is not produced by a refinement: 0 is a quadrant and 1 is the "
        "Cartesian base grid",
    )
    ratio = REFINEMENT_RATIO[resolution]
    assert ratio is not None  # guaranteed by the bound above
    return ratio


def linear_refinement_ratio(resolution: int) -> int:
    """Factor by which the cell *side* shrinks when refining into ``resolution``.

    ``2`` where :func:`refinement_ratio` is 4, ``5`` where it is 25 --
    the exact integer square root, since the cell is a parallelogram whose
    base and height are equal and both divide by the same factor.

    Coordinate descent uses this one; child enumeration uses
    :func:`refinement_ratio`. Table 1 states both, in different columns,
    and confusing them produces a wrong cell in silence rather than an
    exception, which is the most expensive failure mode in this package
    (``D-3.1``).

    Called ``subdivisions_per_axis`` in ``itacart_core/resolutions.py``.
    That version returns ``1`` at resolution 1, treating the base grid as
    an anchor; this one raises instead (``D-3.4``). A descent loop that
    reads ``1`` for resolution 1 divides the side by one and silently
    produces a cell at the wrong level, which is the same failure mode the
    area/linear split exists to prevent.

    Args:
        resolution: Target resolution level, 2 to 13.

    Returns:
        ``2`` for even resolutions, ``5`` for odd ones.

    Raises:
        ResolutionError: If ``resolution`` is below 2 or above 13.
    """
    return math.isqrt(refinement_ratio(resolution))


# --------------------------------------------------------------------------
# Size and area
# --------------------------------------------------------------------------


def cell_size(resolution: int) -> float:
    """Base and height of a cell, in metres.

    Args:
        resolution: Resolution level, 1 to 13.

    Returns:
        Edge length in metres, from 10 000 m down to 0.01 m.

    Raises:
        ResolutionError: If ``resolution`` is 0 or out of range.
    """
    _check_resolution(resolution, 1, "is a global quadrant and has no metric size")
    size = CELL_SIZE_M[resolution]
    assert size is not None  # guaranteed by the bound above
    return size


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
    _check_resolution(resolution, 1, "is a global quadrant and has no nominal area")
    area = CELL_AREA_M2[resolution]
    assert area is not None  # guaranteed by the bound above
    return area


def effective_cell_area(cell: str) -> float | list[float]:
    """True area of a specific cell, accounting for boundary clipping.

    Equals :func:`nominal_cell_area` for parallelogram and triangular
    cells. For trapezoidal cells the longer base is clipped at the
    boundary line, so the area is computed from the actual vertices.

    Cadastral use makes this distinction legally significant: returning
    the nominal area for a clipped cell would misstate a parcel.

    An unclipped cell is answered from the resolution table rather than
    from its own vertices. The plane is equal-area by construction, so
    the two agree to the last bit the projection can carry, and reading
    the table keeps the common case free of a shoelace sum over
    twenty-million-metre coordinates. A clipped cell is measured.

    Args:
        cell: Compositional index string.

    Returns:
        Area in square metres for a single terminal cell, or a
        positionally aligned list when the index holds several.
    """
    from .boundary import absorbs_border, plane_ring, ring_area

    values = []
    for atom in iter_cells(cell):
        if not absorbs_border(atom):
            values.append(nominal_cell_area(get_resolution(atom)))
        else:
            values.append(ring_area(plane_ring(atom)[1]))
    return values[0] if is_atomic(cell) else values


# --------------------------------------------------------------------------
# Cartographic scale
# --------------------------------------------------------------------------


def scale_for_resolution(resolution: int, kind: ScaleKind = "visualization") -> int:
    """Cartographic scale denominator matching a resolution.

    ``visualization`` follows the minimum visible line of 0.1 mm on paper
    maps (Jenny et al., 2008), hence ``denominator = side_m / 1e-4``.
    ``analysis`` follows sampling theory at 0.5 mm (Tobler, 1987), hence
    ``denominator = side_m / 5e-4``. The ratio between the two columns is
    therefore 5 at every resolution.

    Args:
        resolution: Resolution level, 1 to 13.
        kind: Which of the two scale families to report.

    Returns:
        Scale denominator, e.g. ``10_000`` for 1:10 000.

    Raises:
        ResolutionError: If ``resolution`` is 0 or out of range, or if
            ``kind`` is not a known scale family.
    """
    column = _check_scale_kind(kind)
    _check_resolution(resolution, 1, "is a global quadrant and has no scale")
    return _scale_at(column, resolution)


def resolution_for_scale(denominator: int, kind: ScaleKind = "visualization") -> int:
    """Resolution adequate for a target cartographic scale.

    Inverse of :func:`scale_for_resolution`; the practical entry point for
    a surveyor who knows the map scale but not the grid level. On the 13
    denominators printed in Table 1 it returns exactly the row they came
    from, in both families.

    Between table rows the two families bound fineness from opposite
    sides, because the two rules they encode are opposite (``D-3.3``):

    - ``visualization`` asks that a cell still be **drawable**, so the
      answer is the finest level whose cells stay at or above 0.1 mm on
      paper. At 1:2 000 that is resolution 10; resolution 11 would render
      at 0.05 mm and vanish.
    - ``analysis`` asks that a cell be **fine enough to sample**, so the
      answer is the coarsest level at or below the 0.5 mm ground
      distance. At 1:3 000 that is resolution 9.

    Args:
        denominator: Scale denominator, e.g. ``1000`` for 1:1 000.
        kind: Which of the two scale families to interpret against.

    Returns:
        Resolution level, 1 to 13.

    Raises:
        ResolutionError: If ``kind`` is not a known scale family, if
            ``denominator`` is not a positive integer, or if no
            resolution satisfies the target: a visualization scale
            coarser than 1:100 000 000 has no drawable level, and an
            analysis scale finer than 1:20 has no level fine enough.
    """
    column = _check_scale_kind(kind)
    if not isinstance(denominator, int) or isinstance(denominator, bool):
        raise ResolutionError(
            f"denominator must be an int, got {type(denominator).__name__}"
        )
    if denominator <= 0:
        raise ResolutionError(f"denominator must be positive, got {denominator}")

    levels = range(1, RESOLUTION_COUNT)
    if kind == "visualization":
        drawable = [r for r in levels if _scale_at(column, r) >= denominator]
        if not drawable:
            raise ResolutionError(
                f"no resolution is drawable at 1:{denominator}; the coarsest "
                f"level stays visible only up to 1:{_scale_at(column, 1)}"
            )
        return max(drawable)

    fine_enough = [r for r in levels if _scale_at(column, r) <= denominator]
    if not fine_enough:
        raise ResolutionError(
            f"no resolution is fine enough for 1:{denominator}; the finest "
            f"level samples 1:{_scale_at(column, MAX_RESOLUTION)}"
        )
    return min(fine_enough)


# --------------------------------------------------------------------------
# The table itself
# --------------------------------------------------------------------------


def resolution_table() -> list[dict[str, object]]:
    """Full Table 1 of the paper as structured records.

    Each record carries ``resolution``, ``cell_size_m``, ``cell_area_m2``,
    ``refinement``, ``index_alphabet``, ``visualization_scale`` and
    ``analysis_scale``. Supports the OGC ``describe`` operation and gives
    documentation a single source of truth.

    Resolution 0 carries ``None`` in every metric column and the four
    quadrant codes as its alphabet, matching the dashes printed in the
    published table. Resolution 1 carries metric values but no refinement
    ratio and no alphabet, being addressed by an ``XXXX/YYYY`` pair.

    Returns:
        One record per resolution, ordered 0 to 13. Freshly built on each
        call, so a caller may mutate the result without corrupting the
        constants it was read from.
    """
    return [
        {
            "resolution": resolution,
            "cell_size_m": CELL_SIZE_M[resolution],
            "cell_area_m2": CELL_AREA_M2[resolution],
            "refinement": REFINEMENT_RATIO[resolution],
            "index_alphabet": (
                QUADRANTS
                if resolution == MIN_RESOLUTION
                else REFINEMENT_ALPHABET[resolution]
            ),
            "visualization_scale": VISUALIZATION_SCALE[resolution],
            "analysis_scale": ANALYSIS_SCALE[resolution],
        }
        for resolution in range(RESOLUTION_COUNT)
    ]


def is_tokenizable_resolution(resolution: int) -> bool:
    """Whether cells at this resolution carry a whole-number metric area.

    Decimal convergence is what makes one token equal one standard unit
    of area in a blockchain registry. Resolution 13 gives exactly 1 cm2.

    True at resolution 1 and at every odd resolution, whose sides are
    exact powers of ten metres. Even resolutions carry sides of the form
    5 x 10^k, regular but not decimal.

    Resolution 1 is a deliberate divergence from the origin, which
    enumerates only the odd levels. Its side of 10 km is 10^4 m and its
    area 10^8 m2, so it satisfies the stated property; the paper does not
    enumerate the set either way. F0 took that reading as ``D-0.12`` and
    the choice remains F0's to revisit, not this phase's.

    Provenance: ``itacart_core/engine.py`` (``IDGGSEngine``) and
    ``itacart_core/resolutions.py`` (``_TOKENIZABLE_RESOLUTIONS``).

    Args:
        resolution: Resolution level, 1 to 13.

    Returns:
        ``True`` when the nominal area is a whole number in its natural
        metric unit.

    Raises:
        ResolutionError: If ``resolution`` is 0 or out of range.
    """
    _check_resolution(resolution, 1, "is a global quadrant and has no area to tokenize")
    return resolution in TOKENIZABLE_RESOLUTIONS
