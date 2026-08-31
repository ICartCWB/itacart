"""Vector geometry against the grid: filling, vertex mapping, canonical form.

Two ways to represent a vector feature, both described in section 4 of the
paper:

**Cell filling** - the feature is the set of cells it covers. Area follows
from a cell count with no projection distortion, but the index grows
verbose at fine resolutions.

**Vertex representation** - only the cells holding the polygon vertices are
kept, as in conventional vector data. Far more compact, and the form the
binary encodings in :mod:`itacart.serialization` are built on.

The descent
-----------

Filling never tests cells one by one. It projects the geometry onto the
parallels plane, shears the 45-degree tiling onto a square lattice exactly
as :mod:`itacart.cells` does, and walks the hierarchy: a node disjoint from
the geometry is dropped whole, a node wholly inside is accepted whole, and
only a node straddling the outline is subdivided. Cost therefore tracks the
outline, not the area.

The shear is the same one the quantizer uses, ``u = |x| + |y|`` and
``v = |y|``, applied per quadrant as a signed affine map so that the
mirrored quadrants keep their own handedness. Its Jacobian is 1, so a cell
count on the lattice is a cell count on the ellipsoid.

Scope
-----

Filling is defined on the parallelogram interior. Geometry reaching the
prime-meridian column, the last lattice column of its row, or the polar row
is refused rather than filled, because in those three families the cell is
not the sheared square the descent tests and the containment chain of
:func:`polyfill` would no longer be sound. See :func:`polyfill` for the
exact predicate.

Provenance: ``itacart_core/cell_filling.py``, ``densification.py``,
``geometry_blob.py`` (``canonicalize_rings``) and
``cadastral_processor/vertex_extractor.py``. Two of the seven functions
have no counterpart there: ``cells_to_geometry`` is absent from all four
files, and the origin's filling implements one containment mode of the
three.
"""

from __future__ import annotations

import math
import threading
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from typing import TYPE_CHECKING, Callable, Literal, Sequence, TypeVar

from .boundary import crosses_antemeridian, last_lattice_column
from .cells import _ROW_LETTERS, cell_to_anchor, geo_to_cell
from .constants import (
    DESCENT_CLOSE,
    DESCENT_OPEN,
    MAX_RESOLUTION,
    RES1_DIGITS,
    RES1_SEPARATOR,
    SIBLING_SEPARATOR,
)
from .exceptions import (
    AntemeridianError,
    DensificationError,
    DomainError,
    GeometryError,
    NonExistentCellError,
    ResolutionError,
    UnsupportedGeometryTypeError,
)
from .geodesy import direct_geodesic, geodetic_to_sinusoidal, inverse_geodesic
from .resolutions import cell_size, linear_refinement_ratio, refinement_ratio

if TYPE_CHECKING:
    from shapely.geometry import LineString, Point, Polygon
    from shapely.geometry.base import BaseGeometry
    from shapely.prepared import PreparedGeometry

__all__ = [
    "polyfill",
    "count_internal_cells",
    "vertex_to_cell",
    "cells_to_geometry",
    "densify_orthodromic",
    "densify_segment",
    "canonicalize_rings",
]

Containment = Literal["center", "intersects", "contains"]

_L1: float = cell_size(1)
"""Resolution-1 cell side, in metres."""

MAX_FILL_CELLS: int = 200_000_000
"""Ceiling on cells emitted by one :func:`polyfill` call.

A fill is bounded by the geometry's area divided by the cell area, which
at resolution 13 is one square centimetre: a single football pitch holds
seventy billion cells. Without a ceiling the call does not fail, it
simply never returns, and the caller cannot tell the two apart.

The number is measured, not chosen. The index a fill returns costs about
three bytes per cell, so two hundred million cells is a string of some
six hundred megabytes and a peak of roughly twice that while it is being
joined. That is large but survivable on an ordinary machine, and it is
the point at which the output itself, rather than any bookkeeping around
it, becomes the thing that does not fit.

Raising it further is a memory question, not a policy one: measure
``len(polyfill(...))`` against available memory before doing so.
"""

# --------------------------------------------------------------------------
# Plane and lattice helpers
# --------------------------------------------------------------------------


_QUADRANT_SHEAR: dict[str, tuple[float, float, float, float]] = {
    "NE": (1.0, 1.0, 0.0, 1.0),
    "NW": (-1.0, 1.0, 0.0, 1.0),
    "SE": (1.0, -1.0, 0.0, -1.0),
    "SW": (-1.0, -1.0, 0.0, -1.0),
}
"""Per-quadrant affine shear ``(a, b, d, e)`` taking plane ``(x, y)`` to ``(u, v)``.

``u = a*x + b*y`` and ``v = d*x + e*y``. Each is the composition of the
quadrant mirror with the unit shear of :mod:`itacart.cells`, so that
``u = |x| + |y|`` and ``v = |y|`` hold inside the quadrant. The determinant
is 1 in absolute value for all four, which is why a count on the lattice is
a count on the ellipsoid.

Written per quadrant rather than as ``abs()`` because ``abs()`` folds the
four quadrants onto one and loses which one a piece came from. The origin
folds, and recovers the quadrant afterwards; recovering it afterwards is
what makes a cross-quadrant geometry impossible to fill.
"""

_PLANE_LIMIT: float = 4.0e7
"""Half-width of the clipping box, in metres: past any plane coordinate."""


def _check_resolution(resolution: int) -> None:
    """Reject a resolution outside 1..13."""
    if not isinstance(resolution, int) or isinstance(resolution, bool):
        raise ResolutionError(
            f"resolution must be an int, got {type(resolution).__name__}"
        )
    if not 1 <= resolution <= MAX_RESOLUTION:
        raise ResolutionError(f"resolution {resolution} outside 1..{MAX_RESOLUTION}")


def _check_jobs(n_jobs: int) -> None:
    """Reject a worker count below one."""
    if not isinstance(n_jobs, int) or isinstance(n_jobs, bool):
        raise ValueError(f"n_jobs must be an int, got {type(n_jobs).__name__}")
    if n_jobs < 1:
        raise ValueError(f"n_jobs must be >= 1, got {n_jobs}")


def _leaves_between(level: int, target: int) -> int:
    """How many ``target``-resolution cells sit under one ``level`` cell."""
    total = 1
    for step in range(level + 1, target + 1):
        total *= refinement_ratio(step)
    return total


def _child_code(row: int, column: int, level: int) -> str:
    """Refinement code of the child at ``(row, column)`` of a ``level`` grid.

    Rows run south to north and columns west to east *in the cell's own
    quadrant*, which is the frame the descent works in. The quaternary
    code is ``row * 2 + column + 1``; the quinary code is the row letter
    followed by the one-based column. Both are the inverse of the mapping
    :func:`itacart.cells._anchor_on_plane` applies when it ascends, and
    the letters come from that module so the north-south orientation
    keeps a single definition in the package.
    """
    if linear_refinement_ratio(level) == 2:
        return str(row * 2 + column + 1)
    return f"{_ROW_LETTERS[row]}{column + 1}"


def _project(geometry: "BaseGeometry") -> "BaseGeometry":
    """Map a geodetic geometry onto the parallels plane, vertex by vertex.

    ``shapely.ops.transform`` is not used: the projection is applied to
    the coordinates that are there, and no vertex is invented. A caller
    wanting the edges to follow geodesics densifies first, which is what
    :func:`polyfill` does on its behalf.
    """
    from shapely.geometry import LineString, MultiPolygon, Point, Polygon

    if isinstance(geometry, Point):
        return Point(*geodetic_to_sinusoidal(geometry.x, geometry.y))
    if isinstance(geometry, LineString):
        return LineString([geodetic_to_sinusoidal(x, y) for x, y in geometry.coords])
    if isinstance(geometry, Polygon):
        return Polygon(
            [geodetic_to_sinusoidal(x, y) for x, y in geometry.exterior.coords],
            [
                [geodetic_to_sinusoidal(x, y) for x, y in hole.coords]
                for hole in geometry.interiors
            ],
        )
    if isinstance(geometry, MultiPolygon):
        return MultiPolygon([_project(part) for part in geometry.geoms])
    raise UnsupportedGeometryTypeError(
        f"cannot fill a {geometry.geom_type}; supported types are Point, "
        "LineString, Polygon and MultiPolygon"
    )


def _quadrant_pieces(plane: "BaseGeometry") -> list[tuple[str, "BaseGeometry"]]:
    """Split a plane geometry into its four quadrant parts, dropping empties.

    A geometry lying in one quadrant yields one piece, which is the
    ordinary case. A geometry straddling an axis yields two or four, each
    filled independently: the half-open convention already awards every
    axis position to exactly one cell, so the parts do not double count.
    """
    from shapely.geometry import box

    out: list[tuple[str, "BaseGeometry"]] = []
    limit = _PLANE_LIMIT
    windows = {
        "NE": box(0.0, 0.0, limit, limit),
        "NW": box(-limit, 0.0, 0.0, limit),
        "SE": box(0.0, -limit, limit, 0.0),
        "SW": box(-limit, -limit, 0.0, 0.0),
    }
    for quadrant, window in windows.items():
        piece = plane.intersection(window)
        if piece.is_empty:
            continue
        if piece.geom_type in {"Point", "LineString", "MultiPoint", "MultiLineString"}:
            if plane.geom_type in {"Polygon", "MultiPolygon"}:
                continue  # a polygon touching the axis, not crossing it
        out.append((quadrant, piece))
    return out


def _to_lattice(piece: "BaseGeometry", quadrant: str) -> "BaseGeometry":
    """Shear a quadrant piece onto the square lattice."""
    from shapely.affinity import affine_transform

    a, b, d, e = _QUADRANT_SHEAR[quadrant]
    return affine_transform(piece, [a, b, d, e, 0.0, 0.0])


def _check_addressable(quadrant: str, column: int, row: int) -> None:
    """Refuse a resolution-1 cell that is not an ordinary parallelogram.

    Three families are refused, and the reason is the same in all three:
    the descent tests a sheared square, and in these families the cell is
    not that square.

    * ``column == 0`` is the prime-meridian column, whose cells are
      triangles spanning both sides of the line.
    * the last lattice column of a row absorbs the strip between itself
      and the domain border, so it is wider than the square.
    * the polar row is clipped by the pole and carries a fraction of the
      nominal area.

    The screen is applied at resolution 1 and is therefore conservative:
    a geometry inside base column 500 cannot reach a border-absorbing
    cell at any finer resolution, so refusing the whole base column
    refuses more than strictly necessary and never less.

    **Order matters.** The polar row is tested before the last column,
    because the polar row *is* a row whose only column is its last one:
    at resolution 1 it holds exactly one, so a last-column test placed
    first would answer every polar position with the wrong diagnosis and
    the polar branch would never be reached at all. Measured: the branch
    was unreachable until the two were swapped.
    """
    if column <= 0:
        raise NonExistentCellError(
            f"row {row} column {column} of quadrant {quadrant} is the "
            "prime-meridian column, whose cells are triangles; filling is "
            "defined on the parallelogram interior"
        )
    last = last_lattice_column(quadrant, row, _L1)
    if last <= 0:
        raise DomainError(
            f"row {row} of quadrant {quadrant} addresses no cell: its "
            "parallel circle is shorter than one cell side"
        )
    if last_lattice_column(quadrant, row + 1, _L1) <= 0:
        raise DomainError(
            f"row {row} of quadrant {quadrant} is the polar row, which is "
            "clipped by the pole and does not carry the nominal cell area"
        )
    if column >= last:
        raise NonExistentCellError(
            f"row {row} column {column} of quadrant {quadrant} is the last "
            "lattice column of its row, which absorbs the border strip and "
            "is not a parallelogram"
        )


def _anomalous_band(
    quadrant: str, column: int, row: int, side: float
) -> "BaseGeometry":
    """The part of a resolution-1 cell that is still anomalous at ``side``.

    The three families are bounded by lines of constant ``w = u - v`` or
    of constant ``v``, so each one is a band and their union is what the
    screen has to protect. Returned in lattice coordinates.

    Derived, and validated by enumeration against
    :func:`_check_addressable` over 88,880 nodes in four quadrants. The
    lattice column of a node of side ``s`` at ``(u0, v0)`` is
    ``(u0 - v0) / s``, which reduces to ``u_index - row`` at resolution
    1. From it, ``column`` of a child is ``column * d + (c - r)`` with
    ``c`` and ``r`` in ``[0, d)``, so a node with a positive column never
    fathers one at column zero: the meridian family does not spread, and
    neither do the other two.

    Each band is drawn one cell side wide, which is an outer bound on a
    family that is really a staircase of single cells. Refusing the
    bound refuses at most one extra cell of width and never one less,
    and one cell of width is the whole point of the exercise.
    """
    from shapely.geometry import Polygon, box
    from shapely.ops import unary_union

    u0 = (column + row) * _L1
    v0 = row * _L1
    cell = box(u0, v0, u0 + _L1, v0 + _L1)
    reach = 4.0 * _L1
    bands: list["BaseGeometry"] = []

    def by_column(lo: float, hi: float) -> "BaseGeometry":
        """The band ``lo <= u - v <= hi``, as a sheared quadrilateral."""
        return Polygon(
            [
                (v0 - reach + lo, v0 - reach),
                (v0 - reach + hi, v0 - reach),
                (v0 + reach + hi, v0 + reach),
                (v0 + reach + lo, v0 + reach),
            ]
        )

    if column <= 0:
        bands.append(by_column(-reach, side))

    top_row = int(math.floor((v0 + _L1 - side / 2.0) / side))
    last_at_top = last_lattice_column(quadrant, top_row, side)
    if last_at_top > 0:
        bands.append(by_column(last_at_top * side, reach))
    else:
        bands.append(cell)

    polar_v = None
    probe = int(math.floor(v0 / side))
    if last_lattice_column(quadrant, probe + 1, side) <= 0:
        polar_v = v0
    elif last_lattice_column(quadrant, top_row + 1, side) <= 0:
        lo, hi = probe, top_row
        while lo < hi:
            mid = (lo + hi) // 2
            if last_lattice_column(quadrant, mid + 1, side) <= 0:
                hi = mid
            else:
                lo = mid + 1
        polar_v = lo * side
    if polar_v is not None:
        bands.append(box(u0 - reach, polar_v, u0 + reach, v0 + _L1))

    return unary_union([band.intersection(cell) for band in bands])


def _base_cells(
    prepared: "PreparedGeometry",
    lattice: "BaseGeometry",
    quadrant: str,
    resolution: int,
) -> list[tuple[int, int]]:
    """Resolution-1 ``(column, row)`` pairs whose square meets ``lattice``.

    Candidates come from the bounding box, but only those a square
    actually meets are screened. Screening the whole box instead would
    refuse a geometry for a neighbouring column it never touches.

    The screen runs at ``resolution``, not at resolution 1. A base cell
    of one of the three families holds mostly ordinary descendants --
    measured, 90 of 100 for the meridian, 107 of 114 for the last
    column, 139 of 144 for the polar row -- and screening the family at
    resolution 1 refused all of them along with the few that deserved
    it. That is a band up to ten kilometres wide standing in for an
    anomaly one cell wide: a thousand times too much at resolution 7 and
    a million times at resolution 13.

    So an anomalous base cell is not refused for being one. Its
    anomalous band at the target side is computed and the geometry is
    refused only if it reaches it.
    """
    from shapely.geometry import box

    min_u, min_v, max_u, max_v = lattice.bounds
    row_lo = int(math.floor(min_v / _L1))
    row_hi = int(math.floor(max_v / _L1))
    u_lo = int(math.floor(min_u / _L1))
    u_hi = int(math.floor(max_u / _L1))
    out: list[tuple[int, int]] = []
    for row in range(row_lo, row_hi + 1):
        for u_index in range(u_lo, u_hi + 1):
            u0 = u_index * _L1
            v0 = row * _L1
            if not prepared.intersects(box(u0, v0, u0 + _L1, v0 + _L1)):
                continue
            column = u_index - row
            try:
                _check_addressable(quadrant, column, row)
            except (DomainError, NonExistentCellError):
                band = _anomalous_band(quadrant, column, row, cell_size(resolution))
                if band.is_empty or not prepared.intersects(band):
                    out.append((column, row))
                    continue
                raise
            out.append((column, row))
    return out


# --------------------------------------------------------------------------
# Cell filling
# --------------------------------------------------------------------------


def _auto_segment(resolution: int) -> float:
    """Densification threshold implied by a target resolution, in metres.

    A straight line on the parallels plane is not a geodesic, and the two
    part company by a sagitta that grows with the square of the span. A
    chord of length ``d`` departs from its geodesic by roughly
    ``d^2 / (8 R)``; requiring that to stay under half a cell side gives
    ``d = sqrt(4 R l)``. At resolution 7 that is about sixteen kilometres,
    at resolution 13 about five hundred metres.

    The result is capped at one kilometre, the blanket threshold the
    briefing already fixed, so this rule is never looser than the
    project's standing decision and is stricter wherever the cell asks
    for it.
    """
    from .constants import WGS84_A

    return min(math.sqrt(4.0 * WGS84_A * cell_size(resolution)), 1000.0)


def _accept_leaf(
    prepared: "PreparedGeometry",
    containment: Containment,
    u0: float,
    v0: float,
    side: float,
) -> bool:
    """Whether a target-resolution square is kept, given how it was reached.

    Called only for a square that the descent has found to intersect the
    geometry and *not* to be wholly inside it. Both facts are already
    established, so two of the three modes answer without a further
    predicate call, and the containment chain follows from the geometry
    of a square rather than from three independent tests:

    ``contains`` reached here means not wholly inside, so the cell is
    dropped. ``intersects`` reached here means touching, so it is kept.
    ``center`` asks the one question still open. Since the centre is a
    point of the square, wholly-inside implies centre-inside implies
    touching, and therefore ``contains`` is a subset of ``center`` is a
    subset of ``intersects`` by construction, for every square the
    descent visits.
    """
    from shapely.geometry import Point

    if containment == "contains":
        return False
    if containment == "intersects":
        return True
    return bool(prepared.contains(Point(u0 + side / 2.0, v0 + side / 2.0)))


@lru_cache(maxsize=None)
def _expansion(level: int, target: int) -> str:
    """Index suffix that expands a wholly-accepted node down to ``target``.

    Every node accepted whole at the same level expands to the same
    text, so it is built once and the callers share one string object.
    That is what keeps the interior of a fill cheap: a region wholly
    inside the geometry contributes one pointer per node rather than one
    string per cell.

    Empty when the node is already at the target resolution.
    """
    if level == target:
        return ""
    step = level + 1
    divisor = linear_refinement_ratio(step)
    inner = _expansion(step, target)
    parts = [
        _child_code(row, column, step) + inner
        for row in range(divisor)
        for column in range(divisor)
    ]
    return f"{DESCENT_OPEN}{SIBLING_SEPARATOR.join(parts)}{DESCENT_CLOSE}"


class _Budget:
    """Running cell count for one fill, checked before any expansion.

    The check has to come before :func:`_expansion` rather than after
    the fact, because expanding a resolution-1 node to resolution 13
    builds a string of a trillion cells and there is no counting that
    afterwards.
    """

    __slots__ = ("spent",)

    def __init__(self) -> None:
        self.spent = 0

    def charge(self, cells: int) -> None:
        """Add ``cells`` to the tally and refuse to pass the ceiling."""
        self.spent += cells
        if self.spent > MAX_FILL_CELLS:
            raise GeometryError(
                f"fill exceeded {MAX_FILL_CELLS} cells; pass compact=True, "
                "or fill a smaller geometry. A coarser resolution is a last "
                "resort: in cadastral use the resolution is prescribed by the "
                "mapping scale and is not the caller's to lower"
            )


def _fill_node(
    prepared: "PreparedGeometry",
    u0: float,
    v0: float,
    side: float,
    level: int,
    target: int,
    containment: Containment,
    code: str,
    compact: bool,
    budget: _Budget,
) -> str | None:
    """Walk one node and return its index fragment, or ``None`` if empty.

    The fragment is the node's own code followed by its accepted
    descendants in brackets, which is the compositional index of the
    subtree. Building it during the descent rather than collecting
    atomic indices and composing them afterwards is not a
    micro-optimisation: the atomic form of a fill costs about two
    kilobytes of live memory per cell, against three bytes per cell for
    the composed form it is about to be turned into. Measured at 133,018
    cells: 303 MB peak for a 406 kB answer.

    Three outcomes. Disjoint: nothing below can be kept under any mode,
    because every descendant square is a subset of this one. Wholly
    inside: everything below is kept under every mode, so the node is
    emitted as itself when compacting and with a shared expansion
    otherwise. Straddling: subdivide, or apply the leaf predicate if
    this is already the target.
    """
    from shapely.geometry import box

    square = box(u0, v0, u0 + side, v0 + side)
    if not prepared.intersects(square):
        return None
    if prepared.contains(square):
        if compact or level == target:
            budget.charge(1)
            return code
        budget.charge(_leaves_between(level, target))
        return code + _expansion(level, target)
    if level == target:
        if _accept_leaf(prepared, containment, u0, v0, side):
            budget.charge(1)
            return code
        return None
    step = level + 1
    divisor = linear_refinement_ratio(step)
    child = side / divisor
    parts: list[str] = []
    for row in range(divisor):
        for column in range(divisor):
            fragment = _fill_node(
                prepared,
                u0 + column * child,
                v0 + row * child,
                child,
                step,
                target,
                containment,
                _child_code(row, column, step),
                compact,
                budget,
            )
            if fragment is not None:
                parts.append(fragment)
    if not parts:
        return None
    return f"{code}{DESCENT_OPEN}{SIBLING_SEPARATOR.join(parts)}{DESCENT_CLOSE}"


def _count_node(
    prepared: "PreparedGeometry",
    u0: float,
    v0: float,
    side: float,
    level: int,
    target: int,
) -> int:
    """Accumulate the cell count under one node without naming any cell.

    The same three outcomes as :func:`_fill_node`, with the wholly-inside
    case answered by arithmetic rather than by enumeration: the number of
    target cells under a node is the product of the refinement ratios
    between the two levels, so a node inside the geometry contributes its
    whole subtree in one step and is never descended.

    Memory is the recursion stack, whose depth is the resolution
    difference and therefore at most twelve. Nothing accumulates per
    cell, which is what lets the count run at resolution 13 where naming
    the cells could not.
    """
    from shapely.geometry import Point, box

    square = box(u0, v0, u0 + side, v0 + side)
    if not prepared.intersects(square):
        return 0
    if prepared.contains(square):
        return _leaves_between(level, target)
    if level == target:
        centre = Point(u0 + side / 2.0, v0 + side / 2.0)
        return 1 if prepared.contains(centre) else 0
    step = level + 1
    divisor = linear_refinement_ratio(step)
    child = side / divisor
    total = 0
    for row in range(divisor):
        for column in range(divisor):
            total += _count_node(
                prepared,
                u0 + column * child,
                v0 + row * child,
                child,
                step,
                target,
            )
    return total


def _prepare(
    geometry: "BaseGeometry", resolution: int, densify: bool
) -> list[tuple[str, "_LatticeView"]]:
    """Screen, densify, project and shear a geometry, one part per quadrant.

    Returns ``(quadrant, lattice view)`` pairs ready for descent.
    """

    if geometry.is_empty:
        return []
    if crosses_antemeridian(geometry):
        raise AntemeridianError(
            "geometry crosses 180 degrees longitude outside an extension "
            "zone; split it at the antemeridian or express it with "
            "longitudes past 180 inside a defined zone"
        )
    if densify and geometry.geom_type in {"Polygon", "MultiPolygon"}:
        geometry = _densify_any(geometry, _auto_segment(resolution))
    plane = _project(geometry)
    out: list[tuple[str, "_LatticeView"]] = []
    for quadrant, piece in _quadrant_pieces(plane):
        out.append((quadrant, _LatticeView(_to_lattice(piece, quadrant))))
    return out


_T = TypeVar("_T")
_R = TypeVar("_R")


class _LatticeView:
    """A per-thread prepared copy of one quadrant's lattice geometry.

    Shapely's prepared geometry is **not** safe to share across threads.
    A prepared geometry builds its spatial index lazily on the first
    predicate call and mutates the underlying GEOS object while doing so,
    so two threads entering it together corrupt that structure and the
    interpreter dies with a segmentation fault rather than an exception.
    Found by ``test_threads_and_serial_agree``, which is why that test
    compares the two paths instead of only exercising one.

    Each thread therefore rebuilds the geometry from its WKB and prepares
    its own. Nothing is shared but the bytes, which are immutable, and
    the copy costs one deserialisation per thread rather than per cell.
    """

    def __init__(self, lattice: "BaseGeometry") -> None:
        self._wkb = lattice.wkb
        self.geometry = lattice
        self._local = threading.local()

    @property
    def prepared(self) -> "PreparedGeometry":
        """This thread's prepared geometry, built on first use."""
        from shapely import wkb as _wkb
        from shapely.prepared import prep

        cached = getattr(self._local, "prepared", None)
        if cached is None:
            cached = prep(_wkb.loads(self._wkb))
            self._local.prepared = cached
        return cached


def _map_jobs(work: Sequence[_T], call: Callable[[_T], _R], n_jobs: int) -> list[_R]:
    """Run ``call`` over ``work``, on a thread pool when asked for one.

    Threads rather than processes: the geometry predicates are Shapely's,
    which release the interpreter lock, and the work items are small
    enough that pickling a polygon to a child would dominate. Thread
    safety comes from :class:`_LatticeView` giving each worker its own
    prepared copy, not from sharing one.
    """
    if n_jobs == 1:
        return [call(item) for item in work]
    with ThreadPoolExecutor(max_workers=n_jobs) as pool:
        return list(pool.map(call, work))


def polyfill(
    geometry: "BaseGeometry",
    resolution: int,
    containment: Containment = "center",
    compact: bool = False,
    n_jobs: int = 1,
) -> str:
    """Rasterise a geometry into a compositional index.

    Descends the hierarchy over the sinusoidal plane rather than testing
    cells one by one, so the cost tracks the boundary rather than the
    area.

    ``center`` keeps cells whose centre falls inside, ``intersects``
    keeps every cell touching the geometry, and ``contains`` keeps only
    cells wholly inside. The three nest, ``contains`` inside ``center``
    inside ``intersects``, and they nest by construction rather than by
    coincidence: see :func:`_accept_leaf`.

    **The centre, not the anchor.** A cell's anchor is a vertex and lies
    on its own border, so an anchor test would award a cell by a point it
    shares with three neighbours, and ``contains`` would stop being a
    subset of ``center`` for every cell whose anchor sits on the outline.
    The centre is interior to the cell and has neither problem.

    **Densification is applied, not assumed.** A straight edge on the
    plane is not a geodesic, so an undensified long edge fills the wrong
    cells in between and the result is a plausible, wrong cell set with
    no error raised. Areal input is therefore densified first, at the
    threshold :func:`_auto_segment` derives from the target resolution.
    Pass an already densified geometry and the step is idempotent.

    **Restriction.** Three families are refused rather than filled: the
    prime-meridian column, the last lattice column of a row, and the
    polar row. In all three the cell is not the sheared square the
    descent tests, so the containment chain would no longer hold. The
    screen runs at resolution 1 and is conservative.

    Args:
        geometry: A Shapely geometry in EPSG:4326.
        resolution: Target resolution level, 1 to 13.
        containment: Predicate deciding whether a cell is kept.
        compact: Return a mixed-resolution compacted index instead of a
            uniform one.
        n_jobs: Worker count; above 1 spreads base cells over threads.

    Returns:
        A compositional index string covering the geometry.

    Raises:
        AntemeridianError: If the geometry crosses 180 degrees outside an
            extension zone.
        UnsupportedGeometryTypeError: On unsupported geometry types.
        NonExistentCellError: If the geometry reaches the prime-meridian
            column or a border-absorbing column.
        DomainError: If the geometry reaches the polar row.
        GeometryError: If the fill exceeds :data:`MAX_FILL_CELLS`.
    """
    _check_resolution(resolution)
    _check_jobs(n_jobs)
    if containment not in ("center", "intersects", "contains"):
        raise ValueError(
            f"containment must be 'center', 'intersects' or 'contains', "
            f"got {containment!r}"
        )

    budget = _Budget()
    roots: list[str] = []
    for quadrant, view in _prepare(geometry, resolution, densify=True):

        def _one(base: tuple[int, int], _v: "_LatticeView" = view) -> str | None:
            column, row = base
            code = f"{column:0{RES1_DIGITS}d}{RES1_SEPARATOR}{row:0{RES1_DIGITS}d}"
            return _fill_node(
                _v.prepared,
                (column + row) * _L1,
                row * _L1,
                _L1,
                1,
                resolution,
                containment,
                code,
                compact,
                budget,
            )

        # Sorted by column then row so that the siblings come out in the
        # order the index is read in. compose() preserves the order it is
        # given rather than imposing one, so the ordering has to be here.
        cells = sorted(_base_cells(view.prepared, view.geometry, quadrant, resolution))
        parts = [f for f in _map_jobs(cells, _one, n_jobs) if f is not None]
        if parts:
            joined = SIBLING_SEPARATOR.join(parts)
            roots.append(f"{quadrant}{DESCENT_OPEN}{joined}{DESCENT_CLOSE}")

    if not roots:
        raise GeometryError(
            "geometry covers no cell at this resolution; it may be empty "
            "after projection or narrower than one cell under the chosen "
            "containment mode"
        )
    return SIBLING_SEPARATOR.join(roots)


def count_internal_cells(polygon: "Polygon", resolution: int, n_jobs: int = 1) -> int:
    """Count the cells a polygon covers, without naming any of them.

    Fast path over :func:`polyfill` with ``containment="center"``: the
    count is accumulated during descent, so nothing is materialised and
    the call runs at resolution 13 where the index could not be built.

    Since every cell carries the same area, the count times
    :func:`itacart.resolutions.nominal_cell_area` gives a
    distortion-free area, which is the property the paper builds
    tokenization on.

    **Why the centre and not strict containment.** Counting only the
    cells wholly inside undercounts by the ring of cells the outline
    crosses, about ``perimeter * side / 2`` in area terms. For a
    square kilometre at resolution 7 that is two per cent. Under the
    centre rule the cells the outline crosses are kept or dropped
    according to which side their centre falls, and the two errors
    cancel to first order. The origin counts the same way.

    **What the residual error is, and what it is not.** The count is
    exact for every cell the outline does not cross, so the area it
    implies can only be wrong about the crossed ones:

    ``abs(n * a - A) / A <= P * s / A`` for perimeter ``P``, cell side
    ``s`` and area ``A``. That bound is arithmetic and holds for every
    parcel.

    Where inside the bound a given parcel lands is not predictable.
    Measured over 24 phases at each of seven parcel sizes, three
    estimators of the scaling exponent disagree outright: 0.61 from the
    RMS, -0.47 from the median, 0.00 from the ninetieth percentile.
    There is no power law to fit, because the residual depends on how
    that outline happens to sit against that lattice.

    This matters for how the resolution is chosen. It is not chosen from
    the parcel: in cadastral work it follows from the mapping scale and
    the positional tolerance the applicable standard attaches to it, and
    it is fixed before any parcel is measured. The parcel's size and
    perimeter then determine the residual, which is an outcome to be
    reported against the bound above, not a target to refine toward.

    **Restriction.** The same three families :func:`polyfill` refuses are
    refused here, and for the same reason: outside the parallelogram
    interior a cell does not carry the nominal area, so the product would
    not be an area.

    Provenance: ``itacart_core/cell_filling.py``
    (``polygon_to_cells_count``).

    Args:
        polygon: A Shapely polygon in EPSG:4326.
        resolution: Target resolution level, 1 to 13.
        n_jobs: Worker count; above 1 spreads base cells over threads.

    Returns:
        Number of cells whose centre lies inside the polygon.

    Raises:
        AntemeridianError: If the polygon crosses 180 degrees outside an
            extension zone.
        NonExistentCellError: If the polygon reaches the prime-meridian
            column or a border-absorbing column.
        DomainError: If the polygon reaches the polar row.
    """
    _check_resolution(resolution)
    _check_jobs(n_jobs)
    total = 0
    for quadrant, view in _prepare(polygon, resolution, densify=True):

        def _one(base: tuple[int, int], _v: "_LatticeView" = view) -> int:
            column, row = base
            return _count_node(
                _v.prepared, (column + row) * _L1, row * _L1, _L1, 1, resolution
            )

        cells = _base_cells(view.prepared, view.geometry, quadrant, resolution)
        total += sum(_map_jobs(cells, _one, n_jobs))
    return total


# --------------------------------------------------------------------------
# Vertex representation
# --------------------------------------------------------------------------


def _dedupe_consecutive(cells: Sequence[str], cyclic: bool) -> list[str]:
    """Collapse runs of the same index, optionally treating the list as a ring.

    Only runs. A repeat that is not adjacent is left alone: at the target
    resolution two consecutive vertices in the same cell are the same
    surveyed point written twice, while a vertex revisited later in the
    ring is a self-touching outline, which is a genuine pathology and the
    caller has to be able to see it.

    With ``cyclic`` the wrap-around pair is treated as adjacent too,
    since a ring's last vertex is followed by its first.
    """
    if not cells:
        return []
    out = [cells[0]]
    for cell in cells[1:]:
        if cell != out[-1]:
            out.append(cell)
    if cyclic and len(out) >= 2 and out[0] == out[-1]:
        out.pop()
    return out


def _ring_cells(
    coords: Sequence[tuple[float, ...]], resolution: int, dedupe: bool
) -> list[str]:
    """Map one closed ring to cells, dropping the repeated closing vertex."""
    points = list(coords)
    if len(points) >= 2 and points[0] == points[-1]:
        points = points[:-1]
    cells = [geo_to_cell(point[0], point[1], resolution) for point in points]
    return _dedupe_consecutive(cells, cyclic=True) if dedupe else cells


def vertex_to_cell(
    geometry: "Point | LineString | Polygon",
    resolution: int,
    dedupe_consecutive: bool = True,
) -> list[str]:
    """Map geometry vertices to cells, preserving sequence.

    Order is topological, not sorted: it is what allows the original
    geometry to be reconstructed. Rings keep their winding, and holes
    follow the exterior.

    Consecutive vertices landing in the same cell are collapsed by
    default, since at resolution 13 that reflects survey precision rather
    than distinct corners. Non-consecutive repeats are kept, being a
    genuine self-touching pathology the caller should see.

    A ring is treated as cyclic for the purposes of collapsing, so a ring
    whose last vertex falls in the same cell as its first loses the last
    one rather than closing on a duplicate. A LINESTRING is not cyclic
    and keeps both ends even when they coincide.

    Provenance: ``cadastral_processor/vertex_extractor.py``.

    Args:
        geometry: A Shapely geometry in EPSG:4326.
        resolution: Target resolution level, 1 to 13.
        dedupe_consecutive: Collapse consecutive duplicates.

    Returns:
        Atomic index strings in traversal order.

    Raises:
        UnsupportedGeometryTypeError: On unsupported geometry types.
    """
    from shapely.geometry import LineString, Point, Polygon

    _check_resolution(resolution)
    if isinstance(geometry, Point):
        return [geo_to_cell(geometry.x, geometry.y, resolution)]
    if isinstance(geometry, Polygon):
        out = _ring_cells(
            list(geometry.exterior.coords), resolution, dedupe_consecutive
        )
        for hole in geometry.interiors:
            out.extend(_ring_cells(list(hole.coords), resolution, dedupe_consecutive))
        return out
    if isinstance(geometry, LineString):
        cells = [geo_to_cell(x, y, resolution) for x, y in geometry.coords]
        return _dedupe_consecutive(cells, cyclic=False) if dedupe_consecutive else cells
    raise UnsupportedGeometryTypeError(
        f"cannot take vertices of a {geometry.geom_type}; supported types "
        "are Point, LineString and Polygon"
    )


def cells_to_geometry(
    cells: list[str], geometry_type: str = "Polygon"
) -> "BaseGeometry":
    """Rebuild a geometry from an ordered vertex cell list.

    Inverse of :func:`vertex_to_cell`. Reconstruction lands on cell
    anchors, so it is exact only to the resolution used: at resolution 13
    that is 1 cm.

    **The inverse is exact for a hole-free geometry only.**
    :func:`vertex_to_cell` returns one flat sequence with the holes
    appended to the exterior, and a flat sequence does not say where the
    exterior ended. A caller needing holes keeps the rings apart and
    calls this once per ring.

    Args:
        cells: Atomic index strings in traversal order.
        geometry_type: OGC SFA type to build: ``Point``, ``LineString``,
            ``LinearRing``, ``Polygon`` or ``MultiPoint``.

    Returns:
        A Shapely geometry in EPSG:4326.

    Raises:
        UnsupportedGeometryTypeError: On unsupported types.
        GeometryError: If the cell count cannot make the type asked for.
    """
    from shapely.geometry import LinearRing, LineString, MultiPoint, Point, Polygon

    if not cells:
        raise GeometryError("cannot rebuild a geometry from an empty cell list")
    coords: list[tuple[float, float]] = []
    for cell in cells:
        anchor = cell_to_anchor(cell)
        if not isinstance(anchor, tuple):
            raise GeometryError(
                f"{cell!r} names more than one cell; cells_to_geometry takes "
                "atomic indices in traversal order"
            )
        coords.append((float(anchor[0]), float(anchor[1])))

    if geometry_type == "Point":
        if len(coords) != 1:
            raise GeometryError(f"a Point needs exactly one cell, got {len(coords)}")
        return Point(coords[0])
    if geometry_type == "MultiPoint":
        return MultiPoint(coords)
    if geometry_type == "LineString":
        if len(coords) < 2:
            raise GeometryError(
                f"a LineString needs at least two cells, got {len(coords)}"
            )
        return LineString(coords)
    if geometry_type in ("Polygon", "LinearRing"):
        if len(coords) < 3:
            raise GeometryError(
                f"a {geometry_type} needs at least three cells, got {len(coords)}"
            )
        ring = coords + [coords[0]]
        return LinearRing(ring) if geometry_type == "LinearRing" else Polygon(ring)
    raise UnsupportedGeometryTypeError(
        f"cannot rebuild a {geometry_type!r}; supported types are Point, "
        "MultiPoint, LineString, LinearRing and Polygon"
    )


# --------------------------------------------------------------------------
# Densification
# --------------------------------------------------------------------------


EDGE_MODELS: tuple[str, ...] = ("WGS84_GEODESIC",)
"""Edge interpretations :func:`densify_segment` understands."""


def _check_threshold(max_segment_m: float) -> None:
    """Reject a threshold that is not a positive finite length."""
    if not isinstance(max_segment_m, (int, float)) or isinstance(max_segment_m, bool):
        raise DensificationError(
            f"max_segment_m must be a number, got {type(max_segment_m).__name__}"
        )
    if not math.isfinite(max_segment_m) or max_segment_m <= 0.0:
        raise DensificationError(
            f"max_segment_m must be a positive finite length in metres, "
            f"got {max_segment_m!r}"
        )


def densify_segment(
    p1: tuple[float, float],
    p2: tuple[float, float],
    max_segment_m: float,
    edge_model: str = "WGS84_GEODESIC",
) -> list[tuple[float, float]]:
    """Densify a single segment.

    The building block of :func:`densify_orthodromic`, exposed for
    callers working segment by segment such as open LINESTRING handling.

    The span is measured with the Vincenty inverse and the intermediate
    points are placed with the Vincenty direct along the forward azimuth,
    equally spaced in geodesic arc length. Nothing is measured on the
    plane, because the projection preserves area and not direction.

    ``n_segments = floor(d / max_segment_m) + 1``, so every resulting leg
    is strictly shorter than the threshold whenever any subdivision
    happens at all. That is what makes the operation idempotent: a second
    pass measures each leg, finds it under the threshold, and leaves it.

    Provenance: ``itacart_core/geometry_blob.py`` (``densify_segment``).

    Args:
        p1: ``(lon, lat)`` of the start point.
        p2: ``(lon, lat)`` of the end point.
        max_segment_m: Longest segment to leave undensified, in metres.
        edge_model: Edge interpretation, currently ``"WGS84_GEODESIC"``.

    Returns:
        ``(lon, lat)`` pairs including both endpoints.

    Raises:
        DensificationError: If the threshold is not positive and finite,
            or the edge model is unknown.
    """
    _check_threshold(max_segment_m)
    if edge_model not in EDGE_MODELS:
        raise DensificationError(
            f"unknown edge model {edge_model!r}; supported models are "
            f"{', '.join(EDGE_MODELS)}"
        )
    lon1, lat1 = float(p1[0]), float(p1[1])
    lon2, lat2 = float(p2[0]), float(p2[1])
    distance, azimuth = inverse_geodesic(lon1, lat1, lon2, lat2)
    if distance <= 0.0:
        return [(lon1, lat1), (lon2, lat2)]
    pieces = int(math.floor(distance / max_segment_m)) + 1
    if pieces == 1:
        return [(lon1, lat1), (lon2, lat2)]
    step = distance / pieces
    out = [(lon1, lat1)]
    for index in range(1, pieces):
        longitude, latitude = direct_geodesic(lon1, lat1, azimuth, step * index)
        out.append((_on_the_branch_of(lon1, longitude), latitude))
    out.append((lon2, lat2))
    return out


def _on_the_branch_of(reference: float, longitude: float) -> float:
    """Move ``longitude`` onto the 360-degree branch of ``reference``.

    :func:`itacart.geodesy.direct_geodesic` normalises what it returns to
    ``(-180, 180]``, which is the right answer to the question it is
    asked and the wrong one for densification. A segment written from
    179.9 to 180.3 -- the natural way to describe a footprint inside an
    extension zone, where the domain reaches past the antemeridian --
    has interior points at 180.1, and normalising those to -179.9 folds
    the ring back across the globe. The result self-intersects, and the
    quadrant clip downstream fails inside GEOS with a topology error
    rather than anywhere the package can explain.

    Densification therefore keeps the branch the caller wrote. Segments
    that genuinely wrap are refused upstream by
    :func:`itacart.boundary.crosses_antemeridian`, so no segment reaching
    here spans more than half the globe and the branch is unambiguous.
    """
    return longitude - 360.0 * round((longitude - reference) / 360.0)


def _densify_ring(
    coords: Sequence[tuple[float, ...]], max_segment_m: float
) -> list[tuple[float, float]]:
    """Densify one ring, returned open with no repeated closing vertex."""
    points = [(float(point[0]), float(point[1])) for point in coords]
    if len(points) >= 2 and points[0] == points[-1]:
        points = points[:-1]
    if len(points) < 2:
        return points
    out: list[tuple[float, float]] = []
    for index, start in enumerate(points):
        end = points[(index + 1) % len(points)]
        out.extend(densify_segment(start, end, max_segment_m)[:-1])
    return out


def densify_orthodromic(polygon: "Polygon", max_segment_m: float = 1000.0) -> "Polygon":
    """Insert intermediate vertices along geodesics.

    Applied to the exterior ring and every hole, with
    ``n_segments = floor(d_geo / max_segment_m) + 1`` and points equally
    spaced in geodesic distance.

    Needed because a straight line on the sinusoidal plane is not a
    geodesic on the ellipsoid; without densification a long edge would
    fill the wrong cells in between.

    Holes are densified too. The origin drops them, which is silent: the
    filled cell set comes back plausible and includes the hole.

    Provenance: ``itacart_core/densification.py``.

    Args:
        polygon: A Shapely polygon in EPSG:4326.
        max_segment_m: Longest segment to leave undensified, in metres.

    Returns:
        A densified polygon in EPSG:4326.

    Raises:
        DensificationError: If a segment cannot be densified.
        TypeError: If the argument is not a Shapely polygon.
    """
    from shapely.geometry import Polygon

    if not isinstance(polygon, Polygon):
        raise TypeError(f"expected a shapely Polygon, got {type(polygon).__name__}")
    _check_threshold(max_segment_m)
    if polygon.is_empty:
        return polygon
    return Polygon(
        _densify_ring(list(polygon.exterior.coords), max_segment_m),
        [_densify_ring(list(hole.coords), max_segment_m) for hole in polygon.interiors],
    )


def _densify_any(geometry: "BaseGeometry", max_segment_m: float) -> "BaseGeometry":
    """Densify a polygon or a multipolygon, leaving anything else alone."""
    from shapely.geometry import MultiPolygon, Polygon

    if isinstance(geometry, Polygon):
        return densify_orthodromic(geometry, max_segment_m)
    if isinstance(geometry, MultiPolygon):
        return MultiPolygon(
            [densify_orthodromic(part, max_segment_m) for part in geometry.geoms]
        )
    return geometry


# --------------------------------------------------------------------------
# Canonicalization
# --------------------------------------------------------------------------


def _min_rotation(ring: Sequence[str]) -> int:
    """Start index of the lexicographically smallest rotation of ``ring``.

    Booth's algorithm, linear in the ring length. A naive scan over the
    smallest vertex breaks on a ring that repeats its minimum, which a
    self-touching outline does; Booth compares whole rotations and orders
    those cases too.

    Provenance: ``itacart_core/geometry_blob.py``
    (``_booth_min_rotation``).
    """
    length = len(ring)
    if length == 0:
        return 0
    doubled = list(ring) * 2
    failure = [-1] * (2 * length)
    start = 0
    for position in range(1, 2 * length):
        candidate = doubled[position]
        offset = failure[position - start - 1]
        while offset != -1 and candidate != doubled[start + offset + 1]:
            if candidate < doubled[start + offset + 1]:
                start = position - offset - 1
            offset = failure[offset]
        if candidate != doubled[start + offset + 1]:
            if candidate < doubled[start + offset + 1]:
                start = position
            failure[position - start] = -1
        else:
            failure[position - start] = offset + 1
    return start


def canonicalize_rings(rings: list[list[str]]) -> list[list[str]]:
    """Normalise rings to a single spelling per geometry.

    Rings are rotated to their minimum lexicographic cyclic rotation, so
    the same ring starting at a different vertex canonicalises to the
    same sequence. Rotation applies to closed rings only; LINESTRING is
    directional and is left as given, which is why this function takes
    rings and not lines.

    Winding is **not** touched. A ring's direction distinguishes an
    exterior from a hole and is information the caller put there;
    rotating is a change of spelling, reversing would be a change of
    meaning. Ordering the rings themselves is likewise the caller's,
    since the exterior comes first by contract.

    The ordering key is the index string. Two cells of the same
    resolution have strings of the same length, so string order is a
    total order over the ring and the rotation is unique.

    This is what makes a geometry content-addressable, and therefore what
    makes hashing it meaningful.

    Provenance: ``itacart_core/geometry_blob.py``
    (``_canonicalize_polygon_rings``, ``_booth_min_rotation``).

    Args:
        rings: Rings as lists of atomic index strings, exterior first.

    Returns:
        The canonicalised rings, in input order.
    """
    out: list[list[str]] = []
    for ring in rings:
        if len(ring) < 2:
            out.append(list(ring))
            continue
        start = _min_rotation(ring)
        out.append(list(ring[start:]) + list(ring[:start]))
    return out
