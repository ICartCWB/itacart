"""Neighbourhood and adjacency on the parallelogram lattice.

Completes OGC DGGS Core requirement 17. Adjacency is derived from the index
strings alone; coordinates are never consulted to *find* a neighbour, only
:func:`itacart.boundary.is_valid_cell` is asked whether the neighbour exists,
and that predicate is closed form in the cosine of the latitude.

The cell is a parallelogram sheared by 45 degrees, with its base equal to its
height, so a step toward the pole also moves one column toward the prime
meridian. On the parallels plane the step measures ``dx = -10 km`` and
``dy = +10 km`` exactly. In index space the four edge-adjacent steps are the
same arithmetic in all four quadrants, because quadrant mirroring keeps ``X``
growing away from the prime meridian and ``Y`` growing toward the pole::

    N  (X - 1, Y + 1)        E  (X + 1, Y)
    S  (X + 1, Y - 1)        W  (X - 1, Y)

Section 3.1 of the paper attaches the plus and minus one to ``Y`` for the
vertical step. The magnitude and the sign are right; the component is not.
Reading it literally returns the cell diagonally above, which shares a single
vertex rather than an edge, and every symmetry test still passes because a
diagonal is symmetric too.

The four vertex-adjacent steps are the sums of the edge steps, so ``NE`` is
``(X, Y + 1)``, ``NW`` is ``(X - 2, Y + 1)``, ``SE`` is ``(X + 2, Y - 1)`` and
``SW`` is ``(X, Y - 1)``.

Grid distance is measured on the lattice, not on the ellipsoid. Callers coming
from H3 tend to assume otherwise; use :func:`itacart.geodesy.inverse_geodesic`
on the anchors or centroids for metric separation.
"""

from __future__ import annotations

from collections.abc import Callable
from functools import lru_cache
from typing import Literal, cast

from shapely.affinity import translate
from shapely.geometry import LineString, Polygon
from shapely.geometry.base import BaseGeometry
from shapely.ops import snap

from .boundary import (
    ZONE_ROWS,
    absorbs_border,
    cell_shape,
    is_valid_cell,
    last_lattice_column,
)
from .cells import cell_to_boundary
from .exceptions import DomainError, InvalidIndexError, ResolutionError
from .hierarchy import _parent_cell, get_children
from .index import decompose, join_components, split_components
from .resolutions import cell_size, get_resolution

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

#: Lattice displacement of each direction, as ``(dX, dY)``. The cardinal four
#: are edge steps; the diagonal four are their sums and share only a vertex.
LATTICE_STEP: dict[str, tuple[int, int]] = {
    "N": (-1, 1),
    "S": (1, -1),
    "E": (1, 0),
    "W": (-1, 0),
    "NE": (0, 1),
    "NW": (-2, 1),
    "SE": (2, -1),
    "SW": (0, -1),
}

#: In-parent displacement of each direction, as ``(dcolumn, drow)``. The
#: refinement grid of a parent is an ordinary square lattice: the parent's own
#: frame shears with it, so the poleward step keeps the column here, unlike the
#: resolution-1 step, which is expressed in absolute columns.
GRID_STEP: dict[str, tuple[int, int]] = {
    "N": (0, 1),
    "S": (0, -1),
    "E": (1, 0),
    "W": (-1, 0),
    "NE": (1, 1),
    "NW": (-1, 1),
    "SE": (1, -1),
    "SW": (-1, -1),
}

_COMPASS: dict[tuple[int, int], Direction] = {
    step: cast(Direction, name) for name, step in GRID_STEP.items()
}
_ALPHABET = "ABCDE"
_QUATERNARY_SIDE = 2
_QUINARY_SIDE = 5

#: The four edge steps, and the four vertex steps, as ordered tuples.
_CARDINALS: tuple[Direction, ...] = ("N", "S", "E", "W")
_DIAGONALS: tuple[Direction, ...] = ("NE", "NW", "SE", "SW")

_MIRRORED_EAST_WEST = frozenset({"NW", "SW"})
_MIRRORED_NORTH_SOUTH = frozenset({"SE", "SW"})
_FLIP_EAST_WEST = str.maketrans({"E": "W", "W": "E"})
_FLIP_NORTH_SOUTH = str.maketrans({"N": "S", "S": "N"})

_ACROSS_MERIDIAN = {"NE": "NW", "NW": "NE", "SE": "SW", "SW": "SE"}
_ACROSS_EQUATOR = {"NE": "SE", "SE": "NE", "NW": "SW", "SW": "NW"}
_POLAR_ROW = 1000
_BASE_SIDE_M = 10_000.0

#: Snapping distance for ring comparison, in degrees, and the length above
#: which a shared boundary counts as an edge rather than a point.
_SNAP_TOLERANCE = 1e-9
_EDGE_TOLERANCE = 1e-8

#: Columns of slack around the cell's own column, and the width of the
#: bands searched at the prime meridian and at the antemeridian seam.
_WINDOW = 3
_SEAM_REACH = 8


def _res1_parts(cell: str) -> tuple[str, int, int]:
    """Split an atomic resolution-1 index into quadrant, column and row."""
    components = split_components(cell)
    if len(components) != 2:
        raise ResolutionError(
            f"{cell!r} is not a resolution-1 cell; the lattice step is defined "
            "on the XXXX/YYYY pair"
        )
    column, _, row = components[1].partition("/")
    return components[0], int(column), int(row)


def _spell(quadrant: str, column: int, row: int) -> str:
    """Render a resolution-1 index from its parts."""
    return f"{quadrant}({column:04d}/{row:04d})"


def _last_column(quadrant: str, row: int) -> int:
    """Greatest addressable column of a resolution-1 row."""
    return last_lattice_column(quadrant, row, _BASE_SIDE_M)


def _zone_limit_rows() -> frozenset[int]:
    """Rows on either side of an extension zone's latitude limits.

    A zone is realised on whole resolution-1 rows, so its two limits are four
    rows: the last row outside and the first row inside, at each end. Across
    those limits one quadrant is extended past 180 degrees while its mirror is
    truncated, and a cell's poleward or equatorward edge is then shared with
    cells of the *other* quadrant. The neighbour is real but the lexical step
    cannot name it, and the count of cells involved does not follow from the
    difference of the two row limits, so these cells are dispatched here and
    resolved geometrically rather than guessed at.
    """
    rows: set[int] = set()
    for first, last in ZONE_ROWS.values():
        rows.update({first - 1, first, last, last + 1})
    return frozenset(rows)


ZONE_LIMIT_ROWS = _zone_limit_rows()


def _lattice_direction(quadrant: str, direction: Direction) -> Direction:
    """Translate a true cardinal direction into a quadrant's lattice frame.

    The four quadrants reuse one piece of arithmetic by mirroring: the
    lattice counts columns away from the prime meridian and rows away
    from the equator, whichever way that runs on the globe. So a lattice
    ``"E"`` is cardinal east in NE and SE and cardinal **west** in NW and
    SW, and a lattice ``"N"`` is cardinal north in NE and NW and cardinal
    **south** in SE and SW. In SW both axes are mirrored, and a lattice
    ``"NE"`` points cardinally south-west.

    That mirroring is an implementation convenience and not a property of
    the grid, so it stays inside. The public directions are true compass
    points, and this is the one place the two frames meet.

    The translation is its own inverse, which is what makes the public
    function symmetric wherever the lattice one was: applying it at the
    origin and again at the answer's quadrant returns the original step.
    """
    text: str = direction
    if quadrant in _MIRRORED_EAST_WEST:
        text = text.translate(_FLIP_EAST_WEST)
    if quadrant in _MIRRORED_NORTH_SOUTH:
        text = text.translate(_FLIP_NORTH_SOUTH)
    return cast(Direction, text)


def get_neighbor(index: str, direction: Direction) -> str | None | list[str | None]:
    """Single adjacent cell in a given direction.

    Directions are true compass points in every quadrant: ``"N"`` is
    cardinal north and ``"E"`` is cardinal east whether the cell sits east
    or west of the prime meridian, north or south of the equator. The
    cardinal four share an edge with the origin; the diagonal four share
    only a vertex.

    Inside, the lattice counts columns away from the meridian and rows
    away from the equator, which is what lets one piece of arithmetic
    serve four quadrants. :func:`_lattice_direction` is the seam.

    **Reciprocity is not promised.** This is a directional selector, not
    an involution. Where a row's last column shifts between rows, the
    step back from the answer can land one column off the origin.
    Enumerated over the four border families -- meridian, last column,
    equator and polar row, every cell of each, in four quadrants --
    6,338 such pairs out of 102,090. Every one of them was classified,
    and every one is the same case: the answer genuinely touches its
    origin, and the origin is genuinely among the answer's neighbours.
    None is a step that missed, and none is the tie-break of
    :func:`cell_to_edges` choosing differently. What varies is the
    length of the row.

    Geometric adjacency *is* symmetric -- :func:`are_neighbor_cells`
    answers the same from either side -- and that is the invariant to
    compose on.

    **The diagonals are lattice steps.** ``"NE"`` moves one column and
    one row, which is geographically north-east everywhere except at the
    last column of a row, where the shear can carry the answer the other
    way in longitude: 4,382 of 25,572 diagonal pairs in that family, and
    none at all in the other three. The cardinal four are true compass
    points everywhere, with zero exceptions in 102,090 enumerated
    pairs.

    Derivation is lexical: the direction's displacement is added to the
    ``XXXX/YYYY`` pair and nothing else is consulted. Existence is a separate
    question, delegated to :func:`itacart.boundary.is_valid_cell`, and when
    the lexical target does not exist the step is handed to :func:`deflect`.

    Args:
        index: Compositional index string.
        direction: One of the eight lattice directions.

    Returns:
        The neighbour index, or ``None`` where no cell exists on that side,
        for a single cell; or a positionally aligned list.

    Raises:
        ValueError: If ``direction`` is not one of the eight.
        ResolutionError: If the cell is not at resolution 1. Descent through
            the refinement alphabets is not yet implemented here.
        DomainError: If the cell sits in the geometric exception set, where
            the step is real but no lexical rule names it. See :func:`deflect`.
    """
    if direction not in LATTICE_STEP:
        raise ValueError(
            f"{direction!r} is not a lattice direction; expected one of "
            f"{sorted(LATTICE_STEP)}"
        )
    atoms = decompose(index)
    results = [
        _neighbor_of_atom(atom, _lattice_direction(atom[:2], direction))
        for atom in atoms
    ]
    return results[0] if len(results) == 1 else results


def _lexical_target(cell: str, direction: Direction) -> tuple[str, int, int]:
    """The index the direction's displacement spells, and its raw components.

    Pure arithmetic on the ``XXXX/YYYY`` pair. Nothing here asks whether the
    result exists, which is what lets the derivation be tested with every
    existence predicate replaced by one that raises.
    """
    quadrant, column, row = _res1_parts(cell)
    shift_x, shift_y = LATTICE_STEP[direction]
    target_x, target_y = column + shift_x, row + shift_y
    return _spell(quadrant, target_x, target_y), target_x, target_y


def _neighbor_of_atom(cell: str, direction: Direction) -> str | None:
    """Resolve one cell's neighbour: lexically first, measured only if it must.

    Outside the exception set the lexical step answers, and :func:`deflect`
    normalises it when it leaves the quadrant's range. Inside the exception
    set the step is not enough — a trapezoid carries several cells on one
    side, a triangle's base is split between quadrants — so the neighbour sets
    are measured and the one lying in the asked direction is returned. Where
    several do, the one sharing the longest boundary wins, which is stated
    rather than hidden: :func:`cell_to_edges` returns all of them.
    """
    resolution = get_resolution(cell)
    if resolution < 1:
        raise ResolutionError(
            f"{cell!r} is a quadrant; a quadrant has no lattice neighbour"
        )
    if _is_exceptional(cell) or (resolution > 1 and _has_exceptional_ancestor(cell)):
        return _by_contact(cell, direction)
    try:
        return _verified_lexical(cell, direction, resolution)
    except DomainError:
        # The index stopped being able to represent the step: the neighbouring
        # parent does not tile a square grid, so the wrapped code would name
        # the wrong sibling. This is the line the phase drew — structure is
        # resolved lexically, and geometry enters only past it.
        return _by_contact(cell, direction)


def _verified_lexical(cell: str, direction: Direction, resolution: int) -> str | None:
    """The lexical step, refused when it leaves the base cell without touching.

    The single entry point for the index route. Both callers go through it —
    :func:`_neighbor_of_atom` for one direction and :func:`_contacts` for all
    eight — because verifying in one and not the other lets a phantom into the
    neighbour set: a step across the seam can spell a real cell that the
    origin does not touch, and an unverified set would then hold it while
    missing the cell that is really there.
    """
    answer = _lexical_neighbor(cell, direction, resolution)
    if answer is not None and _changes_frame(cell, answer):
        _require_contact(cell, answer, direction)
    return answer


def _lexical_neighbor(cell: str, direction: Direction, resolution: int) -> str | None:
    """The step, taken from the index alone, with no measured fallback.

    Kept apart from :func:`_neighbor_of_atom` so that :func:`_contacts` can
    call it without re-entering itself: the fallback resolves a cell by
    measuring its own neighbourhood, and measuring a neighbourhood asks for
    contacts. Every call here either answers or descends to a strictly coarser
    cell, so the recursion terminates at resolution 1.
    """
    if resolution > 1:
        return _descend_neighbor(cell, direction)
    target, target_x, target_y = _lexical_target(cell, direction)
    if target_x >= 0 and target_y >= 0 and is_valid_cell(target):
        return target
    return _deflect_lattice(cell, direction)


def _changes_frame(cell: str, answer: str) -> bool:
    """Whether a step crossed into another frame rather than moving inside one.

    Only a frame change can make the arithmetic non-invertible, so only a
    frame change is worth the two rings it costs to check. Every boundary the
    normalisation handles changes quadrant — the meridian, the equator, the
    seam, and the overflow into a shorter row, which the seam branch resolves
    — except the step onto a polar cell, which keeps its quadrant and is
    therefore named separately.

    Checking instead whether the step left its resolution-1 cell would be
    true of every step at resolution 1, and would put a geometric test on the
    path of the interior, which is the path the phase exists to keep lexical.
    """
    origin, target = split_components(cell), split_components(answer)
    if origin[0] != target[0]:
        return True
    return int(target[1].partition("/")[2]) == _POLAR_ROW


def _require_contact(cell: str, answer: str, direction: Direction) -> None:
    """Refuse an answer that leaves the base cell without touching it.

    A step that changes frame is not always invertible by arithmetic. The
    column at the prime meridian holds a single triangle spanning both
    quadrants, so read from either side it is two columns wide; a refinement
    of it composed against ordinary parallelograms does not lay out as a grid,
    and a diagonal crossing it is not a lattice direction at all. Rather than
    special-case the triangle, the candidate is checked, and one that turns
    out to touch nothing hands the cell over to measurement — which reports
    the contact even where no lattice step names it.

    Checked only when the step changed frame. Steps that stay in one frame
    are pure arithmetic and pay nothing for this.
    """
    wanted = "edge" if direction in _CARDINALS else "vertex"
    if _contact(_ring(cell), _ring(answer)) != wanted:
        raise DomainError(
            f"stepping {direction} from {cell!r} crosses a frame the index "
            "cannot invert; the neighbour must be measured"
        )


def _by_contact(cell: str, direction: Direction) -> str | None:
    """The neighbour in a direction, taken from the measured contact sets.

    Where a side carries several cells — the base of a meridian triangle is
    split between the quadrants, a trapezoid can face three — the one sharing
    the longest boundary is returned. The tie-break is stated rather than
    hidden: :func:`cell_to_edges` returns every one of them.
    """
    edges, vertices = _contacts(cell)
    cardinal = direction in _CARDINALS
    pool = edges if cardinal else vertices
    allowed = _CARDINALS if cardinal else _DIAGONALS
    matching = [other for other in pool if _classify(cell, other, allowed) == direction]
    if not matching:
        return None
    return max(matching, key=lambda other: _shared_length(cell, other))


def _has_exceptional_ancestor(cell: str) -> bool:
    """Whether any ancestor of a refined cell is itself exceptional.

    A child of a triangle or of an absorbing trapezoid does not tile a square
    refinement grid, so the column-and-row arithmetic would name the wrong
    sibling rather than fail. Checking the chain is cheap: it is at most
    thirteen prefixes.
    """
    components = split_components(cell)
    for depth in range(2, len(components)):
        if _is_exceptional(join_components(components[:depth])):
            return True
    return False


def _grid_side(resolution: int) -> int:
    """Refinement grid width at a resolution: 2 for even, 5 for odd."""
    return _QUATERNARY_SIDE if resolution % 2 == 0 else _QUINARY_SIDE


def _decode(code: str, side: int) -> tuple[int, int]:
    """Refinement code to ``(column, row)``, both zero-based, row poleward.

    Quaternary codes run ``1`` south-west, ``2`` south-east, ``3`` north-west,
    ``4`` north-east. Quinary codes name the row with the letter, ``A``
    equatorward through ``E`` poleward, and the column with the digit, ``1``
    nearest the prime meridian.
    """
    if side == _QUATERNARY_SIDE:
        ordinal = int(code) - 1
        return ordinal % side, ordinal // side
    return int(code[1]) - 1, _ALPHABET.index(code[0])


def _encode(column: int, row: int, side: int) -> str:
    """Inverse of :func:`_decode`."""
    if side == _QUATERNARY_SIDE:
        return str(row * side + column + 1)
    return f"{_ALPHABET[row]}{column + 1}"


def _refinable(cell: str) -> bool:
    """Whether a cell's children tile it as a plain square refinement grid.

    Triangular and border-absorbing cells do not: the meridian triangle packs
    its twenty-five codes into rows of nine, seven, five, three and one, and a
    trapezoid keeps only the codes its clipped area still reaches, so neither
    admits the column-and-row arithmetic used below.
    """
    return cell_shape(cell) == "parallelogram" and not absorbs_border(cell)


def _descend_neighbor(cell: str, direction: Direction) -> str | None:
    """Neighbour of a refined cell: step in the parent, or ascend and descend.

    The refinement grid is square and the step stays inside it whenever the
    result does. When it leaves, the parent takes the same step — recursively,
    so a run of corner cells walks up as far as it must — and the coordinate
    that overflowed wraps to the far side of the neighbouring parent.
    """
    components = split_components(cell)
    resolution = len(components) - 1
    side = _grid_side(resolution)
    parent = _parent_cell(cell)
    # The parent is refinable by construction: an unrefinable one is an
    # exceptional ancestor, and those are routed to measurement before the
    # descent is entered at all.
    column, row = _decode(components[-1], side)
    shift_column, shift_row = GRID_STEP[direction]
    column, row = column + shift_column, row + shift_row

    overflow = (
        (column // side if column < 0 or column >= side else 0),
        (row // side if row < 0 or row >= side else 0),
    )
    outward = (
        max(-1, min(1, overflow[0])),
        max(-1, min(1, overflow[1])),
    )
    if outward == (0, 0):
        return join_components(split_components(parent) + [_encode(column, row, side)])

    neighbour_parent = _neighbor_of_atom(parent, _COMPASS[outward])
    if neighbour_parent is None:
        return None
    if not _refinable(neighbour_parent):
        raise DomainError(
            f"{neighbour_parent!r} is a {cell_shape(neighbour_parent)}; the "
            "refinement grid of a triangle or an absorbing trapezoid is not "
            "square, so the wrapped code would not name the right child"
        )
    code = _encode(column % side, row % side, side)
    return join_components(split_components(neighbour_parent) + [code])


# --------------------------------------------------------------------------
# Boundary normalisation
# --------------------------------------------------------------------------


def _normalize_target(
    quadrant: str, column: int, row: int, shift_row: int
) -> tuple[str, int, int] | None:
    """Bring a lexical target back into range by changing frame, not by cases.

    The lexical step is allowed to spell a target outside its quadrant's range
    in either component or in both. This maps such a target onto the cell it
    actually denotes, by applying the frame transformation of each boundary it
    crossed. Composing them is what lets a diagonal cross two boundaries at
    once, which the earlier one-dimensional dispatch could not express and
    answered ``None`` to.

    The transformations, each measured rather than assumed:

    - **Equator**, ``row < 0``. The hemispheres meet column to column, so the
      row reflects to ``-1 - row``, and the column takes back the vertical
      component of the step, because the shear mirrors: poleward moves one
      column toward the meridian in the north and one away from it in the
      south.
    - **Polar row**, ``row > 999``. One cell per hemisphere, in the eastern
      quadrant, spanning the whole parallel.
    - **Prime meridian**, ``column < 0`` or ``column == 0`` in a western
      quadrant. The triangle at column zero belongs to the eastern quadrant
      alone, so a step into or across it changes quadrant.
    - **Antemeridian seam**, ``column`` past the row's last. The two quadrants
      meet at the last column of each, whatever longitude that happens to be.

    Returns ``None`` when the target lies outside the addressable domain.
    """
    if row < 0:
        quadrant = _ACROSS_EQUATOR[quadrant]
        column, row = column + shift_row, -1 - row
        shift_row = -shift_row

    if row > _POLAR_ROW:
        return None
    if row == _POLAR_ROW:
        return (_eastern(quadrant), 0, _POLAR_ROW)

    if column < 0:
        quadrant = _ACROSS_MERIDIAN[quadrant]
        column = -column
    elif column == 0 and quadrant[1] == "W":
        quadrant = _ACROSS_MERIDIAN[quadrant]

    last = _last_column(quadrant, row)
    if column > last:
        mirror = _ACROSS_MERIDIAN[quadrant]
        overshoot = column - last - 1
        mirrored = _last_column(mirror, row) - overshoot
        if mirrored < 0:
            return None
        return (mirror, mirrored, row)
    return (quadrant, column, row)


def _eastern(quadrant: str) -> str:
    """The eastern quadrant of a quadrant's hemisphere."""
    return quadrant if quadrant[1] == "E" else _ACROSS_MERIDIAN[quadrant]


def deflect(cell: str, direction: Direction) -> str | None:
    """Resolve a step whose lexical target falls outside its quadrant's range.

    Takes a true compass direction, like :func:`get_neighbor`. The
    lattice-frame worker behind it is :func:`_deflect_lattice`.
    """
    return _deflect_lattice(cell, _lattice_direction(cell[:2], direction))


def _deflect_lattice(cell: str, direction: Direction) -> str | None:
    """Resolve a step whose lexical target falls outside its quadrant's range.

    Isolated from :func:`get_neighbor` because it is the subtle part. It does
    not decide between boundary cases; it normalises the target through
    whatever frame changes it crossed, which is the same operation for one
    boundary or two and therefore works for the diagonals as well.

    Args:
        cell: Atomic resolution-1 index at a boundary.
        direction: Lattice direction of the step.

    Returns:
        The deflected neighbour index, or ``None`` when the step leaves the
        addressable domain or when the normalised target names no cell.
    """
    quadrant, column, row = _res1_parts(cell)
    shift_column, shift_row = LATTICE_STEP[direction]
    normalised = _normalize_target(
        quadrant, column + shift_column, row + shift_row, shift_row
    )
    if normalised is None:
        return None
    candidate = _spell(*normalised)
    return candidate if is_valid_cell(candidate) else None


# --------------------------------------------------------------------------
# The geometric exception set
# --------------------------------------------------------------------------


def _is_exceptional(cell: str) -> bool:
    """Whether a cell's adjacency is beyond what the lexical step can name.

    Four families, all recognised without looking at a coordinate:

    - cells that absorb the border, whose sides carry more than one neighbour;
    - triangles, whose base is divided between the two quadrants and whose
      apex meets the cell above at a point;
    - the two polar cells, which cover a whole parallel;
    - ordinary parallelograms sitting on a row where an extension zone starts
      or ends, where one quadrant is extended past 180 degrees while its
      mirror is truncated and the neighbour opposite is in the other quadrant.

    The set is about six thousand cells at resolution 1, out of six million.
    Everything outside it is answered lexically, which is what keeps criterion
    9 meaningful rather than vacuous.
    """
    components = split_components(cell)
    if len(components) < 2:
        return True
    row = int(components[1].partition("/")[2])
    if row in ZONE_LIMIT_ROWS or row >= _POLAR_ROW:
        return True
    return bool(absorbs_border(cell)) or cell_shape(cell) != "parallelogram"


def _ring(cell: str, longitude_offset: float = 0.0) -> Polygon:
    ring = Polygon(cell_to_boundary(cell))
    return translate(ring, xoff=longitude_offset) if longitude_offset else ring


def _contact(first: Polygon, second: Polygon) -> str:
    """Whether two rings share an edge, a single point, or nothing.

    The antemeridian is written as plus or minus 180 but stored a unit in the
    last place either side of it, so an exact intersection reports a gap of
    about six nanometres on roughly a third of the rows. Snapping first makes
    this a statement about the tessellation rather than about binary64.
    """
    for offset in (0.0, 360.0, -360.0):
        shifted = translate(second, xoff=offset) if offset else second
        if first.distance(shifted) > 1.0:
            continue
        overlap = first.intersection(snap(shifted, first, _SNAP_TOLERANCE))
        if overlap.is_empty:
            continue
        if overlap.length > _EDGE_TOLERANCE:
            return "edge"
        return "vertex"
    return "none"


def _geometric_candidates(cell: str) -> list[str]:
    """Cells that could possibly touch ``cell``, as a bounded, complete window.

    Three column bands per neighbouring row, because a cell in the exception
    set can have neighbours in any of them: around its own column, at the
    prime meridian, and at the antemeridian seam. The seam band is what a
    window centred on the column alone would miss, and missing it is silent —
    the cell simply loses a neighbour.

    Three quadrants are searched for most cells: the cell's own, its mirror
    across the meridian, and its mirror across the equator. A fourth, the
    quadrant diagonally opposite, is searched when the cell's own ring reaches
    both axes, because then a single vertex of it lies on the point the four
    quadrants share and the cells meeting there include one on the far side of
    both.

    That fourth quadrant used to be excluded by an assertion in this
    docstring, which said no cell touches it. A brute-force scan of every
    resolution-1 cell within three columns and three rows of the origin, in
    all four quadrants, falsified it: ``SE(0000/0000)`` and
    ``NW(0001/0000)`` meet at longitude -0.08983152841195215, latitude 0,
    which is a vertex both rings carry with identical bits. The scan is
    ``test_the_corner_triangle_touches_the_opposite_quadrant``.

    The reach is decided by measuring the ring rather than by naming the
    corner cell. A cell that reaches neither axis, or only one, pays nothing:
    the fourth quadrant is not searched for it.
    """
    quadrant, _, row = _res1_parts(cell)
    west, south, east, north = _ring(cell).bounds
    quadrants = [quadrant, _ACROSS_MERIDIAN[quadrant], _ACROSS_EQUATOR[quadrant]]
    if west <= 0.0 <= east and south <= 0.0 <= north:
        quadrants.append(_ACROSS_MERIDIAN[_ACROSS_EQUATOR[quadrant]])
    # Two readings of the cell's longitude span, kept apart. An extension
    # cell's ring runs past 180 degrees, and the cells below it on an ordinary
    # row are addressed from the other side of the seam, at 360 minus that.
    # Merging the two into one interval would sweep everything between them —
    # more than twenty degrees of longitude for Chukotka, which is hundreds of
    # columns of candidates that cannot touch anything.
    direct = (min(abs(west), abs(east)), max(abs(west), abs(east)))
    mirrored = (
        min(abs(360.0 - abs(west)), abs(360.0 - abs(east))),
        max(abs(360.0 - abs(west)), abs(360.0 - abs(east))),
    )
    candidates: list[str] = []
    for other in quadrants:
        for other_row in (row - 1, row, row + 1):
            if not 0 <= other_row <= _POLAR_ROW:
                continue
            last = _last_column(other, other_row)
            step = _degrees_per_column(other, other_row)
            columns = {
                *range(0, min(last, _WINDOW) + 1),
                *range(max(0, last - _SEAM_REACH), last + 1),
            }
            for low, high in (direct, mirrored):
                columns |= set(
                    range(
                        max(0, int(low / step) - _WINDOW),
                        min(last, int(high / step) + _WINDOW) + 1,
                    )
                )
            for other_column in sorted(columns):
                candidate = _spell(other, other_column, other_row)
                if candidate != cell and is_valid_cell(candidate):
                    candidates.append(candidate)
    return candidates


@lru_cache(maxsize=8192)
def _degrees_per_column(quadrant: str, row: int) -> float:
    """Longitude a column spans in a row, asked of the row rather than assumed.

    Not ``180 / (last + 1)``. Inside an extension zone one quadrant runs past
    180 degrees and its mirror stops short of it, so that formula misplaces a
    column by tens of degrees exactly where the neighbours are hardest to
    find. Reading the far edge of the row's last cell covers both cases.
    """
    last = _last_column(quadrant, row)
    bounds = _ring(_spell(quadrant, last, row)).bounds
    return float(max(abs(bounds[0]), abs(bounds[2]))) / (last + 1)


def _side_labels(cell: str) -> list[tuple[LineString, Direction]]:
    """Label each side of a cell's ring with the lattice direction it faces.

    Done on the cell's own ring, in its own unwrapped longitudes, so no
    coordinate is ever compared across a quadrant. That is what makes this
    survive the antemeridian: an extension cell's ring runs past 180 degrees
    and its outward side is simply the one further from the prime meridian,
    with nothing to unwrap.

    A side lying along a parallel faces poleward or equatorward; any other
    side faces away from the prime meridian or toward it. Both are decided
    against the cell's own centroid and signed by the quadrant, so the labels
    come out as lattice directions and read the same in all four quadrants.

    A triangle has three sides and is therefore missing one label. The polar
    cell's two non-base sides are the same meridian glued to itself and carry
    no neighbour at all.
    """
    quadrant = split_components(cell)[0]
    east_sign = 1.0 if quadrant[1] == "E" else -1.0
    pole_sign = 1.0 if quadrant[0] == "N" else -1.0
    ring = cast(list[tuple[float, float]], cell_to_boundary(cell))
    centre = Polygon(ring).centroid

    labels: list[tuple[LineString, Direction]] = []
    for index, start_point in enumerate(ring):
        end_point = ring[(index + 1) % len(ring)]
        middle = (
            (start_point[0] + end_point[0]) / 2.0,
            (start_point[1] + end_point[1]) / 2.0,
        )
        if start_point[1] == end_point[1]:
            outward = (middle[1] - centre.y) * pole_sign > 0.0
            direction = "N" if outward else "S"
        else:
            outward = (middle[0] - centre.x) * east_sign > 0.0
            direction = "E" if outward else "W"
        labels.append(
            (LineString([start_point, end_point]), cast(Direction, direction))
        )
    return labels


def _shared_geometry(first: str, second: str) -> BaseGeometry | None:
    """The boundary two cells hold in common, with the seam shift applied."""
    here = _ring(first)
    best: BaseGeometry | None = None
    best_length = -1.0
    for offset in (0.0, 360.0, -360.0):
        there = _ring(second, offset)
        if here.distance(there) > 1.0:
            continue
        overlap = here.intersection(snap(there, here, _SNAP_TOLERANCE))
        if not overlap.is_empty and overlap.length > best_length:
            best, best_length = overlap, overlap.length
    return best


def _classify(
    origin: str, other: str, allowed: tuple[Direction, ...]
) -> Direction | None:
    """Which of ``allowed`` directions ``other`` lies in from ``origin``.

    Found by asking which side of the origin's own ring carries the shared
    boundary, not by comparing positions. A trapezoid with three cells along
    one side gives all three the same label, which is correct and is why
    :func:`cell_to_edges` returns the set while :func:`get_neighbor` returns
    one of them.

    A vertex contact sits where two sides meet, so it is named by the pair,
    and that pair is the diagonal.
    """
    shared = _shared_geometry(origin, other)
    if shared is None:
        return None
    # Against the shared boundary's centre, not against the boundary itself: a
    # shared edge reaches the ring's corners, so it lies at distance zero from
    # the two sides meeting there as well as from its own.
    centre = shared.centroid
    sides = sorted(_side_labels(origin), key=lambda item: item[0].distance(centre))
    if allowed == _CARDINALS:
        return sides[0][1]
    pair = {sides[0][1], sides[1][1]}
    for direction in allowed:
        if set(direction) == pair:
            return direction
    return None


def _shared_length(first: str, second: str) -> float:
    """Length of the boundary two cells hold in common, zero when none."""
    here = _ring(first)
    for offset in (0.0, 360.0, -360.0):
        there = _ring(second, offset)
        if here.distance(there) > 1.0:
            continue
        overlap = here.intersection(snap(there, here, _SNAP_TOLERANCE))
        if overlap.length > _EDGE_TOLERANCE:
            return float(overlap.length)
    return 0.0


def _refined_candidates(cell: str) -> list[str]:
    """Cells at ``cell``'s resolution that could touch it, found by descent.

    Anything adjacent to a refined cell is a child of its parent or a child of
    one of the parent's own neighbours, edge or vertex. Resolving the parent
    first and descending is therefore complete, and it costs one level of
    recursion rather than a scan of the row.
    """
    parent = _parent_cell(cell)
    families = [parent, *_edge_neighbors(parent), *_vertex_neighbors(parent)]
    candidates: list[str] = []
    for family in families:
        for child in get_children(family, flatten=True):
            if child != cell:
                candidates.append(cast(str, child))
    return candidates


def _neighbors_by_contact(cell: str) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Edge neighbours and vertex neighbours of a cell, measured.

    The only place in the module that looks at a coordinate to *find* a
    neighbour rather than to check one. Reached only inside the exception set,
    where the lexical step provably cannot name the answer.
    """
    candidates = (
        _geometric_candidates(cell)
        if get_resolution(cell) == 1
        else _refined_candidates(cell)
    )
    here = _ring(cell)
    edges: list[str] = []
    vertices: list[str] = []
    for candidate in candidates:
        contact = _contact(here, _ring(candidate))
        if contact == "edge":
            edges.append(candidate)
        elif contact == "vertex":
            vertices.append(candidate)
    return tuple(sorted(set(edges))), tuple(sorted(set(vertices)))


@lru_cache(maxsize=4096)
def _contacts(cell: str) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Cached neighbour sets. Lexical where possible, measured where not."""
    if _is_exceptional(cell):
        return _neighbors_by_contact(cell)
    resolution = get_resolution(cell)
    if resolution > 1 and _has_exceptional_ancestor(cell):
        return _neighbors_by_contact(cell)
    try:
        found = {
            direction: _verified_lexical(cell, direction, resolution)
            for direction in _CARDINALS + _DIAGONALS
        }
    except DomainError:
        # The index stopped being able to represent the step. Measure instead.
        return _neighbors_by_contact(cell)
    edges = tuple(sorted(n for d in _CARDINALS if (n := found[d]) is not None))
    vertices = tuple(sorted(n for d in _DIAGONALS if (n := found[d]) is not None))
    return edges, vertices


def _edge_neighbors(cell: str) -> tuple[str, ...]:
    """Every cell sharing an edge with this one."""
    return _contacts(cell)[0]


def _vertex_neighbors(cell: str) -> tuple[str, ...]:
    """Every cell meeting this one at a single vertex."""
    return _contacts(cell)[1]


def _touching(cell: str, metric: Metric) -> tuple[str, ...]:
    """Neighbours a metric may step to: edges only, or edges and vertices."""
    edges, vertices = _contacts(cell)
    return edges if metric == "manhattan" else edges + vertices


def _atoms(index: str) -> list[str]:
    """Cells of an index, in the order :func:`itacart.index.decompose` gives."""
    return decompose(index)


def _directions(metric: Metric) -> tuple[Direction, ...]:
    """Directions a metric may step in.

    Chebyshev admits the diagonals, which meet the origin at a vertex, so its
    unit disk is the nine cells around a lattice point. Manhattan admits only
    the edge steps, so its unit disk is five.
    """
    if metric == "chebyshev":
        return _CARDINALS + _DIAGONALS
    if metric == "manhattan":
        return _CARDINALS
    raise ValueError(
        f"{metric!r} is not a lattice metric; expected chebyshev or manhattan"
    )


def _expand(origin: str, k_distance: int, metric: Metric) -> list[set[str]]:
    """Shells around one origin: index ``n`` holds the cells at exactly ``n`` steps.

    Breadth-first over :func:`get_neighbor`, so the shape of the boundary is
    whatever the deflections say it is rather than something this function
    assumes. Composing steps is also what keeps the result symmetric: a shell
    is a set, and set membership cannot disagree with itself.
    """
    if k_distance < 0:
        raise ValueError(f"k_distance must not be negative, got {k_distance}")
    _directions(metric)
    shells: list[set[str]] = [{origin}]
    seen = {origin}
    for _ in range(k_distance):
        frontier: set[str] = set()
        for cell in shells[-1]:
            frontier |= {n for n in _touching(cell, metric) if n not in seen}
        seen |= frontier
        shells.append(frontier)
    return shells


def _shape_result(
    per_origin: list[list[str]], dedupe: bool, flatten: bool, single: bool
) -> list[str] | list[list[str]]:
    """Apply the dedupe and flatten options to one result per origin."""
    if dedupe:
        counts: dict[str, int] = {}
        for group in per_origin:
            for cell in group:
                counts[cell] = counts.get(cell, 0) + 1
        per_origin = [[c for c in group if counts[c] == 1] for group in per_origin]
    if single:
        return per_origin[0]
    if flatten:
        return sorted({cell for group in per_origin for cell in group})
    return per_origin


def grid_disk(
    index: str,
    k_distance: int = 1,
    metric: Metric = "chebyshev",
    dedupe: bool = False,
    flatten: bool = False,
) -> list[str] | list[list[str]]:
    """Cells within ``k`` grid steps of the origin, the origin included.

    ``chebyshev`` admits diagonal steps, giving a filled rhombus on the
    lattice; ``manhattan`` admits only axis steps, giving a filled diamond. In
    the interior the unit disks hold nine and five cells respectively.

    Args:
        index: Compositional index string.
        k_distance: Ring radius in grid steps.
        metric: Lattice metric to expand under.
        dedupe: When the input holds several cells, drop cells that appear in
            more than one disk. Breaks positional alignment, so it is off by
            default.
        flatten: When the input holds several cells, ``False`` returns one
            list per input cell; ``True`` returns a single flat list.

    Returns:
        A cell list for a single origin, or a list of lists aligned with
        :func:`itacart.index.decompose`.

    Raises:
        ValueError: If ``k_distance`` is negative or ``metric`` is unknown.
        DomainError: If the disk reaches the geometric exception set, where a
            direction is multivalued and composing steps would lose cells
            silently.
    """
    atoms = _atoms(index)
    per_origin = [
        sorted(set().union(*_expand(atom, k_distance, metric))) for atom in atoms
    ]
    return _shape_result(per_origin, dedupe, flatten, len(atoms) == 1)


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
        A cell list for a single origin, or a positionally aligned list of
        lists. A ring of radius zero is the origin alone.

    Raises:
        ValueError: If ``k_distance`` is negative or ``metric`` is unknown.
        DomainError: If the ring reaches the geometric exception set.
    """
    atoms = _atoms(index)
    per_origin = [sorted(_expand(atom, k_distance, metric)[-1]) for atom in atoms]
    return _shape_result(per_origin, dedupe, flatten, len(atoms) == 1)


def _lattice_descent(cell: str) -> tuple[int, int, bool]:
    """A cell's place on the lattice of its own resolution, by descent.

    The reference is the resolution-1 cell and every level below it only
    subdivides the unit it inherits. Three things vary with the parent and
    all three are read from the descent itself rather than from the index
    prefix, which names a quadrant and not a side of the line.

    The shear. A poleward step moves one column towards the meridian at
    every scale, so a child's column is its parent's scaled column plus
    ``code_column - code_row``. Folding the code in without that term
    treats the refinement grid as unsheared and measurably misplaces every
    child on an odd row.

    The orientation. West of the meridian the lattice mirrors and the same
    term changes sign. A cell of an eastern quadrant can sit wholly west --
    ``NE(0000/0110(3))`` runs from -0.09119 degrees to zero -- so the sign
    comes from the running column, not the quadrant.

    The shape. A triangle's twenty-five children pack into rows of nine,
    seven, five, three and one rather than a five-by-five grid. Measured,
    the row of a child is ``min(code_row, code_column)`` instead of
    ``code_row``, the column term is unchanged, and the one child that is
    itself a triangle is the one whose code has ``code_row`` equal to
    ``code_column``. That makes the shape arithmetic too: a triangle at
    resolution 1 is column zero, and below it the diagonal code.
    """
    components = split_components(cell)
    quadrant = components[0]
    column, _, row = components[1].partition("/")
    x, y = int(column), int(row)
    if quadrant[1] == "W":
        x = -x
    triangular = x == 0
    for depth, code in enumerate(components[2:], start=2):
        side = _grid_side(depth)
        code_column, code_row = _decode(code, side)
        towards = -1 if x < 0 else 1
        y = y * side + (min(code_row, code_column) if triangular else code_row)
        x = x * side + towards * (code_column - code_row)
        triangular = triangular and code_row == code_column
    return x, y, quadrant[0] == "S"


def _rectified(cell: str) -> tuple[int, int]:
    """A cell's position in the frame where the lattice is a square grid.

    The lattice is spanned by the eastward step ``(1, 0)`` and the poleward
    step ``(-1, 1)``, and the unimodular map ``(X, Y) -> (X + Y, Y)`` sends
    those to ``(1, 0)`` and ``(0, 1)``. Integer in and integer out, with an
    integer inverse: no scale is lost and no float is introduced.

    Two departures from that plain statement, both measured. The eastward
    coordinate uses the distance from the meridian rather than a signed
    column, because west of the line the shear mirrors and ``X + Y`` would
    read it backwards. And the southern rows continue the northern ones
    downwards through ``-Y - 1`` rather than through ``-Y``: the equator is
    a mirror lying between two rows, not through one, so there is no shared
    row to count twice. Signing the row without the offset leaves an error
    that grows with the distance from the equator instead of a fixed one.
    """
    x, y, southern = _lattice_descent(cell)
    return abs(x) + y, (-y - 1) if southern else y


def _span(here: tuple[int, int], there: tuple[int, int], metric: Metric) -> int:
    """The metric, read in the rectified frame."""
    east, pole = abs(there[0] - here[0]), abs(there[1] - here[1])
    return max(east, pole) if metric == "chebyshev" else east + pole


def _seam_row_limit(resolution: int) -> int:
    """Rows between the equator and the pole at a resolution."""
    rows = _POLAR_ROW
    for depth in range(2, resolution + 1):
        rows *= _grid_side(depth)
    return rows


def _is_polar_cap(cell: str) -> bool:
    """Whether a cell is the single triangle that closes a hemisphere.

    Resolution 1 only, and the restriction is the measurement rather than a
    simplification. What makes the cap worth a case of its own is that its
    other two contacts are trapezoids of the last column, leaving one way
    out; refine it once and that stops being true, with three contacts free
    of trapezoids at resolution 2 and eight at resolution 3. The reduction
    below is a cut of the graph at a single vertex, and there is no single
    vertex to cut below the top level.
    """
    _, row, _ = _lattice_descent(cell)
    return get_resolution(cell) == 1 and row >= _seam_row_limit(1)


def _polar_equatorward_triangle(cell: str) -> str:
    """The one cell a polar cap can be left by.

    The cap touches three cells and two of them are trapezoids of the last
    column, which this measure does not step on. What remains is the meridian
    triangle of the row below, and because it is the only way out, the
    distance from the cap is one step more than the distance from it -- a cut
    of the graph at a single vertex rather than a shortcut chosen among
    several.
    """
    components = split_components(cell)
    quadrant = components[0]
    equatorward = _seam_row_limit(get_resolution(cell)) - 1
    descent = [quadrant, f"0000/{equatorward:04d}"]
    return join_components(descent)


def _antimeridian_reach(cell: str) -> int:
    """Columns between a cell and the last addressable one of its row.

    Rows narrow towards the pole, so this is read from the row the cell
    actually sits on rather than from a global width. Measured, the
    descended column never exceeds the limit the row declares, at every
    resolution checked.
    """
    column, row, _ = _lattice_descent(cell)
    quadrant = split_components(cell)[0]
    resolution = get_resolution(cell)
    return last_lattice_column(quadrant, row, cell_size(resolution)) - abs(column)


def _crosses_antimeridian(origin: str, destination: str, through_meridian: int) -> bool:
    """Whether the shorter way between two cells leaves through the seam.

    The measure knows one seam, the prime meridian, and reads a distance as
    the path that reaches it and crosses. Two cells far to the east and west
    are also neighbours the other way round, across the antimeridian, and
    the contact sets report that contact: measured, cells two columns inside
    the last one face each other three steps apart while the prime-meridian
    route runs to 1547.

    That other route is not a path this measure can express. ITACaRT
    addresses the land beyond the line through an extension zone rather than
    by stepping across it, so a number returned here would name a walk the
    index does not admit. The pair is refused instead.

    The test is a comparison and not a proximity, which is what keeps it
    from refusing pairs it should answer: the seam route costs each cell its
    reach to the last column of its own row plus the step across, and only
    when that total does not exceed the prime-meridian reading is the
    reading the wrong one to return.
    """
    through_seam = _antimeridian_reach(origin) + _antimeridian_reach(destination) + 1
    return through_seam <= through_meridian


def _seam(rectified_row: int) -> tuple[int, int]:
    """The meridian triangle of a rectified row, in the rectified frame.

    Column zero carries one triangle per row and it is the cell both sides
    of the line share. Its rectified position is ``(row, row)`` north of
    the equator and ``(-row - 1, row)`` south of it, which is the same
    statement once the southern rows are counted downwards.
    """
    return (rectified_row if rectified_row >= 0 else -rectified_row - 1), rectified_row


def _apex(rectified_row: int) -> tuple[int, int]:
    """The pair of cells that meet at a meridian triangle's apex.

    They face each other across the line one column out from the seam, and
    the rectified frame gives them the same coordinates because it measures
    the distance from the meridian rather than a signed column. Enumerated
    rather than assumed: every contact that crosses the meridian without
    touching a triangle has the eastern cell at column ``+1``, the western
    at ``-1`` and no row difference, in both hemispheres and at every
    resolution measured.
    """
    column, row = _seam(rectified_row)
    return column + 1, row


def _rows_to_search(
    here: tuple[int, int], there: tuple[int, int], resolution: int
) -> list[tuple[int, int]]:
    """Bounds of the crossing row, one interval per hemisphere.

    The seam bends at the equator: the column of a seam cell rises with the
    row to the north and falls with it to the south. Each branch on its own
    is affine, which is what the search below needs, so the two are searched
    separately and the better one wins.

    Each interval is clipped to the rows that exist and widened by two
    around the coordinates that can hold a breakpoint. The minimum of a sum
    of distances sits at a breakpoint, and every breakpoint is one of the
    four coordinates.
    """
    limit = _seam_row_limit(resolution)
    reach = [here[0], here[1], there[0], there[1]]
    low, high = min(reach) - 2, max(reach) + 2
    intervals = []
    for first, last in ((0, limit - 1), (-limit, -1)):
        lower, upper = max(first, low), min(last, high)
        if lower > upper:
            lower = upper = first if low > last else last
        intervals.append((lower, upper))
    return intervals


def _minimise(cost: Callable[[int], int], lower: int, upper: int) -> int:
    """Least value of a discretely convex cost over an integer interval.

    Each cost below is a sum of maxima of absolute affine functions of the
    row, so it is convex, its breakpoints are integers, and a ternary search
    finds the minimum in a logarithmic number of evaluations. Measured
    against exhaustive enumeration over the same intervals, the two agree
    everywhere.
    """
    while upper - lower > 2:
        left = lower + (upper - lower) // 3
        right = upper - (upper - lower) // 3
        if cost(left) <= cost(right):
            upper = right
        else:
            lower = left
    return min(cost(row) for row in range(lower, upper + 1))


def _meridian_distance(
    origin: str, destination: str, metric: Metric, resolution: int
) -> int:
    """Steps between cells on opposite sides of the meridian.

    Two structures cross the line and both are needed. One goes through the
    shared triangle of column zero; the other steps directly between the
    cells that meet at its apex. Measured on every pair of an exhaustive
    neighbourhood at three resolutions and in both hemispheres, neither
    family alone reproduces the shortest path and their minimum always
    does: pairs exist that only the triangle reaches, and pairs exist where
    the triangle is unavoidable, so the apex cannot stand alone either.

    Under Manhattan the apex is not available at all. The two cells there
    share a single point, and a vertex is not a Manhattan step, so charging
    one for it would report a path that cannot be walked.
    """
    here, there = _rectified(origin), _rectified(destination)

    def through_seam(row: int) -> int:
        point = _seam(row)
        return _span(here, point, metric) + _span(point, there, metric)

    def through_apex(row: int) -> int:
        point = _apex(row)
        return _span(here, point, metric) + 1 + _span(point, there, metric)

    best = None
    for lower, upper in _rows_to_search(here, there, resolution):
        candidates = [_minimise(through_seam, lower, upper)]
        if metric == "chebyshev":
            candidates.append(_minimise(through_apex, lower, upper))
        found = min(candidates)
        best = found if best is None else min(best, found)
    assert best is not None
    return best


def grid_distance(origin: str, destination: str, metric: Metric = "chebyshev") -> int:
    """Steps between two cells on the lattice.

    Not a geodesic distance. Use :func:`itacart.geodesy.inverse_geodesic` on
    the anchors or centroids for metric separation.

    Defined at every resolution. The lattice mirrors at the prime meridian
    and again at the equator, so no single basis measures both sides of
    either: the position is built by descent from the resolution-1 cell and
    read in a rectified frame, and the meridian crossing is resolved by
    minimising over the row at which the line is crossed.

    Which side of the meridian a cell lies on comes from that descent and
    never from its quadrant prefix. A cell of an eastern quadrant can sit
    wholly west of the line, and two cells sharing a prefix can face each
    other across it.

    Args:
        origin: Atomic index of the first cell.
        destination: Atomic index of the second cell.
        metric: Lattice metric to measure under.

    Returns:
        Step count, zero when the cells coincide.

    Raises:
        ResolutionError: If the cells sit at different resolutions, or are
            quadrants.
        DomainError: If a path cannot be resolved across the boundary
            between the cells, which includes either cell being a trapezoid
            of the last addressable column.
        ValueError: If ``metric`` is unknown.
    """
    if metric not in ("chebyshev", "manhattan"):
        raise ValueError(
            f"{metric!r} is not a lattice metric; expected chebyshev or manhattan"
        )
    if get_resolution(origin) != get_resolution(destination):
        raise ResolutionError(
            f"{origin!r} and {destination!r} sit at different resolutions; grid "
            "distance is defined between cells of one resolution"
        )
    if get_resolution(origin) < 1:
        raise ResolutionError("a quadrant has no position on the lattice")
    for cell in (origin, destination):
        if cell_shape(cell) == "trapezoid":
            raise DomainError(
                f"{cell!r} is a trapezoid; the lattice step is not defined at "
                "the last addressable column of a row, where the row above "
                "narrows and the poleward neighbour lies several columns in. "
                "Land beyond the line is addressed through an extension zone "
                "rather than by stepping across it"
            )
    if origin == destination:
        return 0
    for cell in (origin, destination):
        if _is_polar_cap(cell) and metric == "manhattan":
            raise DomainError(
                f"{cell!r} closes its hemisphere and its one remaining "
                "contact is a shared vertex; a vertex is not a Manhattan "
                "step, so no Manhattan path leaves it"
            )
    if _is_polar_cap(origin):
        return 1 + grid_distance(
            _polar_equatorward_triangle(origin), destination, metric
        )
    if _is_polar_cap(destination):
        return 1 + grid_distance(
            origin, _polar_equatorward_triangle(destination), metric
        )
    resolution = get_resolution(origin)
    here, there = _lattice_descent(origin), _lattice_descent(destination)
    if (here[0] > 0) == (there[0] > 0) or here[0] == 0 or there[0] == 0:
        return _span(_rectified(origin), _rectified(destination), metric)
    crossing = _meridian_distance(origin, destination, metric, resolution)
    if _crosses_antimeridian(origin, destination, crossing):
        raise DomainError(
            "grid_distance is not defined across the antimeridian seam; "
            f"{origin!r} and {destination!r} are no further apart around the "
            "far side of the lattice than through the prime meridian, and "
            "ITACaRT reaches the land beyond the line through an extension "
            "zone rather than by stepping across it"
        )
    return crossing


def are_neighbor_cells(origin: str, destination: str) -> bool:
    """Whether two cells share an edge.

    Edge adjacency only; cells meeting at a single vertex are not neighbours,
    which is what separates criterion 6 from the Chebyshev disk of criterion 7.

    Asked from both sides, because a cell in the geometric exception set can
    name a neighbour that cannot name it back with a single step: a trapezoid
    has up to three cells along one edge and only one of them is *the* step.
    Adjacency itself is symmetric, so answering from either side is enough.

    Args:
        origin: Atomic index of the first cell.
        destination: Atomic index of the second cell.

    Returns:
        ``True`` when the cells are edge-adjacent.

    Raises:
        ResolutionError: If the cells sit at different resolutions.
    """
    if get_resolution(origin) != get_resolution(destination):
        raise ResolutionError(
            f"{origin!r} and {destination!r} sit at different resolutions; "
            "adjacency is defined between cells of one resolution"
        )
    return destination in _edge_neighbors(origin) or origin in _edge_neighbors(
        destination
    )


# --------------------------------------------------------------------------
# Directed edges
# --------------------------------------------------------------------------

#: Separator between the two endpoints of a directed edge identifier. Chosen
#: because the index grammar never produces it, so an edge cannot be mistaken
#: for a cell and splitting one is unambiguous.
EDGE_SEPARATOR = ">"


def cells_to_directed_edge(origin: str, destination: str) -> str | list[str]:
    """Identifier of the directed edge between two adjacent cells.

    Validates edge adjacency first. Origins and destinations are paired by
    position, so N origins and N destinations give N edges; broadcast a single
    origin by repeating it.

    Args:
        origin: Compositional index of the origin cell or cells.
        destination: Compositional index of the destination cell or cells.

    Returns:
        An edge identifier for a single pair, or a positionally aligned list.

    Raises:
        ValueError: If the two indices hold different cell counts.
        DomainError: If any pair is not edge-adjacent.
    """
    origins, destinations = _atoms(origin), _atoms(destination)
    if len(origins) != len(destinations):
        raise ValueError(
            f"{len(origins)} origins against {len(destinations)} destinations; "
            "the two indices must hold the same number of cells"
        )
    edges: list[str] = []
    for tail, head in zip(origins, destinations):
        if not are_neighbor_cells(tail, head):
            raise DomainError(
                f"{tail!r} and {head!r} do not share an edge, so no directed "
                "edge joins them"
            )
        edges.append(f"{tail}{EDGE_SEPARATOR}{head}")
    return edges[0] if len(edges) == 1 else edges


def directed_edge_to_cells(edge: str) -> tuple[str, str] | list[tuple[str, str]]:
    """Recover the origin and destination of a directed edge.

    Args:
        edge: Directed edge identifier, or several separated by commas.

    Returns:
        ``(origin, destination)`` for a single edge, or a positionally aligned
        list.

    Raises:
        InvalidIndexError: If the identifier does not hold exactly two cells.
    """
    pairs: list[tuple[str, str]] = []
    for token in edge.split(","):
        tail, separator, head = token.partition(EDGE_SEPARATOR)
        if not separator or not tail or not head:
            raise InvalidIndexError(
                f"{token!r} is not a directed edge; expected "
                f"'origin{EDGE_SEPARATOR}destination'"
            )
        pairs.append((tail, head))
    return pairs[0] if len(pairs) == 1 else pairs


def cell_to_edges(cell: str) -> list[str] | list[list[str]]:
    """Every directed edge leaving a cell.

    Four for a cell whose four cardinal steps all land, which is every cell
    outside the geometric exception set. The paper's count of three for a
    prime-meridian triangle is wrong on both halves: the triangle has three
    sides but four edge neighbours, because its base is divided between the
    two quadrants.

    Args:
        cell: Compositional index string.

    Returns:
        An edge list for a single cell, or a positionally aligned list of
        lists.

    Raises:
        DomainError: If the cell is in the geometric exception set, where the
            edge set is not recoverable from the lexical step alone.
    """
    per_cell = [
        [f"{atom}{EDGE_SEPARATOR}{neighbour}" for neighbour in _edge_neighbors(atom)]
        for atom in _atoms(cell)
    ]
    return per_cell[0] if len(per_cell) == 1 else per_cell
