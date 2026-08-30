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
from .hierarchy import _parent_cell, _render, get_children
from .index import decompose, split_components
from .resolutions import get_resolution

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


def get_neighbor(index: str, direction: Direction) -> str | None | list[str | None]:
    """Single adjacent cell in a given direction.

    Directions are lattice directions under quadrant mirroring, so ``"N"``
    means one step toward the pole and ``"E"`` one step away from the prime
    meridian, in every quadrant. The cardinal four share an edge with the
    origin; the diagonal four share only a vertex.

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
    results = [_neighbor_of_atom(atom, direction) for atom in atoms]
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
    if answer is not None and _crosses_base_cell(cell, answer):
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
    return deflect(cell, direction)


def _crosses_base_cell(cell: str, answer: str) -> bool:
    """Whether a step left the resolution-1 cell it started in."""
    return split_components(cell)[:2] != split_components(answer)[:2]


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

    Checked only when the base cell changed, which is the only place a frame
    can have changed. Steps inside one resolution-1 cell are pure arithmetic
    and pay nothing for this.
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
        if _is_exceptional(_render(components[:depth])):
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
        return _render(split_components(parent) + [_encode(column, row, side)])

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
    return _render(split_components(neighbour_parent) + [code])


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

    Three quadrants are searched: the cell's own, its mirror across the
    meridian, and its mirror across the equator. No cell touches the quadrant
    diagonally opposite.
    """
    quadrant, _, row = _res1_parts(cell)
    west, _, east, _ = _ring(cell).bounds
    # Both readings of the cell's longitude span: as written, and as the
    # mirrored quadrant writes it. An extension cell's ring runs past 180, and
    # the cells below it on an ordinary row are addressed from the other side
    # of the seam, at 360 minus that.
    reach = [abs(west), abs(east), abs(360.0 - abs(west)), abs(360.0 - abs(east))]
    span = (min(reach), max(reach))
    candidates: list[str] = []
    for other in (quadrant, _ACROSS_MERIDIAN[quadrant], _ACROSS_EQUATOR[quadrant]):
        for other_row in (row - 1, row, row + 1):
            if not 0 <= other_row <= _POLAR_ROW:
                continue
            last = _last_column(other, other_row)
            step = _degrees_per_column(other, other_row)
            columns = {
                *range(0, min(last, _WINDOW) + 1),
                *range(max(0, last - _SEAM_REACH), last + 1),
                *range(
                    max(0, int(span[0] / step) - _WINDOW),
                    min(last, int(span[1] / step) + _WINDOW) + 1,
                ),
            }
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


def _deep_lattice_ij(cell: str) -> tuple[int, int, int]:
    """A cell's column and row on the lattice of its own resolution.

    Built by folding each refinement code into the resolution-1 pair: a
    quaternary level multiplies by two and a quinary one by five, and the
    code's own column and row are added. The product of those factors is
    returned as well, because it is how many cells of this resolution fill one
    resolution-1 cell along each axis, and the distance formula needs it.

    The lattice is meridian-anchored at resolution 1 and parent-anchored below
    it, which is why the poleward step moves a column at the top level and
    none inside a parent. Both are the same geometry; the frames differ.
    """
    components = split_components(cell)
    _, column, row = _res1_parts(_render(components[:2]))
    factor = 1
    for depth, code in enumerate(components[2:], start=2):
        side = _grid_side(depth)
        code_column, code_row = _decode(code, side)
        column, row = column * side + code_column, row * side + code_row
        factor *= side
    return column, row, factor


def _basis_coordinates(cell: str) -> tuple[int, int]:
    """Displacement of a cell in the lattice's own basis.

    The lattice is spanned by the eastward step ``(1, 0)`` and the poleward
    step ``(-1, 1)``, so a pair decomposes as ``column + row`` eastward steps
    and ``row`` poleward ones. Chebyshev distance is the larger of the two
    magnitudes, Manhattan their sum.

    Below resolution 1 the shear sits at the top level only: a poleward step
    keeps the column inside a parent and moves one whole resolution-1 column
    when it leaves. The correction adds that column back for each
    resolution-1 row crossed, which makes the poleward step cost exactly one
    at every resolution.
    """
    column, row, factor = _deep_lattice_ij(cell)
    return column + factor * (row // factor), row


def grid_distance(origin: str, destination: str, metric: Metric = "chebyshev") -> int:
    """Steps between two cells on the lattice.

    Not a geodesic distance. Use :func:`itacart.geodesy.inverse_geodesic` on
    the anchors or centroids for metric separation.

    Defined at every resolution, within one quadrant. Quadrant mirroring flips
    the sense of the shear, so the poleward step moves one column west in the
    northern hemisphere and one column east in the southern: no single basis
    measures both, and a formula that ignored this would be wrong by twice the
    row difference without ever breaking symmetry.

    Args:
        origin: Atomic index of the first cell.
        destination: Atomic index of the second cell.
        metric: Lattice metric to measure under.

    Returns:
        Step count, zero when the cells coincide.

    Raises:
        ResolutionError: If the cells sit at different resolutions, or are
            quadrants.
        DomainError: If the two cells sit in different quadrants.
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
    if split_components(origin)[0] != split_components(destination)[0]:
        raise DomainError(
            f"{origin!r} and {destination!r} sit in different quadrants; the "
            "lattice basis mirrors across the meridian and the equator, so no "
            "single basis measures both"
        )
    here = _basis_coordinates(origin)
    there = _basis_coordinates(destination)
    east, pole = abs(there[0] - here[0]), abs(there[1] - here[1])
    return max(east, pole) if metric == "chebyshev" else east + pole


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
