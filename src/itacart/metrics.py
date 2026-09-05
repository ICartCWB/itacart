"""Cell shape metrics, for comparison against other DGGS implementations.

The two metrics of :func:`compactness` and :func:`normalized_cell_area`
are transcribed from Kmoch, Vasilyev, Virro and Uuemaa (2022), "Area and
shape distortions in open-source discrete global grid systems", *Big
Earth Data* 6(3), 256-275, doi:10.1080/20964471.2022.2094926. That paper
measures H3, S2, OpenEAGGR, rHEALPix and DGGRID on exactly these two
numbers, and the point of reproducing its conventions rather than better
ones is that the result has to be comparable with its figures. Where a
geodesic quantity would have been more accurate, the paper's planar
convention is kept and the difference is stated.

:func:`cell_base_angle` is not from that paper. It answers the question
the ITACaRT paper leaves qualitative: the cells are identical
parallelograms with a 45-degree acute angle on the sinusoidal plane, and
Figure 8 of Silva, Dietzsch and Shiguemori (2025) shows in pictures what
becomes of that angle on the ellipsoid. This function is that picture as
a number.

Three properties of the results are worth knowing before reading them.
Compactness is not monotone in latitude: it rises from 0.539 at the
equator to 0.785 at row 900, where the shear of the index cancels the
convergence of the meridians and the cell is a square to five decimal
places, and then falls again. The angle is a function of the lattice
position, not of latitude alone, and it is not monotone either. And
compactness is undefined on the polar cap cells, for a reason documented
in :func:`compactness`; the angle is defined there and is zero.
"""

from __future__ import annotations

import math
from typing import Literal, cast

from .cells import cell_to_anchor, cell_to_boundary, cell_to_centroid
from .constants import WGS84_A, WGS84_E2
from .exceptions import GeometryError
from .geodesy import geodetic_to_sinusoidal, inverse_geodesic
from .index import is_atomic, iter_cells, quadrant_of
from .resolutions import effective_cell_area, get_resolution, nominal_cell_area

__all__ = [
    "compactness",
    "cell_base_angle",
    "normalized_cell_area",
]

AngleOutput = Literal["degrees", "radians", "sin"]

_E = math.sqrt(WGS84_E2)

#: Longitudes closer than this are the same meridian. One micro-degree is
#: about 0.11 m at the equator, four orders below the finest cell side.
_LON_TOL = 1e-6

#: Latitudes closer than this are the same parallel, same reasoning.
_LAT_TOL = 1e-6


# --------------------------------------------------------------------------
# Lambert Azimuthal Equal Area, oblique aspect, on the ellipsoid
# --------------------------------------------------------------------------


def _authalic_q(sin_phi: float) -> float:
    """Snyder's ``q``, the authalic area function.

    Equation (3-12) of Snyder, *Map Projections - A Working Manual*, USGS
    Professional Paper 1395 (1987), page 16. ``q`` is proportional to the
    area of the zone between the equator and the given latitude, which is
    what makes the projection built on it equal-area.
    """
    t = _E * sin_phi
    return (1.0 - WGS84_E2) * (
        sin_phi / (1.0 - t * t) - (1.0 / (2.0 * _E)) * math.log((1.0 - t) / (1.0 + t))
    )


_QP = _authalic_q(1.0)

#: Radius of the sphere of equal surface area, Snyder equation (3-13).
_RQ = WGS84_A * math.sqrt(_QP / 2.0)


def _authalic_latitude(phi: float) -> float:
    """Authalic latitude in radians, Snyder equation (3-11)."""
    ratio = _authalic_q(math.sin(phi)) / _QP
    return math.asin(min(1.0, max(-1.0, ratio)))


def _laea(
    lon_deg: float, lat_deg: float, lon0_deg: float, lat0_deg: float
) -> tuple[float, float]:
    """Project one point into the oblique ellipsoidal LAEA of an origin.

    Snyder equations (24-13) through (24-20), pages 187-190. The origin
    is the cell's own centroid, which is what step 3 of section 2.2 of
    Kmoch et al. (2022) prescribes: each cell is reprojected individually
    with its centroid as the origin of the projection.

    The projection is equal-area, so the planar area of a cell's ring
    agrees with :func:`itacart.resolutions.effective_cell_area` to within
    the error of drawing the cell's sides as straight lines in the plane.
    That agreement is asserted in the test suite, and it is the control
    that says this transcription is the projection it claims to be.
    """
    lam = math.radians(lon_deg - lon0_deg)
    beta = _authalic_latitude(math.radians(lat_deg))
    phi0 = math.radians(lat0_deg)
    beta0 = _authalic_latitude(phi0)
    d = (WGS84_A * math.cos(phi0)) / (
        math.sqrt(1.0 - WGS84_E2 * math.sin(phi0) ** 2) * _RQ * math.cos(beta0)
    )
    denominator = (
        1.0
        + math.sin(beta0) * math.sin(beta)
        + math.cos(beta0) * math.cos(beta) * math.cos(lam)
    )
    b = _RQ * math.sqrt(2.0 / denominator)
    x = b * d * math.cos(beta) * math.sin(lam)
    y = (b / d) * (
        math.cos(beta0) * math.sin(beta)
        - math.sin(beta0) * math.cos(beta) * math.cos(lam)
    )
    return x, y


# --------------------------------------------------------------------------
# Ring predicates and planar measures
# --------------------------------------------------------------------------


def _normalized_longitude(lon_deg: float) -> float:
    """Fold a longitude into ``[-180, 180)`` so that 180 and -180 agree."""
    return ((lon_deg + 180.0) % 360.0) - 180.0


def _repeated_vertex(ring: list[tuple[float, float]]) -> tuple[int, int] | None:
    """Return the first pair of ring positions holding the same point.

    A ring that visits a point twice is not a simple closed curve, and
    the area enclosed by it is not the area of the cell. The polar cap
    cells are in exactly that state: their ring runs pole, -180, +180 at
    one latitude, and the two ends of the base are the same point, so the
    limiting parallel is missing from the boundary altogether.
    """
    for i, (lon_i, lat_i) in enumerate(ring):
        for j in range(i + 1, len(ring)):
            lon_j, lat_j = ring[j]
            same_parallel = abs(lat_i - lat_j) <= _LAT_TOL
            same_meridian = abs(_normalized_longitude(lon_i - lon_j)) <= _LON_TOL
            if same_parallel and same_meridian:
                return i, j
    return None


def _planar_area(points: list[tuple[float, float]]) -> float:
    """Shoelace area of a closed planar ring, in square metres."""
    total = 0.0
    for (x1, y1), (x2, y2) in zip(points, points[1:] + points[:1]):
        total += x1 * y2 - x2 * y1
    return abs(total) / 2.0


def _planar_perimeter(points: list[tuple[float, float]]) -> float:
    """Perimeter of a closed planar ring, in metres."""
    return sum(
        math.hypot(x2 - x1, y2 - y1)
        for (x1, y1), (x2, y2) in zip(points, points[1:] + points[:1])
    )


def _ring_of(cell: str) -> list[tuple[float, float]]:
    return cast("list[tuple[float, float]]", cell_to_boundary(cell))


def _projected_ring(cell: str) -> list[tuple[float, float]]:
    """The cell's ring in its own LAEA, refusing a non-simple ring."""
    ring = _ring_of(cell)
    repeated = _repeated_vertex(ring)
    if repeated is not None:
        i, j = repeated
        raise GeometryError(
            f"{cell!r} has no compactness: its boundary visits the same point "
            f"at positions {i} and {j}, so the ring is not a simple closed "
            f"curve and the area it encloses is not the area of the cell"
        )
    lon0, lat0 = cast("tuple[float, float]", cell_to_centroid(cell))
    return [_laea(lon, lat, lon0, lat0) for lon, lat in ring]


# --------------------------------------------------------------------------
# Public metrics
# --------------------------------------------------------------------------


def compactness(cell: str) -> float | list[float]:
    """Isoperimetric quotient of a cell, ``4 * pi * A / p ** 2``.

    A transcription of equation (1) of section 2.2 of Kmoch, Vasilyev,
    Virro and Uuemaa (2022), *Big Earth Data* 6(3), 256-275. The metric
    is unitless, a circle scores 1, and less compact shapes fall towards
    zero.

    ``A`` and ``p`` are **planar**, measured in the oblique ellipsoidal
    Lambert Azimuthal Equal Area projection centred on the cell's own
    centroid, which is step 3 of the paper's workflow. They are not
    geodesic and not measured on the sinusoidal plane. For a cell of 10
    km the difference against the geodesic quantities is of the order of
    ``(d / R) ** 2``, about ``2e-6``; the paper's convention is kept
    anyway, because the purpose of the number is to be comparable with
    the paper's figures.

    Values to expect, all measured: a regular hexagon scores 0.9069, a
    square 0.7854 and an equilateral triangle 0.6046, which is the ladder
    the paper's Figure 11(b) puts H3, S2 and the triangular grids on.
    ITACaRT cells fall between the triangular and the square grids, from
    0.539 at the equator to 0.785 at row 900, where the cell is a square
    to five decimal places. That is the price of the 45-degree shear and
    it is a number to declare, not a defect to hide.

    Args:
        cell: Compositional index string.

    Returns:
        The quotient for a single terminal cell, or a positionally
        aligned list when the index holds several.

    Raises:
        GeometryError: If the cell's boundary is not a simple closed
            curve. The polar cap cells are the only members of that
            class: their ring visits the same point twice and omits the
            limiting parallel, so the enclosed area is not the cell's.
            Step 2 of the paper's own workflow removes cells whose
            geometry is invalid for exactly this kind of reason.
    """
    values = []
    for atom in iter_cells(cell):
        points = _projected_ring(atom)
        perimeter = _planar_perimeter(points)
        values.append(4.0 * math.pi * _planar_area(points) / (perimeter * perimeter))
    return values[0] if is_atomic(cell) else values


def cell_base_angle(
    cell: str,
    output: AngleOutput = "degrees",
) -> float | list[float]:
    """Ellipsoidal image of the 45-degree angle of the sinusoidal plane.

    On the sinusoidal plane every ITACaRT cell is the same parallelogram:
    Figure 1 of Silva, Dietzsch and Shiguemori (2025) places its vertices
    at ``{x, y}``, ``{x + l, y}``, ``{x, y + l}`` and ``{x - l, y + l}``,
    so the base equals the height and the acute angle is exactly 45
    degrees. This function reports what that angle becomes once the cell
    is carried to the WGS84 ellipsoid.

    It is measured at the **anchor**, between two lines. One is the
    **parallel**, which is not a geodesic: its direction at the anchor is
    the due east and due west tangent, and the arm taken is the one the
    45 degrees opens against, opposite the cell's base. The other is the
    **geodesic** from the anchor to the vertex at the most negative X in
    the cell's own frame -- the corner the shear displaces -- in a
    parallelogram, a triangle or a trapezoid alike.

    Reading the base as a geodesic chord instead of as a parallel biases
    the result by half the convergence of the meridians across the base.
    That bias is zero at the equator, where a parallel is a geodesic,
    and it grows with latitude: 0.023 degree at row 300 and 0.284 at row
    900, where it is the whole difference between 89.426479 and
    89.710055.

    The angle is a function of ``(X, Y)``, the lattice position, not of
    latitude. Along row 0 it grows with longitude at 3.95e-4 degree per
    degree, from 45.000035 at column 1 to 45.070984 at column 2000, at
    constant latitude throughout. It is also not monotone: it reaches
    89.426479 at row 900, where the cell is a square, and keeps opening
    past a right angle, to 134.375799 at row 950 and 142.6159 in
    Kamchatka, which is the "severely obtuse angle exceeding 135
    degrees" that Figure 8 of the paper shows as a picture. The four
    quadrants agree to nine decimals at every one of those latitudes.

    Triangular cells return 0, and the zero is measured rather than
    declared by a clause. The anchor of a meridian triangle is not a
    vertex of its ring: it sits at longitude 0 in the middle of the base,
    and the vertex at the most negative X is the far end of that base, at
    the same latitude. The side leaving the anchor towards it *is* the
    parallel, so the angle against the parallel is zero. The polar caps
    are triangles whose anchor is on Greenwich at the middle of their
    base, and they fall under the same rule without a clause of their
    own.

    Args:
        cell: Compositional index string.
        output: ``"degrees"`` for the angle, ``"radians"`` for the same
            angle in radians, or ``"sin"`` for its sine. The sine is the
            multiplicand of ``base * side * sin(angle)``, so it is the
            sine of this angle and of no other quantity.

    Returns:
        The angle for a single terminal cell, or a positionally aligned
        list when the index holds several.

    Raises:
        GeometryError: If ``output`` is not one of the three named forms.
    """
    if output not in ("degrees", "radians", "sin"):
        raise GeometryError(
            f"output must be 'degrees', 'radians' or 'sin', not {output!r}"
        )
    values = []
    for atom in iter_cells(cell):
        values.append(_angle_in(_atom_base_angle(atom), output))
    return values[0] if is_atomic(cell) else values


def _local_offsets(
    cell: str, ring: list[tuple[float, float]], anchor: tuple[float, float]
) -> list[float]:
    """Each ring vertex's X offset from the anchor, in the cell's own frame.

    Measured on the sinusoidal plane and not in longitude, because the
    45 degrees is a property of that plane: Figure 1 of the paper puts
    the vertices at ``{x, y}``, ``{x + l, y}``, ``{x, y + l}`` and
    ``{x - l, y + l}``, so the offsets are exactly ``0``, ``+l``, ``0``
    and ``-l`` at every latitude. In longitude they are not: the
    meridians converge, and above roughly row 900 the vertex that the
    shear puts to the west of the anchor comes out east of it, which
    would pick the wrong vertex for the very cells where the angle is
    most interesting.

    The sign follows the quadrant, so ``+X`` always runs the way the
    lattice counts columns: east in the eastern quadrants, west in the
    western ones, where the shear falls out negative all the same.
    """
    east = 1.0 if quadrant_of(cell)[1] == "E" else -1.0
    anchor_x, _ = geodetic_to_sinusoidal(*anchor)
    return [east * (geodetic_to_sinusoidal(*v)[0] - anchor_x) for v in ring]


def _atom_base_angle(cell: str) -> float:
    """Angle at the anchor, between the parallel and the leaning side."""
    ring = _ring_of(cell)
    anchor = cast("tuple[float, float]", cell_to_anchor(cell))
    offsets = _local_offsets(cell, ring, anchor)
    target = ring[min(range(len(ring)), key=lambda k: offsets[k])]
    if abs(target[1] - anchor[1]) <= _LAT_TOL:
        return 0.0
    east = 1.0 if quadrant_of(cell)[1] == "E" else -1.0
    reference = 270.0 if east > 0.0 else 90.0
    _, azimuth = inverse_geodesic(*anchor, *target)
    separation = abs(azimuth - reference) % 360.0
    return separation if separation <= 180.0 else 360.0 - separation


def _angle_in(degrees: float, output: AngleOutput) -> float:
    if output == "degrees":
        return degrees
    radians = math.radians(degrees)
    return radians if output == "radians" else math.sin(radians)


def normalized_cell_area(cell: str) -> float | list[float]:
    """Cell area divided by the nominal area of its resolution.

    The second metric of section 2.2 of Kmoch et al. (2022), with one
    documented substitution. The paper divides each cell's area by the
    **mean area of all cells at the same resolution**; that mean is not
    computable here, because resolution 1 alone holds about four million
    cells and every finer level multiplies that. The denominator used
    instead is :func:`itacart.resolutions.nominal_cell_area`, the area
    the resolution table assigns to the level.

    The substitution is safe to state and unsafe to leave unstated. It is
    safe because the sinusoidal plane is equal-area by construction, so
    every unclipped cell has exactly the nominal area and scores exactly
    1.0 -- against plus or minus 50 percent for H3 and S2, whose
    normalized area has a standard deviation of 0.13 and 0.15 in the
    paper. It has to be stated because it is a different denominator, and
    a claim that names its own substitution is worth more than one that
    hides it.

    The exceptions are the boundary families, and they are large.
    Enumerated over all four quadrants at resolution 1: 4000 trapezoids
    range from 0.032263 to 2.067292, and 4004 triangles from 0.121394 to
    1.0, the minimum of the latter being the four polar cap cells.

    Args:
        cell: Compositional index string.

    Returns:
        The ratio for a single terminal cell, or a positionally aligned
        list when the index holds several.
    """
    values = []
    for atom in iter_cells(cell):
        area = cast(float, effective_cell_area(atom))
        values.append(area / nominal_cell_area(get_resolution(atom)))
    return values[0] if is_atomic(cell) else values
