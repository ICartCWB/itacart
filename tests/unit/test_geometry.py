"""Tests for :mod:`itacart.geometry`.

Nothing is portable from ``itacart_core``: it has no coverage for this
module. The 406 figure that used to stand here was the size of its whole
suite, not a count of anything reusable.

Two families of test carry most of the weight. The containment chain is
asserted as a set inclusion rather than as a count comparison, because
three counts can nest while the sets do not. And the border families are
enumerated rather than sampled: the last column of a row is one cell in
two thousand, and a random draw finds it only by luck.
"""

from __future__ import annotations

import itertools
import math
from typing import Callable, Iterable

import pytest
from shapely.geometry import (
    GeometryCollection,
    LinearRing,
    LineString,
    MultiPoint,
    MultiPolygon,
    Point,
    Polygon,
    box,
)

import itacart
from itacart import boundary, geometry, hierarchy
from itacart.constants import QUADRANTS
from itacart.exceptions import (
    AntemeridianError,
    DensificationError,
    DomainError,
    GeometryError,
    NonExistentCellError,
    ResolutionError,
    UnsupportedGeometryTypeError,
)

# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------

FIGURE_7A = (
    "NW(0625/0451(1(E1(3(B2(4(A2,B2,B3,B4,C2,C3,C4,C5,D1,D2,D3,D4,"
    "D5,E1,E2,E3,E4,E5)),B3(3(C1,D1,D2,E1,E2,E3,E4)),C2(1(A5,B5,C5,"
    "D4,D5,E4,E5),2,3(A4,A5,B3,B4,B5,C3,C4,C5,D2,D3,D4,D5,E2,E3,E4,"
    "E5),4),C3(1,2(A1,B1,B2,B3,C1,C2,C3,C4,D1,D2,D3,D4,E1,E2,E3,E4,"
    "E5),3,4),D2(1(A1,A2,A3,A4,A5,B2,B3,B4,B5,C3,C4,C5,D5),2,4(A3,"
    "A4,A5,B5)),D3(1,2,3(A1,A2,A3,A4,A5,B1,B2,B3,B4,B5,C2,C3,C4,C5,"
    "D4,D5),4),D4(1(D1,E1),3(A1,B1,C1,D1,E1)))))))"
)
"""Alternative (a) of Figure 7: a compositional fill at resolutions 6 and 7."""

FIGURE_7B = (
    "NW(0625/0451(1(E1(3(C3(2(C5(1(A2(3(C3(3(C4)))))))),"
    "B2(4(A2(1(D5(4(A3(1(D4)))))))),"
    "D2(1(B1(1(B2(4(C2(4(E5))))))),4(B3(2(A4(4(B4(2(D5)))))))),"
    "D3(3(E5(2(D5(1(D5(2(E3)))))))),"
    "D4(3(E1(4(D2(3(C3(4(B1))))))),1(E1(2(D2(4(C5(1(B1)))))))))))))"
)
"""Alternative (b) of Figure 7: vertex representation at resolution 13.

Seven vertices, each a chain from resolution 1 to 13 without a gap, so
each walks all twelve quaternary/quinary alternations. It is the only
artefact in the paper that does, which is why it is the cheapest
end-to-end check of the alphabet the package has.

Transcribed across a page break: the chain opens on the last line of
page 12 and closes on page 13, and a reader stopping at the page
boundary gets a truncated string that still balances its parentheses.
"""


@pytest.fixture
def parcel() -> Polygon:
    """A parcel inside the NW quadrant, clear of every border family."""
    return Polygon(
        [
            (-73.9812, 40.7681),
            (-73.9581, 40.8005),
            (-73.9497, 40.7968),
            (-73.9730, 40.7644),
        ]
    )


@pytest.fixture
def small_square() -> Polygon:
    """A square of a few hundred metres, cheap to fill at resolution 5."""
    return Polygon([(10.00, 45.00), (10.01, 45.00), (10.01, 45.01), (10.00, 45.01)])


# --------------------------------------------------------------------------
# Criterion 1 and 6: the paper's Figure 7
# --------------------------------------------------------------------------


def test_figure_7a_is_a_valid_mixed_resolution_index() -> None:
    """Alternative (a) parses and holds cells at resolutions 6 and 7 only."""
    assert itacart.is_valid_index(FIGURE_7A)
    cells = itacart.decompose(FIGURE_7A)
    assert len(cells) == 114
    assert sorted({itacart.get_resolution(cell) for cell in cells}) == [6, 7]


def test_figure_7a_refills_from_its_own_footprint() -> None:
    """Filling the footprint of Figure 7(a) returns Figure 7(a)'s cells.

    The paper publishes an index, not the polygon behind it, so this is
    a fixed point -- ``fill(union(published)) == published`` -- and not
    a check that the author's parcel produces the author's index. The
    stronger claim is not available from the paper, and saying so is
    part of the evidence.
    """
    from shapely.ops import unary_union

    cells = itacart.decompose(FIGURE_7A)
    footprint = unary_union([itacart.cell_to_polygon(cell) for cell in cells])
    refilled = set(itacart.decompose(itacart.polyfill(footprint, 7)))
    expected = set(itacart.uncompact_cells(FIGURE_7A, 7))
    assert refilled == expected


def test_figure_7a_is_reproduced_as_a_string() -> None:
    """The compacted fill reproduces the published index character for character.

    Stronger than the set comparison above, and the strongest check the
    paper allows: the fill, run through the package's own compaction,
    lands on the exact 415-character string the article prints, so the
    package agrees with the authors on where to compact and not merely
    on which ground is covered.
    """
    from shapely.ops import unary_union

    cells = itacart.decompose(FIGURE_7A)
    footprint = unary_union([itacart.cell_to_polygon(cell) for cell in cells])
    assert itacart.compact_cells(itacart.polyfill(footprint, 7)) == FIGURE_7A


def test_geometric_compaction_is_weaker_than_lexical_compaction() -> None:
    """``compact=True`` collapses less than ``compact_cells`` does.

    The fill collapses a node when its plane square is contained in the
    projected geometry. ``compact_cells`` asks whether all the children
    are present, and the leaf predicate is the centre, so it answers yes
    where the containment test answers no. The ground covered is the
    same either way; which of the two ``compact=True`` should mean is an
    open contract question and not something this test settles.

    Measured on a parcel rather than on Figure 7(a). The figure diverged
    at exactly one node, and the part of that node lying outside the
    geometry was 9.3e-09 square metres -- a sliver two ten-billionths of
    a metre wide, below the ulp of the plane coordinates. A test built
    on it was pinning floating-point noise, and it broke when an
    unrelated change perturbed the descent. What is pinned here instead
    survives refining the densification by four orders of magnitude:
    the excess stays at 495 square metres over 24 nodes, between 2.1e-02
    and 59 square metres each, none of them near the noise floor.
    """
    parcel = Polygon(
        [
            (-73.9812, 40.7681),
            (-73.9581, 40.8005),
            (-73.9497, 40.7968),
            (-73.9730, 40.7644),
        ]
    )
    plain = itacart.polyfill(parcel, 7)
    geometric = itacart.polyfill(parcel, 7, compact=True)
    lexical = itacart.compact_cells(plain)

    assert itacart.count_cells(lexical) < itacart.count_cells(geometric)
    assert set(itacart.uncompact_cells(geometric, 7)) == set(itacart.decompose(plain))
    assert set(itacart.uncompact_cells(lexical, 7)) == set(itacart.decompose(plain))

    divergent = set(itacart.decompose(lexical)) - set(itacart.decompose(geometric))
    assert divergent
    leaves = set(itacart.decompose(plain))
    for node in divergent:
        assert set(itacart.uncompact_cells(node, 7)) <= leaves, node


def test_figure_7a_lies_in_the_nw_quadrant() -> None:
    """The fixture is mirrored, and the fill has to honour that."""
    assert all(cell.startswith("NW(") for cell in itacart.decompose(FIGURE_7A))


def test_figure_7b_walks_all_twelve_alternations() -> None:
    """Alternative (b) is a single chain from resolution 1 to 13."""
    cells = itacart.decompose(FIGURE_7B)
    assert cells
    assert all(itacart.get_resolution(cell) == 13 for cell in cells)


def test_figure_7b_alternates_quaternary_and_quinary() -> None:
    """Each level of the chain uses the alphabet its resolution declares."""
    for cell in itacart.decompose(FIGURE_7B):
        codes = itacart.split_components(cell)[2:]
        for offset, code in enumerate(codes):
            level = 2 + offset
            if itacart.linear_refinement_ratio(level) == 2:
                assert code in {"1", "2", "3", "4"}, (cell, level, code)
            else:
                assert code[0] in "ABCDE" and code[1] in "12345", (cell, level, code)


def test_figure_7b_round_trips_through_vertex_to_cell() -> None:
    """Rebuilding and requantizing (b)'s vertices returns the same cells."""
    cells = itacart.decompose(FIGURE_7B)
    ring = itacart.cells_to_geometry(cells, "Polygon")
    assert itacart.vertex_to_cell(ring, 13) == cells


# --------------------------------------------------------------------------
# Criterion 2: the containment chain
# --------------------------------------------------------------------------


def test_containment_chain_is_a_set_inclusion(parcel: Polygon) -> None:
    """contains is inside center is inside intersects, as sets."""
    sets = {
        mode: set(itacart.decompose(itacart.polyfill(parcel, 6, containment=mode)))
        for mode in ("contains", "center", "intersects")
    }
    assert sets["contains"] < sets["center"] < sets["intersects"]


@pytest.mark.parametrize("resolution", [4, 5, 6])
def test_containment_chain_holds_at_several_resolutions(
    parcel: Polygon, resolution: int
) -> None:
    """The chain is a property of the descent, not of one resolution."""
    inner = set(
        itacart.decompose(itacart.polyfill(parcel, resolution, containment="contains"))
    )
    middle = set(
        itacart.decompose(itacart.polyfill(parcel, resolution, containment="center"))
    )
    outer = set(
        itacart.decompose(
            itacart.polyfill(parcel, resolution, containment="intersects")
        )
    )
    assert inner <= middle <= outer


def test_unknown_containment_mode_is_rejected(parcel: Polygon) -> None:
    with pytest.raises(ValueError, match="containment must be"):
        itacart.polyfill(parcel, 5, containment="touches")  # type: ignore[arg-type]


# --------------------------------------------------------------------------
# Criterion 3: count times nominal area against the equal-area plane
# --------------------------------------------------------------------------


def test_area_agreement_converges_under_refinement(parcel: Polygon) -> None:
    """The area error falls as the cell shrinks, for a fixed parcel.

    The reference is the polygon's area on the parallels plane, which is
    the equal-area image of its ellipsoidal area: that identity is the
    paper's central claim and the reason a cell count is an area at all.

    This is a regression property of the package, not a statement about
    how a resolution should be chosen. It varies the resolution against
    a fixed parcel, which is the opposite of the cadastral situation,
    where the resolution follows from the mapping scale and its
    positional tolerance and is fixed before the parcel is measured. The
    operational property is the bound in
    ``test_area_error_is_bounded_by_the_outline_band``.
    """
    reference = geometry._project(parcel).area
    errors = []
    for resolution in (5, 6, 7, 8):
        counted = itacart.count_internal_cells(parcel, resolution)
        area = counted * itacart.nominal_cell_area(resolution)
        errors.append(abs(area - reference) / reference)
    assert errors[-1] < errors[0] / 100.0


def test_area_error_is_bounded_by_the_outline_band(parcel: Polygon) -> None:
    """The error never exceeds the area of the cells the outline crosses.

    A cell the outline crosses is counted or dropped by which side its
    centre falls; every other cell is decided exactly. The total error is
    therefore bounded by the band of crossed cells, whose area is at most
    the perimeter times the cell side.

    This is the property a cadastral user needs, because it is the one
    that holds at a resolution they did not choose.
    """
    plane = geometry._project(parcel)
    reference = plane.area
    for resolution in (5, 6, 7):
        side = itacart.cell_size(resolution)
        band = plane.exterior.length * side
        counted = itacart.count_internal_cells(parcel, resolution)
        area = counted * itacart.nominal_cell_area(resolution)
        assert abs(area - reference) <= band


def test_the_bound_holds_at_a_prescribed_resolution_across_parcel_sizes() -> None:
    """At a fixed resolution the bound holds for every parcel size.

    The cadastral direction: the resolution is prescribed by the mapping
    scale, and the parcel's size and perimeter then determine the
    residual. Enumerated over five sizes spanning a factor of sixteen,
    at one resolution, because that is the axis the practitioner varies
    and the earlier tests do not.

    Only the bound is asserted. Where inside it a parcel lands is not
    predictable: measured over 24 phases at each of seven sizes, the
    RMS, median and ninetieth-percentile estimators of the scaling
    exponent give 0.61, -0.47 and 0.00. There is no power law to pin.
    """
    resolution = 7
    side = itacart.cell_size(resolution)
    for length in (200, 400, 800, 1600, 3200):
        delta_lat = length / 111_320.0
        delta_lon = delta_lat / math.cos(math.radians(45.0))
        parcel = Polygon(
            [
                (10.0, 45.0),
                (10.0 + delta_lon, 45.0),
                (10.0 + delta_lon, 45.0 + delta_lat),
                (10.0, 45.0 + delta_lat),
            ]
        )
        plane = geometry._project(parcel)
        counted = itacart.count_internal_cells(parcel, resolution)
        area = counted * itacart.nominal_cell_area(resolution)
        assert abs(area - plane.area) <= plane.exterior.length * side


def _plane_perimeter(plane: object) -> float:
    """Total ring length of a projected polygon, holes included."""
    return plane.exterior.length + sum(  # type: ignore[attr-defined]
        hole.length for hole in plane.interiors  # type: ignore[attr-defined]
    )


_BOUND_RADIUS_IN_SIDES = math.sqrt(5.0)
"""Radius of the outline neighbourhood, in cell sides.

The radius is the diameter of a cell, and the cell is not a square. On
the plane of the parallels it has base ``(s, 0)`` and side ``(-s, s)``,
so its longer diagonal is ``s*sqrt(5)`` and not the ``s*sqrt(2)`` a
square would give. The value is asserted against a measured cell by
``test_the_bound_radius_is_the_measured_cell_diameter`` so that it
cannot drift back to a hand-written number.
"""


def _residual_within_derived_bound(parcel: Polygon, resolution: int) -> bool:
    """Whether the count's area error respects the proved bound.

    The proof, in two lines. A cell whose closure does not meet the
    outline is decided exactly, so only cells meeting the outline can be
    misclassified. Every such cell lies inside the neighbourhood of the
    outline whose radius ``r`` is the cell diameter, and that
    neighbourhood has area at most ``2*r*P + pi*r*r``. Substituting
    ``r = s*sqrt(5)`` gives

        |n*a - A| <= 2*sqrt(5)*P*s + 5*pi*s*s

    which holds for any rectifiable outline, convex or not, holed or not.
    """
    plane = geometry._project(parcel)
    side = itacart.cell_size(resolution)
    perimeter = _plane_perimeter(plane)
    radius = _BOUND_RADIUS_IN_SIDES * side
    bound = 2.0 * radius * perimeter + math.pi * radius * radius
    counted = itacart.count_internal_cells(parcel, resolution)
    area = counted * itacart.nominal_cell_area(resolution)
    return abs(area - plane.area) <= bound


def _metres(value: float) -> float:
    """Degrees of latitude for a length in metres."""
    return value / 111_320.0


def _degrees_lon(value: float, latitude: float = 45.0) -> float:
    """Degrees of longitude for a length in metres at ``latitude``."""
    return _metres(value) / math.cos(math.radians(latitude))


def _extended_parcel_family() -> list[tuple[str, Polygon]]:
    """Parcels the rectangle family does not cover.

    The 168 parcels behind the residual measurement were all convex
    axis-aligned rectangles far larger than a cell, which is a narrow
    family to draw a bound from. These add a hole, a reentrant outline, a
    sliver, a parcel smaller than one cell, a triangle and a star.
    """
    shell = [
        (10.0, 45.0),
        (10.0 + _degrees_lon(800), 45.0),
        (10.0 + _degrees_lon(800), 45.0 + _metres(800)),
        (10.0, 45.0 + _metres(800)),
    ]
    hole = [
        (10.0 + _degrees_lon(300), 45.0 + _metres(300)),
        (10.0 + _degrees_lon(500), 45.0 + _metres(300)),
        (10.0 + _degrees_lon(500), 45.0 + _metres(500)),
        (10.0 + _degrees_lon(300), 45.0 + _metres(500)),
    ]
    star = []
    for step in range(12):
        angle = 2.0 * math.pi * step / 12.0
        radius = 400.0 if step % 2 == 0 else 150.0
        star.append(
            (
                10.0 + _degrees_lon(radius * math.cos(angle)),
                45.0 + _metres(radius * math.sin(angle)),
            )
        )
    return [
        ("holed", Polygon(shell, [hole])),
        (
            "reentrant",
            Polygon(
                [
                    (10.0, 45.0),
                    (10.0 + _degrees_lon(800), 45.0),
                    (10.0 + _degrees_lon(800), 45.0 + _metres(300)),
                    (10.0 + _degrees_lon(300), 45.0 + _metres(300)),
                    (10.0 + _degrees_lon(300), 45.0 + _metres(800)),
                    (10.0, 45.0 + _metres(800)),
                ]
            ),
        ),
        (
            "sliver",
            Polygon(
                [
                    (10.0, 45.0),
                    (10.0 + _degrees_lon(2000), 45.0),
                    (10.0 + _degrees_lon(2000), 45.0 + _metres(20)),
                    (10.0, 45.0 + _metres(20)),
                ]
            ),
        ),
        (
            "smaller than a cell",
            Polygon(
                [
                    (10.0, 45.0),
                    (10.0 + _degrees_lon(3), 45.0),
                    (10.0 + _degrees_lon(3), 45.0 + _metres(3)),
                    (10.0, 45.0 + _metres(3)),
                ]
            ),
        ),
        (
            "triangle",
            Polygon(
                [
                    (10.0, 45.0),
                    (10.0 + _degrees_lon(600), 45.0),
                    (10.0, 45.0 + _metres(600)),
                ]
            ),
        ),
        ("star", Polygon(star)),
    ]


@pytest.mark.parametrize(
    "name, parcel",
    _extended_parcel_family(),
    ids=[name for name, _ in _extended_parcel_family()],
)
def test_the_derived_bound_holds_beyond_the_rectangle_family(
    name: str, parcel: Polygon
) -> None:
    """The proved bound holds for holes, reentrant outlines and slivers."""
    assert _residual_within_derived_bound(parcel, 7), name


@pytest.mark.parametrize(
    "name, parcel",
    _extended_parcel_family(),
    ids=[name for name, _ in _extended_parcel_family()],
)
def test_the_tighter_observed_bound_also_holds(name: str, parcel: Polygon) -> None:
    """``P*s`` is tighter than the proof gives, and is not a theorem.

    The derived bound carries a factor of ``2*sqrt(5)`` and an ``s*s``
    term. ``P*s`` is about 4.5 times tighter and every parcel measured
    respects it -- rectangles, holes, slivers, and combs of four to
    sixty-four teeth one cell wide, worst ratio 0.34. No counterexample
    was found, which is not the same as a proof, and the distinction is
    recorded rather than blurred.
    """
    plane = geometry._project(parcel)
    side = itacart.cell_size(7)
    counted = itacart.count_internal_cells(parcel, 7)
    area = counted * itacart.nominal_cell_area(7)
    assert abs(area - plane.area) <= _plane_perimeter(plane) * side, name


def test_the_bound_radius_is_the_measured_cell_diameter() -> None:
    """The neighbourhood radius is measured off a cell, not off a square.

    ``s*sqrt(2)`` is the diagonal of a square and the cell is not one, so
    a bound built on it is understated by a factor of 1.58 and is not the
    bound the proof gives. Measuring the diameter here keeps the constant
    honest: writing it by hand again fails this test.
    """
    plane = geometry._project(Polygon(itacart.cell_to_boundary("NE(0500/0300)")))
    corners = list(plane.exterior.coords)[:-1]
    diameter = max(math.dist(a, b) for a, b in itertools.combinations(corners, 2))
    assert diameter / itacart.cell_size(1) == pytest.approx(
        _BOUND_RADIUS_IN_SIDES, rel=1e-9
    )
    assert _BOUND_RADIUS_IN_SIDES > math.sqrt(2.0)


def test_the_bound_says_nothing_about_a_parcel_narrower_than_a_cell() -> None:
    """Below about one cell across, the bound is true and useless.

    A three-metre parcel at ten-metre cells has a bound of thirteen times
    its own area, so it constrains nothing. The bound informs only while
    ``P*s < A``, which is the same as saying the parcel is more than
    roughly one cell wide. Worth pinning, because a bound quoted outside
    its useful range is how a specification becomes decoration.
    """
    tiny = dict(_extended_parcel_family())["smaller than a cell"]
    plane = geometry._project(tiny)
    side = itacart.cell_size(7)
    assert _plane_perimeter(plane) * side > plane.area
    assert _residual_within_derived_bound(tiny, 7)


def test_a_positional_tolerance_selects_one_resolution() -> None:
    """The ladder a tolerance picks from is decimal and unambiguous.

    The selection rule is the *coarsest* cell that does not exceed the
    tolerance, not the finest: every resolution below the tolerance
    satisfies it, and taking the finest would always land on resolution
    13 and over-resolve the survey. The tolerance itself comes from the
    applicable standard and is deliberately not encoded here.
    """
    for tolerance, expected in (
        (0.10, 11),
        (0.25, 11),
        (0.50, 10),
        (1.00, 9),
        (2.50, 9),
        (5.00, 8),
    ):
        chosen = min(
            r
            for r in range(1, itacart.MAX_RESOLUTION + 1)
            if itacart.cell_size(r) <= tolerance
        )
        assert chosen == expected, (tolerance, chosen)
        assert itacart.cell_size(chosen) <= tolerance
        assert itacart.cell_size(chosen - 1) > tolerance


def test_count_matches_the_centre_mode_fill(parcel: Polygon) -> None:
    """The fast path counts exactly what the slow path names.

    This is the test that keeps ``count_internal_cells`` honest: it is a
    different code path with a different accumulator, and the only thing
    tying it to ``polyfill`` is that they agree.
    """
    for resolution in (4, 5, 6):
        named = itacart.count_cells(
            itacart.polyfill(parcel, resolution, containment="center")
        )
        assert itacart.count_internal_cells(parcel, resolution) == named


# --------------------------------------------------------------------------
# Criterion 4: the count does not materialise indices
# --------------------------------------------------------------------------


_PREVIOUSLY_LISTED_DOORS: tuple[tuple[str, str], ...] = (
    ("itacart.cells", "_from_path"),
    ("itacart.cells", "_quantize"),
    ("itacart.cells", "geo_to_cell"),
    ("itacart.cells", "sinusoidal_to_cell"),
    ("itacart.hierarchy", "_render"),
    ("itacart.hierarchy", "_descend"),
    ("itacart.hierarchy", "_shift_column"),
    ("itacart.hierarchy", "_parent_cell"),
    ("itacart.hierarchy", "common_ancestor"),
    ("itacart.hierarchy", "compact_cells"),
    ("itacart.index", "_canonical_base"),
    ("itacart.index", "_render_path"),
    ("itacart.index", "_render_node"),
    ("itacart.index", "_render_tree"),
    ("itacart.index", "compose"),
    ("itacart.index", "normalize"),
    ("itacart.topology", "_contact"),
    ("itacart.topology", "_eastern"),
    ("itacart.topology", "_spell"),
    ("itacart.topology", "_encode"),
)
"""The twenty doors the previous rule found, kept as a regression floor.

The rule that produced them matched a bare ``str`` return annotation in
four named modules. Widening it must not lose any of these, which is
what ``test_the_widening_dropped_no_previous_door`` checks.
"""


_DOORS_ON_THE_COUNT_PATH: frozenset[tuple[str, str]] = frozenset(
    {
        ("itacart.boundary", "_zone_of_row"),
        ("itacart.geometry", "_prepare"),
        ("itacart.geometry", "_quadrant_pieces"),
    }
)
"""The only routes to a ``str`` the count legitimately takes.

Measured rather than reasoned: everything the rule finds was armed at
once, and whatever fired was moved here until the count survived. All
three return a label -- an extension zone name or a quadrant name --
and none of them names a cell, so the property under test survives
their being called.

Pinned in both directions, because a list of exemptions is exactly the
kind of thing that grows quietly. ``test_count_never_builds_an_index_fragment``
arms everything else and survives, so nothing else is on the path;
``test_every_exempted_door_is_really_on_the_count_path`` arms each of
these alone and requires it to fire, so nothing here is idle.
"""


def _index_rendering_doors() -> list[tuple[str, str]]:
    """Every route in the package from components to an index string.

    Discovered, not listed, and widened on the three axes where the
    previous rule was narrow.

    First, it matched ``hints.get("return") in ("str", str)``, so only a
    bare ``str`` counted. ``list[str]``, ``Iterator[str]``,
    ``str | None`` and ``tuple[str, ...]`` all escaped it, and among the
    escapees were ``index.decompose``, which materialises the whole list
    of index strings, and ``index.iter_cells``, which is by name the
    door a descent would use.

    Second, it skipped anything that was a class, and with the class
    went every method it holds.

    Third, the universe it swept was the set of module names appearing
    in its own list, so a module missing from the list was invisible to
    the check whose only job was to find omissions. The universe is now
    the package, walked with ``pkgutil``.
    """
    import importlib
    import inspect
    import pkgutil

    def mentions_str(annotation: object) -> bool:
        if annotation is str:
            return True
        return isinstance(annotation, str) and "str" in annotation

    modules = [itacart]
    for info in pkgutil.walk_packages(itacart.__path__, "itacart."):
        modules.append(importlib.import_module(info.name))

    found: set[tuple[str, str]] = set()
    for module in modules:
        module_name = module.__name__
        for attribute, value in list(vars(module).items()):
            members: list[tuple[str, object]] = []
            if isinstance(value, type):
                if getattr(value, "__module__", None) != module_name:
                    continue
                members.extend(
                    (f"{attribute}.{name}", member)
                    for name, member in vars(value).items()
                )
            elif callable(value) and getattr(value, "__module__", None) == module_name:
                members.append((attribute, value))
            for label, member in members:
                target: object = member
                if isinstance(member, property):
                    target = member.fget
                elif isinstance(member, (staticmethod, classmethod)):
                    target = member.__func__
                if not callable(target):
                    continue
                try:
                    hints = inspect.get_annotations(target, eval_str=False)
                except (TypeError, ValueError):  # pragma: no cover - defensive
                    continue
                if mentions_str(hints.get("return")):
                    found.add((module_name, label))
    return sorted(found)


def _arm_doors(doors: Iterable[tuple[str, str]], patch: pytest.MonkeyPatch) -> None:
    """Replace every named door with something that fails if it is called."""
    import importlib

    def detonator(module_name: str, label: str) -> Callable[..., str]:
        def fired(*args: object, **kwargs: object) -> str:
            raise AssertionError(f"the count reached {module_name}.{label}")

        return fired

    for module_name, label in doors:
        module = importlib.import_module(module_name)
        if "." in label:
            class_name, attribute = label.split(".", 1)
            owner = getattr(module, class_name)
            if attribute not in vars(owner):  # pragma: no cover - defensive
                continue
            patch.setattr(owner, attribute, detonator(module_name, label))
        else:
            patch.setattr(module, label, detonator(module_name, label))


def test_the_widening_dropped_no_previous_door() -> None:
    """Widening the rule may add doors; it may not lose one."""
    doors = set(_index_rendering_doors())
    missing = set(_PREVIOUSLY_LISTED_DOORS) - doors
    assert not missing, sorted(missing)


def test_the_door_rule_sees_what_a_bare_str_rule_cannot() -> None:
    """The rule fails if it narrows back to a bare ``str`` return.

    Each of these escaped the previous rule, and between them they pin
    all three widenings: the annotation shapes, the module universe, and
    the sheer size of the result. Twenty was never the number.
    """
    doors = set(_index_rendering_doors())
    for door in (
        ("itacart.index", "decompose"),
        ("itacart.index", "iter_cells"),
        ("itacart.index", "split_components"),
        ("itacart.hierarchy", "get_children"),
        ("itacart.hierarchy", "uncompact_cells"),
        ("itacart.hierarchy", "_descendants_of"),
        ("itacart.topology", "grid_disk"),
        ("itacart.topology", "_atoms"),
        ("itacart.boundary", "child_code"),
        ("itacart.boundary", "_side_quadrant"),
    ):
        assert door in doors, door
    assert len(doors) >= 100, len(doors)


def test_count_never_builds_an_index_fragment(
    parcel: Polygon, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Counting names no cell, proved against every route to a string.

    Patching the descent's own helpers would prove nothing: the property
    is that no *package* route to an index string is reached, and a
    descent that rendered through another module would satisfy a guard
    aimed at itself.

    So every door the rule finds is armed at once, less the three that
    the count is measured to take and that return a label rather than a
    cell name. That is 106 detonators against the 20 this test used to
    arm, and the count has to survive with all of them live.
    """
    doors = [
        door
        for door in _index_rendering_doors()
        if door not in _DOORS_ON_THE_COUNT_PATH
    ]
    assert len(doors) >= 100, len(doors)
    _arm_doors(doors, monkeypatch)
    assert itacart.count_internal_cells(parcel, 7) == 33_258


def test_every_exempted_door_is_really_on_the_count_path(parcel: Polygon) -> None:
    """No exemption is idle: armed alone, each one fires.

    Without this, the exemption list is a place to hide a door that the
    count should not be taking.
    """
    for door in sorted(_DOORS_ON_THE_COUNT_PATH):
        with pytest.MonkeyPatch.context() as patch:
            _arm_doors([door], patch)
            with pytest.raises(AssertionError, match="the count reached"):
                itacart.count_internal_cells(parcel, 7)


def test_count_stack_depth_is_the_resolution_difference(parcel: Polygon) -> None:
    """Recursion depth is bounded by the resolution, not by the cell count."""
    seen: list[int] = []
    original = geometry._count_node

    def _traced(
        prepared: object,
        u0: float,
        v0: float,
        side: float,
        level: int,
        target: int,
    ) -> int:
        seen.append(level)
        return original(prepared, u0, v0, side, level, target)  # type: ignore[arg-type]

    geometry._count_node = _traced  # type: ignore[assignment]
    try:
        itacart.count_internal_cells(parcel, 6)
    finally:
        geometry._count_node = original  # type: ignore[assignment]
    assert max(seen) <= 6
    assert min(seen) == 1


# --------------------------------------------------------------------------
# Criterion 5: densification
# --------------------------------------------------------------------------


def test_densify_segment_respects_the_threshold() -> None:
    points = itacart.densify_segment((0.0, 0.0), (1.0, 0.0), 10_000.0)
    for start, end in zip(points, points[1:]):
        assert itacart.inverse_geodesic(*start, *end)[0] <= 10_000.0


def test_densify_segment_spaces_points_equally_in_geodesic_distance() -> None:
    points = itacart.densify_segment((0.0, 0.0), (2.0, 0.0), 50_000.0)
    legs = [
        itacart.inverse_geodesic(*start, *end)[0]
        for start, end in zip(points, points[1:])
    ]
    assert len(legs) > 1
    assert max(legs) - min(legs) < 1e-6


def test_densify_segment_leaves_a_short_span_alone() -> None:
    span = ((0.0, 0.0), (0.001, 0.0))
    assert itacart.densify_segment(*span, 10_000.0) == [span[0], span[1]]


def test_densify_segment_handles_coincident_endpoints() -> None:
    assert itacart.densify_segment((3.0, 4.0), (3.0, 4.0), 10.0) == [
        (3.0, 4.0),
        (3.0, 4.0),
    ]


def test_densify_segment_rejects_an_unknown_edge_model() -> None:
    with pytest.raises(DensificationError, match="unknown edge model"):
        itacart.densify_segment((0.0, 0.0), (1.0, 0.0), 100.0, edge_model="LOXODROME")


@pytest.mark.parametrize("threshold", [0.0, -1.0, float("inf"), float("nan")])
def test_densify_segment_rejects_a_bad_threshold(threshold: float) -> None:
    with pytest.raises(DensificationError, match="positive finite"):
        itacart.densify_segment((0.0, 0.0), (1.0, 0.0), threshold)


def test_densify_segment_rejects_a_non_numeric_threshold() -> None:
    with pytest.raises(DensificationError, match="must be a number"):
        itacart.densify_segment(
            (0.0, 0.0), (1.0, 0.0), "1000"  # type: ignore[arg-type]
        )


def test_densify_orthodromic_is_idempotent() -> None:
    polygon = Polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)])
    once = itacart.densify_orthodromic(polygon, 20_000.0)
    twice = itacart.densify_orthodromic(once, 20_000.0)
    assert list(once.exterior.coords) == list(twice.exterior.coords)


def test_densify_orthodromic_densifies_holes_too() -> None:
    """The origin drops interior rings, and drops them silently."""
    shell = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]
    hole = [(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)]
    polygon = Polygon(shell, [hole])
    dense = itacart.densify_orthodromic(polygon, 20_000.0)
    assert len(dense.interiors) == 1
    assert len(dense.interiors[0].coords) > len(polygon.interiors[0].coords)


def test_densify_orthodromic_keeps_every_side_under_the_threshold() -> None:
    polygon = Polygon([(0.0, 0.0), (3.0, 0.0), (3.0, 2.0), (0.0, 2.0)])
    dense = itacart.densify_orthodromic(polygon, 30_000.0)
    coords = list(dense.exterior.coords)
    for start, end in zip(coords, coords[1:]):
        assert itacart.inverse_geodesic(*start, *end)[0] <= 30_000.0


def test_densify_orthodromic_rejects_a_non_polygon() -> None:
    with pytest.raises(TypeError, match="expected a shapely Polygon"):
        itacart.densify_orthodromic(
            LineString([(0, 0), (1, 1)])  # type: ignore[arg-type]
        )


def test_densify_orthodromic_passes_an_empty_polygon_through() -> None:
    empty = Polygon()
    assert itacart.densify_orthodromic(empty).is_empty


def test_densify_any_leaves_a_line_alone() -> None:
    line = LineString([(0.0, 0.0), (10.0, 0.0)])
    assert geometry._densify_any(line, 1000.0) is line


def test_densify_any_handles_a_multipolygon() -> None:
    parts = MultiPolygon(
        [
            Polygon([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]),
            Polygon([(5.0, 5.0), (6.0, 5.0), (6.0, 6.0)]),
        ]
    )
    dense = geometry._densify_any(parts, 20_000.0)
    assert isinstance(dense, MultiPolygon)
    assert len(dense.geoms) == 2


def test_densify_ring_passes_a_degenerate_ring_through() -> None:
    assert geometry._densify_ring([(1.0, 2.0)], 100.0) == [(1.0, 2.0)]


def test_auto_segment_is_monotone_and_capped() -> None:
    thresholds = [geometry._auto_segment(r) for r in range(1, 14)]
    assert thresholds == sorted(thresholds, reverse=True)
    assert max(thresholds) <= 1000.0
    assert geometry._auto_segment(13) < 1000.0


# --------------------------------------------------------------------------
# Criterion 6: vertex order and winding
# --------------------------------------------------------------------------


def test_vertex_to_cell_preserves_ring_winding() -> None:
    """A ring and its reverse produce reversed cell sequences."""
    ring = [(10.0, 45.0), (10.001, 45.0), (10.001, 45.001), (10.0, 45.001)]
    forward = itacart.vertex_to_cell(Polygon(ring), 9)
    backward = itacart.vertex_to_cell(Polygon(ring[::-1]), 9)
    assert forward != backward
    assert set(forward) == set(backward)


def test_vertex_to_cell_puts_holes_after_the_exterior() -> None:
    shell = [(10.0, 45.0), (10.01, 45.0), (10.01, 45.01), (10.0, 45.01)]
    hole = [(10.004, 45.004), (10.006, 45.004), (10.006, 45.006), (10.004, 45.006)]
    cells = itacart.vertex_to_cell(Polygon(shell, [hole]), 9)
    assert len(cells) == 8
    assert cells[:4] == itacart.vertex_to_cell(Polygon(shell), 9)


def test_vertex_to_cell_handles_a_point() -> None:
    assert itacart.vertex_to_cell(Point(10.0, 45.0), 9) == [
        itacart.geo_to_cell(10.0, 45.0, 9)
    ]


def test_vertex_to_cell_handles_a_linestring() -> None:
    line = LineString([(10.0, 45.0), (10.01, 45.0), (10.02, 45.0)])
    assert len(itacart.vertex_to_cell(line, 9)) == 3


def test_vertex_to_cell_rejects_an_unsupported_type() -> None:
    with pytest.raises(UnsupportedGeometryTypeError, match="cannot take vertices"):
        itacart.vertex_to_cell(
            MultiPoint([(0, 0), (1, 1)]), 9  # type: ignore[arg-type]
        )


def test_vertex_to_cell_rejects_a_bad_resolution() -> None:
    with pytest.raises(ResolutionError):
        itacart.vertex_to_cell(Point(10.0, 45.0), 0)


def test_cells_to_geometry_inverts_vertex_to_cell_for_a_ring() -> None:
    ring = [(10.0, 45.0), (10.001, 45.0), (10.001, 45.001), (10.0, 45.001)]
    cells = itacart.vertex_to_cell(Polygon(ring), 13)
    rebuilt = itacart.cells_to_geometry(cells, "Polygon")
    assert itacart.vertex_to_cell(rebuilt, 13) == cells


@pytest.mark.parametrize(
    "kind, count, expected",
    [
        ("Point", 1, Point),
        ("MultiPoint", 3, MultiPoint),
        ("LineString", 3, LineString),
        ("LinearRing", 3, LinearRing),
        ("Polygon", 3, Polygon),
    ],
)
def test_cells_to_geometry_builds_each_supported_type(
    kind: str, count: int, expected: type
) -> None:
    cells = [itacart.geo_to_cell(10.0 + i * 0.001, 45.0, 9) for i in range(count)]
    assert isinstance(itacart.cells_to_geometry(cells, kind), expected)


def test_cells_to_geometry_rejects_an_empty_list() -> None:
    with pytest.raises(GeometryError, match="empty cell list"):
        itacart.cells_to_geometry([])


def test_cells_to_geometry_rejects_a_composed_index() -> None:
    cells = [itacart.geo_to_cell(10.0 + i * 0.001, 45.0, 9) for i in range(2)]
    with pytest.raises(GeometryError, match="more than one cell"):
        itacart.cells_to_geometry([itacart.compose(cells)], "Point")


@pytest.mark.parametrize(
    "kind, count",
    [("Point", 2), ("LineString", 1), ("Polygon", 2), ("LinearRing", 2)],
)
def test_cells_to_geometry_rejects_a_count_the_type_cannot_use(
    kind: str, count: int
) -> None:
    cells = [itacart.geo_to_cell(10.0 + i * 0.001, 45.0, 9) for i in range(count)]
    with pytest.raises(GeometryError):
        itacart.cells_to_geometry(cells, kind)


def test_cells_to_geometry_rejects_an_unsupported_type() -> None:
    cells = [itacart.geo_to_cell(10.0, 45.0, 9)]
    with pytest.raises(UnsupportedGeometryTypeError, match="cannot rebuild"):
        itacart.cells_to_geometry(cells, "Tetrahedron")


# --------------------------------------------------------------------------
# Criterion 7: consecutive collapse, non-consecutive kept
# --------------------------------------------------------------------------


def test_consecutive_duplicates_collapse() -> None:
    line = LineString([(10.0, 45.0), (10.0000001, 45.0), (10.01, 45.0)])
    assert len(itacart.vertex_to_cell(line, 9)) == 2


def test_non_consecutive_repeat_is_preserved() -> None:
    """A revisited vertex is a self-touching outline, not survey noise."""
    line = LineString([(10.0, 45.0), (10.01, 45.0), (10.0, 45.0)])
    cells = itacart.vertex_to_cell(line, 9)
    assert len(cells) == 3
    assert cells[0] == cells[2]


def test_dedupe_can_be_switched_off() -> None:
    line = LineString([(10.0, 45.0), (10.0000001, 45.0), (10.01, 45.0)])
    kept = itacart.vertex_to_cell(line, 9, dedupe_consecutive=False)
    assert len(kept) == 3


def test_ring_closing_duplicate_is_dropped_cyclically() -> None:
    """A ring whose last vertex shares the first cell loses the last one."""
    ring = [
        (10.0, 45.0),
        (10.01, 45.0),
        (10.01, 45.01),
        (10.0, 45.01),
        (10.0000001, 45.0),
    ]
    assert len(itacart.vertex_to_cell(Polygon(ring), 9)) == 4


def test_dedupe_consecutive_handles_an_empty_sequence() -> None:
    assert geometry._dedupe_consecutive([], cyclic=True) == []


def test_ring_cells_handles_a_degenerate_ring() -> None:
    assert geometry._ring_cells([], 9, True) == []


# --------------------------------------------------------------------------
# Criterion 8: canonical rings
# --------------------------------------------------------------------------


def _ring_cells_for_canonical() -> list[str]:
    ring = [(10.0, 45.0), (10.001, 45.0), (10.001, 45.001), (10.0, 45.001)]
    return itacart.vertex_to_cell(Polygon(ring), 9)


def test_canonicalize_rings_is_idempotent() -> None:
    rings = [_ring_cells_for_canonical()]
    once = itacart.canonicalize_rings(rings)
    assert itacart.canonicalize_rings(once) == once


def test_canonicalize_rings_agrees_across_starting_vertices() -> None:
    ring = _ring_cells_for_canonical()
    rotations = [ring[k:] + ring[:k] for k in range(len(ring))]
    canonical = {tuple(itacart.canonicalize_rings([r])[0]) for r in rotations}
    assert len(canonical) == 1


def test_canonicalize_rings_starts_at_the_smallest_index() -> None:
    ring = _ring_cells_for_canonical()
    assert itacart.canonicalize_rings([ring])[0][0] == min(ring)


def test_canonicalize_rings_does_not_reverse_winding() -> None:
    """Direction distinguishes an exterior from a hole; rotation does not."""
    ring = _ring_cells_for_canonical()
    forward = itacart.canonicalize_rings([ring])[0]
    backward = itacart.canonicalize_rings([ring[::-1]])[0]
    assert forward != backward


def test_canonicalize_rings_orders_a_ring_that_repeats_its_minimum() -> None:
    """Booth compares whole rotations, so a repeated minimum still orders."""
    ring = ["A", "B", "A", "C"]
    assert itacart.canonicalize_rings([ring])[0] == ["A", "B", "A", "C"]


def test_canonicalize_rings_leaves_a_short_ring_alone() -> None:
    assert itacart.canonicalize_rings([["A"], []]) == [["A"], []]


def test_min_rotation_of_an_empty_ring_is_zero() -> None:
    assert geometry._min_rotation([]) == 0


def test_canonicalize_rings_preserves_ring_order() -> None:
    a = _ring_cells_for_canonical()
    b = list(reversed(a))
    out = itacart.canonicalize_rings([a, b])
    assert len(out) == 2
    assert set(out[0]) == set(a)


# --------------------------------------------------------------------------
# Criterion 9: the antemeridian
# --------------------------------------------------------------------------


def test_geometry_crossing_the_antemeridian_is_refused() -> None:
    crossing = Polygon([(179.9, 10.0), (-179.9, 10.0), (-179.9, 10.1), (179.9, 10.1)])
    with pytest.raises(AntemeridianError, match="180 degrees"):
        itacart.polyfill(crossing, 5)


def test_a_geometry_near_but_not_across_the_line_is_not_refused() -> None:
    """Strictly astride, not merely close."""
    near = Polygon([(179.0, 10.0), (179.2, 10.0), (179.2, 10.2), (179.0, 10.2)])
    assert not itacart.crosses_antemeridian(near)


def test_the_antemeridian_test_is_strict_and_not_a_tolerance() -> None:
    """A geometry a thousandth of a degree from the line still does not cross.

    This is the case that separates the strict reading from a tolerance
    reading: a predicate that refused anything within some epsilon of 180
    would refuse this polygon, whose far edge sits 0.001 degrees short of
    the line. The strict predicate asks whether the outline is astride
    the meridian at the ordinate in question, and it is not.

    The earlier test does not discriminate: at 179.0 both readings agree.
    """
    hugging = Polygon(
        [
            (179.990, 10.0),
            (179.999, 10.0),
            (179.999, 10.01),
            (179.990, 10.01),
        ]
    )
    assert not itacart.crosses_antemeridian(hugging)

    astride = Polygon(
        [(179.999, 10.0), (-179.999, 10.0), (-179.999, 10.01), (179.999, 10.01)]
    )
    assert itacart.crosses_antemeridian(astride)


# --------------------------------------------------------------------------
# The refused border families, enumerated
# --------------------------------------------------------------------------


def test_prime_meridian_column_is_refused() -> None:
    """Column zero holds triangles, not the square the descent tests."""
    astride = Polygon([(-0.01, 45.0), (0.01, 45.0), (0.01, 45.01), (-0.01, 45.01)])
    with pytest.raises(NonExistentCellError, match="prime-meridian column"):
        itacart.polyfill(astride, 5)


def test_every_row_refuses_its_own_last_column() -> None:
    """Enumerated, not sampled: one cell in two thousand per row.

    The polar row is skipped because it is a row whose only column is
    its last one, so it has no interior column to accept and it raises
    for a different and more specific reason.
    """
    side = itacart.cell_size(1)
    checked = 0
    for row in range(0, 1000):
        last = itacart.last_lattice_column("NE", row, side)
        if last <= 1 or itacart.last_lattice_column("NE", row + 1, side) <= 0:
            continue
        with pytest.raises(NonExistentCellError, match="last lattice column"):
            geometry._check_addressable("NE", last, row)
        with pytest.raises(NonExistentCellError, match="last lattice column"):
            geometry._check_addressable("NE", last + 5, row)
        geometry._check_addressable("NE", last - 1, row)
        checked += 1
    assert checked == 999


def test_the_polar_row_is_refused() -> None:
    side = itacart.cell_size(1)
    polar = max(
        row for row in range(1100) if itacart.last_lattice_column("NE", row, side) > 0
    )
    with pytest.raises(DomainError, match="polar row"):
        geometry._check_addressable("NE", 1, polar)


def test_a_row_above_the_pole_addresses_nothing() -> None:
    side = itacart.cell_size(1)
    beyond = (
        max(
            row
            for row in range(1100)
            if itacart.last_lattice_column("NE", row, side) > 0
        )
        + 1
    )
    with pytest.raises(DomainError, match="addresses no cell"):
        geometry._check_addressable("NE", 1, beyond)


@pytest.mark.parametrize("quadrant", ["NE", "NW", "SE", "SW"])
def test_the_screen_is_the_same_in_all_four_quadrants(quadrant: str) -> None:
    """The mirror flips the shear; it does not move the border families."""
    with pytest.raises(NonExistentCellError, match="prime-meridian column"):
        geometry._check_addressable(quadrant, 0, 100)


# --------------------------------------------------------------------------
# Quadrants and the shear
# --------------------------------------------------------------------------


@pytest.mark.parametrize("quadrant", ["NE", "NW", "SE", "SW"])
def test_the_shear_has_unit_jacobian(quadrant: str) -> None:
    """Area on the lattice is area on the ellipsoid, in all four quadrants."""
    a, b, d, e = geometry._QUADRANT_SHEAR[quadrant]
    assert abs(a * e - b * d) == pytest.approx(1.0)


@pytest.mark.parametrize(
    "lon, lat, quadrant",
    [
        (10.0, 45.0, "NE"),
        (-10.0, 45.0, "NW"),
        (10.0, -45.0, "SE"),
        (-10.0, -45.0, "SW"),
    ],
)
def test_a_fill_lands_in_its_own_quadrant(
    lon: float, lat: float, quadrant: str
) -> None:
    parcel = Polygon(
        [(lon, lat), (lon + 0.01, lat), (lon + 0.01, lat + 0.01), (lon, lat + 0.01)]
    )
    cells = itacart.decompose(itacart.polyfill(parcel, 5))
    assert all(cell.startswith(f"{quadrant}(") for cell in cells)


def test_the_lattice_shear_agrees_with_the_quantizer() -> None:
    """The shear used for filling is the one ``cells`` uses for addressing."""
    lon, lat = -10.0, 45.0
    x, y = itacart.geodetic_to_sinusoidal(lon, lat)
    a, b, d, e = geometry._QUADRANT_SHEAR["NW"]
    assert a * x + b * y == pytest.approx(abs(x) + abs(y))
    assert d * x + e * y == pytest.approx(abs(y))


def test_quadrant_pieces_drops_a_polygon_that_only_touches_an_axis() -> None:
    plane = geometry._project(
        Polygon([(0.0, 45.0), (0.01, 45.0), (0.01, 45.01), (0.0, 45.01)])
    )
    quadrants = {quadrant for quadrant, _ in geometry._quadrant_pieces(plane)}
    assert quadrants == {"NE"}


# --------------------------------------------------------------------------
# Descent internals
# --------------------------------------------------------------------------


def test_child_code_follows_the_quaternary_and_quinary_alphabets() -> None:
    assert geometry._child_code(0, 0, 2) == "1"
    assert geometry._child_code(1, 1, 2) == "4"
    assert geometry._child_code(0, 0, 3) == "A1"
    assert geometry._child_code(4, 4, 3) == "E5"


def test_child_code_inverts_the_quantizer_descent() -> None:
    """The code the fill writes is the code the quantizer would read."""
    cell = itacart.geo_to_cell(10.0, 45.0, 7)
    codes = itacart.split_components(cell)[2:]
    for offset, code in enumerate(codes):
        level = 2 + offset
        divisor = itacart.linear_refinement_ratio(level)
        found = [
            (row, column)
            for row in range(divisor)
            for column in range(divisor)
            if geometry._child_code(row, column, level) == code
        ]
        assert len(found) == 1


def test_leaves_between_is_the_product_of_the_ratios() -> None:
    assert geometry._leaves_between(6, 6) == 1
    assert geometry._leaves_between(6, 7) == 25
    assert geometry._leaves_between(5, 7) == 100


def test_expansion_names_every_descendant_once() -> None:
    """The shared expansion holds each descendant of a whole node once."""
    suffix = geometry._expansion(5, 7)
    codes = suffix.strip("()").replace("(", ",").replace(")", ",").split(",")
    leaves = [c for c in codes if c]
    assert len(leaves) == geometry._leaves_between(5, 7) + 4


def test_expansion_is_empty_at_the_target() -> None:
    assert geometry._expansion(7, 7) == ""


def test_expansion_is_shared_between_callers() -> None:
    """Two whole nodes at the same level reuse one string object.

    This is what keeps the interior of a fill cheap: a region wholly
    inside the geometry costs one pointer per node rather than one
    string per cell. Identity, not equality, is the property.
    """
    assert geometry._expansion(9, 11) is geometry._expansion(9, 11)


def test_budget_refuses_past_the_ceiling(monkeypatch: pytest.MonkeyPatch) -> None:
    """The tally is charged before expansion, never after the fact.

    Charging afterwards cannot work: expanding a resolution-1 node to
    resolution 13 builds a string of a trillion cells, and there is no
    counting that once it exists.
    """
    monkeypatch.setattr(geometry, "MAX_FILL_CELLS", 5)
    budget = geometry._Budget()
    budget.charge(5)
    assert budget.spent == 5
    with pytest.raises(GeometryError, match="exceeded 5 cells"):
        budget.charge(1)


def test_project_handles_each_supported_type() -> None:
    assert geometry._project(Point(10.0, 45.0)).geom_type == "Point"
    assert geometry._project(LineString([(10.0, 45.0), (11.0, 45.0)])).geom_type == (
        "LineString"
    )
    shell = [(10.0, 45.0), (11.0, 45.0), (11.0, 46.0)]
    hole = [(10.4, 45.2), (10.6, 45.2), (10.6, 45.4)]
    assert len(geometry._project(Polygon(shell, [hole])).interiors) == 1
    assert (
        geometry._project(
            MultiPolygon(
                [Polygon(shell), Polygon([(20.0, 45.0), (21.0, 45.0), (21.0, 46.0)])]
            )
        ).geom_type
        == "MultiPolygon"
    )


def test_project_rejects_an_unsupported_type() -> None:
    with pytest.raises(UnsupportedGeometryTypeError, match="cannot fill"):
        geometry._project(GeometryCollection([Point(0, 0)]))


# --------------------------------------------------------------------------
# Guards, arguments and edge cases
# --------------------------------------------------------------------------


@pytest.mark.parametrize("resolution", [0, 14, -1])
def test_polyfill_rejects_a_resolution_out_of_range(
    parcel: Polygon, resolution: int
) -> None:
    with pytest.raises(ResolutionError, match="outside"):
        itacart.polyfill(parcel, resolution)


def test_polyfill_rejects_a_non_integer_resolution(parcel: Polygon) -> None:
    with pytest.raises(ResolutionError, match="must be an int"):
        itacart.polyfill(parcel, 5.0)  # type: ignore[arg-type]


def test_polyfill_rejects_a_boolean_resolution(parcel: Polygon) -> None:
    with pytest.raises(ResolutionError, match="must be an int"):
        itacart.polyfill(parcel, True)  # type: ignore[arg-type]


@pytest.mark.parametrize("jobs", [0, -3])
def test_polyfill_rejects_a_bad_worker_count(parcel: Polygon, jobs: int) -> None:
    with pytest.raises(ValueError, match="n_jobs must be >= 1"):
        itacart.polyfill(parcel, 5, n_jobs=jobs)


def test_polyfill_rejects_a_non_integer_worker_count(parcel: Polygon) -> None:
    with pytest.raises(ValueError, match="n_jobs must be an int"):
        itacart.polyfill(parcel, 5, n_jobs=2.0)  # type: ignore[arg-type]


def test_polyfill_rejects_a_boolean_worker_count(parcel: Polygon) -> None:
    with pytest.raises(ValueError, match="n_jobs must be an int"):
        itacart.polyfill(parcel, 5, n_jobs=True)  # type: ignore[arg-type]


def test_polyfill_refuses_an_empty_geometry() -> None:
    with pytest.raises(GeometryError, match="covers no cell"):
        itacart.polyfill(Polygon(), 5)


def test_polyfill_refuses_a_geometry_narrower_than_a_cell() -> None:
    """A sliver under ``contains`` names nothing, and says so."""
    sliver = Polygon(
        [(10.0, 45.0), (10.000001, 45.0), (10.000001, 45.000001), (10.0, 45.000001)]
    )
    with pytest.raises(GeometryError, match="covers no cell"):
        itacart.polyfill(sliver, 3, containment="contains")


def test_polyfill_honours_the_cell_ceiling(
    parcel: Polygon, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(geometry, "MAX_FILL_CELLS", 10)
    with pytest.raises(GeometryError, match="exceeded 10 cells"):
        itacart.polyfill(parcel, 7)


def test_compact_returns_a_mixed_resolution_index(parcel: Polygon) -> None:
    compacted = itacart.polyfill(parcel, 7, compact=True)
    uniform = itacart.polyfill(parcel, 7)
    resolutions = {
        itacart.get_resolution(cell) for cell in itacart.decompose(compacted)
    }
    assert len(resolutions) > 1
    assert itacart.count_cells(compacted) < itacart.count_cells(uniform)


def test_compact_covers_the_same_ground(parcel: Polygon) -> None:
    compacted = itacart.polyfill(parcel, 7, compact=True)
    uniform = set(itacart.decompose(itacart.polyfill(parcel, 7)))
    expanded = set(itacart.uncompact_cells(compacted, 7))
    assert expanded == uniform


def test_threads_and_serial_agree(parcel: Polygon) -> None:
    """The worker count changes the schedule, never the answer."""
    serial = itacart.polyfill(parcel, 6, n_jobs=1)
    threaded = itacart.polyfill(parcel, 6, n_jobs=4)
    assert serial == threaded
    assert itacart.count_internal_cells(parcel, 6, n_jobs=4) == (
        itacart.count_internal_cells(parcel, 6, n_jobs=1)
    )


def test_count_rejects_a_bad_resolution(parcel: Polygon) -> None:
    with pytest.raises(ResolutionError):
        itacart.count_internal_cells(parcel, 0)


def test_count_rejects_a_bad_worker_count(parcel: Polygon) -> None:
    with pytest.raises(ValueError, match="n_jobs"):
        itacart.count_internal_cells(parcel, 5, n_jobs=0)


def test_count_of_an_empty_polygon_is_zero() -> None:
    assert itacart.count_internal_cells(Polygon(), 5) == 0


def test_polyfill_accepts_a_point(small_square: Polygon) -> None:
    """A point is a degenerate fill and still names its cell."""
    filled = itacart.polyfill(Point(10.0, 45.0), 5, containment="intersects")
    assert itacart.count_cells(filled) >= 1


def test_polyfill_accepts_a_linestring() -> None:
    line = LineString([(10.0, 45.0), (10.02, 45.02)])
    filled = itacart.polyfill(line, 5, containment="intersects")
    assert itacart.count_cells(filled) > 1


def test_polyfill_accepts_a_multipolygon() -> None:
    parts = MultiPolygon(
        [
            box(10.00, 45.00, 10.01, 45.01),
            box(10.05, 45.05, 10.06, 45.06),
        ]
    )
    assert itacart.count_cells(itacart.polyfill(parts, 5)) > 1


def test_polyfill_rejects_an_unsupported_geometry() -> None:
    with pytest.raises(UnsupportedGeometryTypeError):
        itacart.polyfill(GeometryCollection([Point(10.0, 45.0)]), 5)


def test_polyfill_output_parses_and_is_atomic_per_cell(parcel: Polygon) -> None:
    filled = itacart.polyfill(parcel, 6)
    assert itacart.is_valid_index(filled)
    for cell in itacart.decompose(filled):
        assert itacart.is_atomic(cell)
        assert itacart.get_resolution(cell) == 6


def test_polyfill_cells_are_real_cells(parcel: Polygon) -> None:
    """Every cell the fill names survives the boundary validator."""
    for cell in itacart.decompose(itacart.polyfill(parcel, 5)):
        assert itacart.is_valid_cell(cell)
        assert not itacart.absorbs_border(cell)


def test_filled_cells_requantize_to_themselves(parcel: Polygon) -> None:
    """A cell's own centroid addresses that cell, at the fill resolution."""
    for cell in itacart.decompose(itacart.polyfill(parcel, 6)):
        lon, lat = itacart.cell_to_centroid(cell)  # type: ignore[misc]
        assert itacart.geo_to_cell(lon, lat, 6) == cell


def test_intersects_covers_the_whole_geometry(small_square: Polygon) -> None:
    """Nothing of the polygon falls outside the intersects fill."""
    from shapely.ops import unary_union

    filled = itacart.polyfill(small_square, 5, containment="intersects")
    cover = unary_union(
        [itacart.cell_to_polygon(cell) for cell in itacart.decompose(filled)]
    )
    assert cover.buffer(1e-9).contains(small_square)


def test_contains_stays_inside_the_geometry(small_square: Polygon) -> None:
    filled = itacart.polyfill(small_square, 5, containment="contains")
    for cell in itacart.decompose(filled):
        assert small_square.buffer(1e-9).contains(itacart.cell_to_polygon(cell))


# --------------------------------------------------------------------------
# Public surface
# --------------------------------------------------------------------------


def test_every_public_name_is_exported() -> None:
    for name in geometry.__all__:
        assert hasattr(itacart, name), name
        assert name in itacart.__all__, name


def test_module_exports_exactly_what_the_phase_promised() -> None:
    assert sorted(geometry.__all__) == sorted(
        [
            "polyfill",
            "count_internal_cells",
            "vertex_to_cell",
            "cells_to_geometry",
            "densify_orthodromic",
            "densify_segment",
            "canonicalize_rings",
        ]
    )


def test_a_row_with_one_column_is_the_polar_row() -> None:
    """The polar row and the single-column row are the same row.

    This is why the polar test comes first in the screen: at resolution
    1 the polar row holds exactly one column, so a last-column test
    placed ahead of it would answer every polar position with the wrong
    diagnosis and the polar branch would be unreachable.
    """
    side = itacart.cell_size(1)
    single = [
        row for row in range(1100) if itacart.last_lattice_column("NE", row, side) == 1
    ]
    polar = [
        row
        for row in range(1100)
        if itacart.last_lattice_column("NE", row, side) > 0
        and itacart.last_lattice_column("NE", row + 1, side) <= 0
    ]
    assert single == polar


def test_min_rotation_orders_a_ring_with_a_repeated_prefix() -> None:
    """Booth's failure function is exercised by a ring with a period."""
    assert geometry._min_rotation(["B", "A", "B", "A"]) == 1
    assert geometry._min_rotation(["A", "A", "B"]) == 0
    assert geometry._min_rotation(["C", "A", "A", "B"]) == 1
    assert geometry._min_rotation(["B", "A", "A", "A"]) == 1


def test_min_rotation_restarts_mid_comparison() -> None:
    """A smaller vertex found part way through a match restarts the scan.

    ``ABA`` and ``BBA`` are the shortest rings that reach that branch:
    the scan is inside a partial match when it meets a vertex smaller
    than the one it was matching against, and has to rewind to it.
    """
    assert geometry._min_rotation(["A", "B", "A"]) == 2
    assert geometry._min_rotation(["B", "B", "A"]) == 2


def test_fill_is_canonical_without_a_sort(parcel: Polygon) -> None:
    """The descent emits siblings in the order the index is read in.

    ``compose`` preserves the order it is handed rather than imposing
    one, so the canonical ordering has to come from the descent itself.
    Asserted against the sorted-and-composed form the fill used to build
    by materialising every atomic index first.
    """
    for resolution in (4, 5, 6):
        for mode in ("contains", "center", "intersects"):
            filled = itacart.polyfill(parcel, resolution, containment=mode)
            assert filled == itacart.compose(sorted(set(itacart.decompose(filled))))


def test_fill_is_canonical_across_base_cells_and_quadrants() -> None:
    """Ordering holds where it is easiest to get wrong: between roots."""
    wide = Polygon([(9.99, 44.99), (10.30, 44.99), (10.30, 45.30), (9.99, 45.30)])
    filled = itacart.polyfill(wide, 3)
    assert filled == itacart.compose(sorted(set(itacart.decompose(filled))))
    assert filled.count("/") > 1  # more than one resolution-1 base cell


def test_a_wholly_contained_region_costs_far_less_than_its_cells(
    parcel: Polygon,
) -> None:
    """The index of a fill is a fraction of the size of its cell list.

    The measured ratio is what sets :data:`MAX_FILL_CELLS`: the composed
    index runs about three bytes per cell, while naming every cell as an
    atomic string runs about sixty. A fill that materialised the atomic
    form before composing it would exhaust memory two orders of
    magnitude earlier.
    """
    filled = itacart.polyfill(parcel, 7)
    cells = itacart.decompose(filled)
    atomic_bytes = sum(len(cell) for cell in cells)
    assert len(filled) * 10 < atomic_bytes


# --------------------------------------------------------------------------
# Sweep: measurements this phase owes the ones that come after
# --------------------------------------------------------------------------


def test_polyfill_never_consults_the_topology_caches(parcel: Polygon) -> None:
    """Filling asks topology nothing, at any resolution (``P-6.16``).

    The pendency asked whether the 4,096-entry contact cache is large
    enough under a resolution-13 fill. Measured, the question is empty:
    the descent subdivides in the sheared lattice and never asks for a
    neighbour, so both caches stay at zero hits, zero misses and zero
    entries.

    That is also why this phase never needed ``P-6.20`` resolved. The
    ``get_neighbor`` asymmetry cannot reach a fill that composes no step.
    """
    from itacart import topology

    caches = [
        value
        for value in vars(topology).values()
        if hasattr(value, "cache_info") and hasattr(value, "cache_clear")
    ]
    assert caches, "topology has no caches to check; the test has gone stale"
    for cache in caches:
        cache.cache_clear()
    itacart.polyfill(parcel, 7)
    itacart.count_internal_cells(parcel, 9)
    for cache in caches:
        info = cache.cache_info()
        assert (info.hits, info.misses, info.currsize) == (0, 0, 0), info


def test_children_union_matches_the_parent_only_away_from_the_border() -> None:
    """The union of children agrees with the parent in the interior only.

    ``P-4.11`` reported a 1.17e-3 excess without saying where. Enumerated
    by family rather than sampled, the four separate cleanly: an interior
    parent sits at about -1.5e-7, a deficit rather than an excess,
    because each child approximates the curved parallel by the chord of
    its own height and the chords fall inside the parent's. The
    meridian column and the border-absorbing column run positive and
    grow toward the pole.

    The relevance to criterion 2 is that :func:`polyfill` refuses all
    three non-interior families, so the containment chain is asserted
    over squares that do carry this agreement.
    """
    from shapely.ops import unary_union

    from itacart import hierarchy

    side = itacart.cell_size(1)
    for row in (0, 10, 100, 300, 500):
        last = itacart.last_lattice_column("NE", row, side)
        interior = f"NE({last // 2:04d}/{row:04d})"
        assert not itacart.absorbs_border(interior)
        parent = itacart.cell_to_polygon(interior)
        children = list(itacart.get_children(interior))[0]
        union = unary_union([itacart.cell_to_polygon(k) for k in children])
        relative = (union.area - parent.area) / parent.area
        assert -1e-5 < relative <= 0.0, (interior, relative)
        assert len(children) == itacart.refinement_ratio(2)
    assert hierarchy is not None


def test_border_children_survive_every_refinement() -> None:
    """``get_children`` refines a border cell down to resolution 13.

    ``P-5.2`` asked whether the overlap filter behind the validity gate
    is redundant. It is not, but the threshold it compared against was a
    fixed square metre, which is the whole nominal area of a
    resolution-9 cell and larger than every cell below it. The
    comparison is strict, so a child wholly inside its parent was
    rejected for being the size it is supposed to be, and the eastern
    border of the grid had no refinement below resolution 8 at all. The
    threshold is now a fraction of the nominal area of the level and
    shrinks with the cell.

    Two assertions, and the second is the guard against overcorrection.
    Every level has to yield children, and the counts have to keep
    deviating from the canonical refinement ratio where the geometry
    says they should: a threshold pushed too far down would admit
    candidates the parent does not hold and iron that deviation flat.
    The measured series for row 100, resolutions 1 to 12, is
    ``2, 21, 4, 24, 4, 23, 4, 25, 4, 25, 4, 24``. Its first seven
    entries are exactly the ones measured before the fix, and the eighth
    used to be zero.
    """
    side = itacart.cell_size(1)
    row = 100
    cell = f"NE({itacart.last_lattice_column('NE', row, side):04d}/{row:04d})"
    assert itacart.absorbs_border(cell)
    counts: dict[int, int] = {}
    while itacart.get_resolution(cell) < 13:
        resolution = itacart.get_resolution(cell)
        children = list(itacart.get_children(cell))[0]
        counts[resolution] = len(children)
        absorbing = [k for k in children if itacart.absorbs_border(k)]
        assert absorbing, (resolution, children)
        cell = absorbing[-1]
    assert sorted(counts) == list(range(1, 13)), counts
    assert [counts[r] for r in sorted(counts)] == [
        2,
        21,
        4,
        24,
        4,
        23,
        4,
        25,
        4,
        25,
        4,
        24,
    ], counts
    assert all(count > 0 for count in counts.values()), counts
    deviating = [r for r in counts if counts[r] != itacart.refinement_ratio(r + 1)]
    assert deviating, counts


# --------------------------------------------------------------------------
# Filling inside an extension zone
# --------------------------------------------------------------------------


_ZONE_FOOTPRINTS: tuple[tuple[str, float, float], ...] = (
    ("Fiji, southern hemisphere", -18.0, -17.8),
    ("Chukotka, northern hemisphere", 65.0, 65.2),
)
"""One latitude band inside each extension zone, as ``(name, low, high)``.

Fiji is realised on resolution-1 rows 171 to 237 and reaches 182 degrees;
Chukotka on rows 709 to 799, reaching 190.5. Both bands sit well inside
their zone, so a footprint written past 180 there is inside the domain
and has to fill rather than be refused.
"""


@pytest.mark.parametrize("name, low, high", _ZONE_FOOTPRINTS)
def test_a_footprint_astride_the_antemeridian_fills_inside_a_zone(
    name: str, low: float, high: float
) -> None:
    """The positive case criterion 9 of the previous phase never had.

    A ring written from 179.9 to 180.3 lies on one side of the line
    already, so it does not cross: inside an extension zone the domain
    reaches past 180 and the whole footprint is addressable. It used to
    raise a bare ``GEOSException`` from the geometry engine instead,
    because densification normalised its interior vertices to
    ``(-180, 180]`` and folded the ring back across the globe.
    """
    footprint = Polygon([(179.9, low), (180.3, low), (180.3, high), (179.9, high)])
    assert not itacart.crosses_antemeridian(footprint)
    cells = itacart.decompose(itacart.polyfill(footprint, 3))
    assert cells, name
    assert all(itacart.get_resolution(cell) == 3 for cell in cells)


@pytest.mark.parametrize("name, low, high", _ZONE_FOOTPRINTS)
def test_a_footprint_wholly_past_180_fills_inside_a_zone(
    name: str, low: float, high: float
) -> None:
    """Written entirely past the line, and still inside the domain.

    The strict reading of the antemeridian test takes the longitude
    literally, so this footprint does not cross either. It broke for the
    same reason and is fixed by the same change.
    """
    footprint = Polygon([(180.2, low), (180.4, low), (180.4, high), (180.2, high)])
    assert not itacart.crosses_antemeridian(footprint)
    assert itacart.decompose(itacart.polyfill(footprint, 3)), name


@pytest.mark.parametrize("name, low, high", _ZONE_FOOTPRINTS)
def test_the_wrapped_spelling_of_the_same_footprint_is_still_refused(
    name: str, low: float, high: float
) -> None:
    """The other direction of travel, which is a crossing and stays one.

    Spelled with a jump from 179.9 to -179.7 the ring does cross, and the
    refusal is the package's own exception. Filling one spelling must not
    quietly start filling the other: they describe different rings.
    """
    wrapped = Polygon([(179.9, low), (-179.7, low), (-179.7, high), (179.9, high)])
    assert itacart.crosses_antemeridian(wrapped)
    with pytest.raises(AntemeridianError):
        itacart.polyfill(wrapped, 3)


def test_a_footprint_past_180_outside_every_zone_is_refused_by_name() -> None:
    """Outside a zone the domain stops at the line, and so does the fill.

    What matters here is not that it is refused but that the refusal is
    the package's, not the geometry engine's. A raw ``GEOSException``
    tells the caller nothing they can act on.
    """
    footprint = Polygon([(179.9, 5.0), (180.3, 5.0), (180.3, 5.2), (179.9, 5.2)])
    assert not itacart.crosses_antemeridian(footprint)
    with pytest.raises(NonExistentCellError):
        itacart.polyfill(footprint, 3)


def test_densification_keeps_the_longitude_branch_it_was_given() -> None:
    """The mechanism behind the four tests above, isolated.

    ``direct_geodesic`` normalises to ``(-180, 180]``, which is correct
    for the question it answers and wrong for densification: the
    interior vertices of a segment crossing the line came back on the
    other branch and the ring self-intersected. Measured on the ring
    itself rather than through the fill, so that a future change to the
    fill cannot hide a regression here.
    """
    ring = Polygon([(179.9, -18.0), (180.3, -18.0), (180.3, -17.8), (179.9, -17.8)])
    dense = itacart.densify_orthodromic(ring, 1000.0)
    longitudes = [point[0] for point in dense.exterior.coords]
    assert min(longitudes) > 179.0, min(longitudes)
    assert max(longitudes) <= 180.3 + 1e-9, max(longitudes)
    assert dense.is_valid


# --------------------------------------------------------------------------
# The radius of the boundary screen
# --------------------------------------------------------------------------


def _parcel_east_of_greenwich(degrees: float, metres: float = 60.0) -> Polygon:
    """A small square parcel that many degrees east of the prime meridian."""
    side = metres / 111_319.0
    lat = 51.4776
    return Polygon(
        [
            (degrees, lat),
            (degrees + side, lat),
            (degrees + side, lat + side),
            (degrees, lat + side),
        ]
    )


def test_the_screen_refuses_at_the_resolution_that_was_asked_for() -> None:
    """A parcel beside the meridian is not refused for its neighbour's sake.

    The screen used to run at resolution 1, so it refused a band up to
    ten kilometres wide -- the whole base cell -- to protect against a
    family one cell wide. At resolution 7 that is a thousand times too
    much, and at resolution 13 a million times.

    The band now scales with the target, so the same parcel is refused
    at a coarse resolution and filled at a fine one. Nothing about
    ``D-7.3`` changes: the meridian cells are refused either way.
    """
    parcel = _parcel_east_of_greenwich(0.01)
    with pytest.raises(NonExistentCellError):
        itacart.polyfill(parcel, 3)
    cells = itacart.decompose(itacart.polyfill(parcel, 7))
    assert cells
    assert all(itacart.get_resolution(cell) == 7 for cell in cells)


def test_what_is_filled_beside_the_meridian_is_ordinary_and_nominal() -> None:
    """The cells the narrower screen lets through are real cells.

    Inside an anomalous base cell the index descent and the square
    descent are not the same tree, so admitting geometry there is only
    safe if the indices the fill emits name the squares it accepted.
    Measured on the round trip rather than assumed.
    """
    parcel = _parcel_east_of_greenwich(0.001)
    cells = itacart.decompose(itacart.polyfill(parcel, 7))
    assert cells
    nominal = itacart.nominal_cell_area(7)
    for cell in cells:
        assert itacart.is_valid_index(cell), cell
        assert not itacart.absorbs_border(cell), cell
        assert itacart.cell_shape(cell) == "parallelogram", cell
        area = Polygon(boundary.plane_ring(cell)[1]).area
        assert area == pytest.approx(nominal, rel=1e-9), cell
        assert itacart.cell_to_polygon(cell).intersects(parcel), cell


def test_the_three_families_are_still_refused_when_actually_touched() -> None:
    """Narrowing the radius must not narrow the refusal itself.

    ``D-7.3`` is untouched by the change of radius: geometry that reaches
    a cell of one of the three families is refused exactly as before.
    """
    with pytest.raises(NonExistentCellError):
        itacart.polyfill(_parcel_east_of_greenwich(0.0), 13)
    side = itacart.cell_size(1)
    for quadrant in QUADRANTS:
        polar = max(
            row
            for row in range(700, 1000)
            if itacart.last_lattice_column(quadrant, row, side) > 0
        )
        for column, row in (
            (0, 100),
            (itacart.last_lattice_column(quadrant, 100, side), 100),
            (1, polar),
        ):
            with pytest.raises((NonExistentCellError, DomainError)):
                geometry._check_addressable(quadrant, column, row)


def _screen_refuses(quadrant: str, u0: float, v0: float, side: float) -> bool:
    """The geometric screen, expressed only on the node's own coordinates.

    The lattice column of a node of side ``s`` at ``(u0, v0)`` is
    ``(u0 - v0) / s``, which reduces to ``u_index - row`` at resolution
    1. Nothing here names a cell, which is the point: the fill may not
    reach an index-rendering door, so it cannot ask ``absorbs_border``.
    """
    column = round((u0 - v0) / side)
    row = round(v0 / side)
    last_here = itacart.last_lattice_column(quadrant, row, side)
    beyond = itacart.last_lattice_column(quadrant, row + 1, side)
    return column <= 0 or last_here <= 0 or beyond <= 0 or column >= last_here


def _is_anomalous(cell: str) -> bool:
    """The index-side truth the screen has to agree with."""
    try:
        if itacart.absorbs_border(cell):
            return True
        return itacart.cell_shape(cell) != "parallelogram"
    except Exception:  # pragma: no cover - defensive
        return True


def _children_in_lockstep(
    cell: str, u0: float, v0: float, side: float, level: int
) -> list[tuple[str, float, float, float, int]]:
    """Every child of a node, as index and as square, paired by position."""
    step = level + 1
    divisor = itacart.linear_refinement_ratio(step)
    child = side / divisor
    return [
        (
            hierarchy._descend(cell, geometry._child_code(row, column, step)),
            u0 + column * child,
            v0 + row * child,
            child,
            step,
        )
        for row in range(divisor)
        for column in range(divisor)
    ]


def test_the_geometric_screen_never_admits_an_anomalous_cell() -> None:
    """The screen is validated against the index side, by enumeration.

    Complete to level 3, over the three families and the ordinary
    interior in four quadrants. Every node is compared both ways, and
    the geometric answer is never allowed to be the more permissive of
    the two.

    Divergences do exist and are all in the safe direction. Inside the
    meridian family they sit at a negative column -- west of the line,
    outside the quadrant, unreachable because the quadrant clip leaves
    ``u - v`` positive -- and that much is asserted. Inside the
    border-absorbing family they sit at positive columns, because
    absorption gives the parent more children than the square
    subdivision has and the two are not the same tree there. Only
    conservativeness is claimed for that family; equivalence is not.

    Three levels is not thirteen, and the depth is the whole point of
    the screen, so the shallow sweep is the floor and not the argument.
    ``test_the_screen_agrees_all_the_way_down_the_frontier`` carries it
    to resolution 13.
    """
    side = itacart.cell_size(1)
    checked = 0
    for quadrant in QUADRANTS:
        for row in (0, 1, 100, 451, 700):
            last = itacart.last_lattice_column(quadrant, row, side)
            for column in sorted({0, 1, 2, last // 2, last - 1, last, last + 1}):
                stack = [
                    (
                        f"{quadrant}({column:04d}/{row:04d})",
                        (column + row) * side,
                        row * side,
                        side,
                        1,
                    )
                ]
                while stack:
                    cell, u0, v0, cell_side, level = stack.pop()
                    refused = _screen_refuses(quadrant, u0, v0, cell_side)
                    anomalous = _is_anomalous(cell)
                    checked += 1
                    if anomalous and not refused:
                        raise AssertionError(f"screen admits {cell}")
                    if refused and not anomalous and column == 0:
                        assert round((u0 - v0) / cell_side) < 0, cell
                    if level < 3:
                        stack.extend(
                            _children_in_lockstep(cell, u0, v0, cell_side, level)
                        )
    assert checked > 14_000, checked


def test_the_screen_agrees_all_the_way_down_the_frontier() -> None:
    """The same agreement, carried from resolution 1 to resolution 13.

    A complete sweep to resolution 13 is not enumerable -- a single base
    cell holds a million cells per axis down there -- so this walks the
    frontier instead, which is where the derivation says the only
    anomalies can be. At every step *all* children of the current node
    are compared both ways, and the walk then continues into an
    anomalous one. Ordinary nodes are covered exhaustively by the
    shallow sweep; what this adds is depth on the one path where the
    families actually live.

    The shape of the shortcut is deliberate. The refinement is a fixed
    point only where the enumeration is complete, and the enumeration is
    complete at every level of the walk -- it is the breadth across
    levels that is traded away, not the breadth within one.

    This test exists because a bound measured at two levels and asserted
    for thirteen is exactly how the fixed overlap threshold survived
    into this phase.
    """
    base = itacart.cell_size(1)
    deepest = 0
    checked = 0
    for quadrant in QUADRANTS:
        last = itacart.last_lattice_column(quadrant, 451, base)
        polar = max(
            row
            for row in range(700, 1000)
            if itacart.last_lattice_column(quadrant, row, base) > 0
        )
        for column, row in ((0, 451), (last, 451), (1, polar)):
            cell = f"{quadrant}({column:04d}/{row:04d})"
            u0, v0, side, level = (column + row) * base, row * base, base, 1
            while level < 13:
                children = _children_in_lockstep(cell, u0, v0, side, level)
                anomalous_children = []
                for name, cu, cv, cs, step in children:
                    refused = _screen_refuses(quadrant, cu, cv, cs)
                    anomalous = _is_anomalous(name)
                    checked += 1
                    if anomalous and not refused:
                        raise AssertionError(f"screen admits {name}")
                    if anomalous:
                        anomalous_children.append((name, cu, cv, cs, step))
                if not anomalous_children:
                    break
                cell, u0, v0, side, level = anomalous_children[-1]
                deepest = max(deepest, level)
            assert level == 13, (quadrant, column, row, level)
    assert deepest == 13, deepest
    assert checked > 1_500, checked


def test_the_band_of_a_polar_base_cell_is_the_whole_cell() -> None:
    """At the pole the anomaly is not a strip, and the band says so.

    ``_anomalous_band`` draws the three families as bands one cell side
    wide, which is the point of narrowing the screen. The polar family
    is the one case where that narrowing buys nothing: every target row
    in the last base row of a quadrant is clipped by the pole, so the
    band covers the base cell entirely and the screen refuses exactly
    what it refused before.

    Two rows and three resolutions, because the branch that detects the
    polar reach differs between them. Row 999 is the last row that
    addresses a cell at resolution 1; row 1000 is straddled by the pole,
    so it addresses none, and its topmost target rows lie beyond the
    meridian quadrant.
    """
    base = itacart.cell_size(1)
    for row in (999, 1000):
        for resolution in (1, 3, 5):
            band = geometry._anomalous_band("NE", 1, row, itacart.cell_size(resolution))
            assert not band.is_empty, (row, resolution)
            if row == 1000:
                assert band.area == pytest.approx(base * base, rel=1e-9), (
                    row,
                    resolution,
                )
            else:
                assert band.area > 0.7 * base * base, (row, resolution)


def test_the_polar_reach_is_found_by_search_when_it_starts_mid_cell() -> None:
    """The pole can fall inside a base cell, not only below it.

    Row 1000 at resolution 5 is the case: its lowest target rows still
    address cells, its highest do not, and the first row that does not
    has to be located rather than assumed. A linear walk would visit up
    to a million rows at resolution 13, so it is a search, and the
    search is what this covers.
    """
    base = itacart.cell_size(1)
    side = itacart.cell_size(5)
    lowest = int(math.floor(1000 * base / side))
    highest = int(math.floor((1000 * base + base - side / 2.0) / side))
    assert itacart.last_lattice_column("NE", lowest + 1, side) > 0
    assert itacart.last_lattice_column("NE", highest + 1, side) <= 0
    crossings = [
        row
        for row in range(lowest, highest + 1)
        if (itacart.last_lattice_column("NE", row + 1, side) <= 0)
        != (itacart.last_lattice_column("NE", row, side) <= 0)
    ]
    assert len(crossings) == 1, crossings
    band = geometry._anomalous_band("NE", 1, 1000, side)
    assert band.area == pytest.approx(base * base, rel=1e-9)


def test_the_polar_row_narrows_with_the_target_like_the_others() -> None:
    """The polar family narrows too, and the first version of this lied.

    It asserted that a parcel in the last base row is refused at every
    resolution, on the reasoning that the pole clips the whole row. That
    is true of base row 1000, which the pole straddles, and false of row
    999: most of its target rows address cells perfectly well, and only
    the topmost do not. The screen refuses the top and admits the rest,
    which is the whole intent of measuring the band at the target side.

    So the assertion is the one that survives measurement. Base row 999
    is refused entirely at resolution 1, because at that side the row is
    the polar row and holds one cell. Below that side the band shrinks,
    and it never grows back.
    """
    base = itacart.cell_size(1)
    whole = base * base
    assert geometry._anomalous_band("NE", 1, 999, base).area == pytest.approx(
        whole, rel=1e-9
    )
    areas = [
        geometry._anomalous_band("NE", 1, 999, itacart.cell_size(r)).area
        for r in (1, 3, 5, 7)
    ]
    assert areas[0] == pytest.approx(whole, rel=1e-9)
    assert all(area < whole for area in areas[1:]), areas
    assert max(areas[1:]) < 0.9 * whole, areas
    for row in (999, 1000):
        for resolution in (1, 3, 5):
            band = geometry._anomalous_band("NE", 1, row, itacart.cell_size(resolution))
            assert not band.is_empty, (row, resolution)
