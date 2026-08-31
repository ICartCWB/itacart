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

import math

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
from itacart import geometry
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
    """``compact=True`` compacts less than ``compact_cells`` does.

    The fill compacts a node when the plane square is contained in the
    projected geometry, which is a geometric test and sensitive to the
    geodetic round trip that built the footprint. ``compact_cells`` asks
    a lexical question -- are all the children present -- and answers yes
    where the geometric test answers no. Measured on Figure 7(a): the
    published index and the lexical compaction both hold 114 cells, the
    geometric one holds 138, and they differ at exactly one node.

    The ground covered is the same. This pins the divergence rather than
    hiding it; which of the two ``compact=True`` should mean is an open
    question for the bridge, not something this test settles.
    """
    from shapely.ops import unary_union

    cells = itacart.decompose(FIGURE_7A)
    footprint = unary_union([itacart.cell_to_polygon(cell) for cell in cells])
    geometric = itacart.polyfill(footprint, 7, compact=True)
    lexical = itacart.compact_cells(itacart.polyfill(footprint, 7))
    assert itacart.count_cells(lexical) < itacart.count_cells(geometric)
    assert set(itacart.uncompact_cells(geometric, 7)) == set(
        itacart.uncompact_cells(lexical, 7)
    )


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


def _residual_within_derived_bound(parcel: Polygon, resolution: int) -> bool:
    """Whether the count's area error respects the proved bound.

    The proof, in two lines. A cell whose closure does not meet the
    outline is decided exactly, so only cells meeting the outline can be
    misclassified. Every such cell lies inside the ``s*sqrt(2)``
    neighbourhood of the outline, whose area is at most ``2*r*P +
    pi*r*r``. Substituting ``r = s*sqrt(2)`` gives

        |n*a - A| <= 2*sqrt(2)*P*s + 2*pi*s*s

    which holds for any rectifiable outline, convex or not, holed or not.
    """
    plane = geometry._project(parcel)
    side = itacart.cell_size(resolution)
    perimeter = _plane_perimeter(plane)
    bound = 2.0 * math.sqrt(2.0) * perimeter * side + 2.0 * math.pi * side * side
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

    The derived bound carries a factor of ``2*sqrt(2)`` and an ``s*s``
    term. ``P*s`` is about 2.8 times tighter and every parcel measured
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


INDEX_RENDERING_DOORS: tuple[tuple[str, str], ...] = (
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
"""Every route in the package from components to an index string.

Enumerated rather than sampled, and the enumeration is itself the
measurement: a proof that counting names no cell is only as good as the
list of ways a cell could be named. Four of these are the duplicated
renderers the package already tracks; the rest reach a string by other
paths.

Two of them, the contact kind and the quadrant name, do not build an
index at all. They are armed anyway. The list is derived from a
mechanical rule -- every function in these modules annotated to return
``str`` -- rather than curated by judgement about which ones "really"
render, because a curated list is exactly what goes stale without
anyone noticing.

Hand-written first, and the completeness check below found seven more.
"""


def test_the_door_list_is_complete() -> None:
    """Every enumerated door exists, and none was missed.

    The second half is the part that matters. A list of doors that has
    fallen behind the code makes the guard below prove less than it
    claims while still passing, which is the failure mode this whole
    test pair exists to avoid.
    """
    import importlib
    import inspect

    for module_name, attribute in INDEX_RENDERING_DOORS:
        module = importlib.import_module(module_name)
        assert hasattr(module, attribute), f"{module_name}.{attribute} is gone"

    listed = set(INDEX_RENDERING_DOORS)
    found: set[tuple[str, str]] = set()
    for module_name in {name for name, _ in INDEX_RENDERING_DOORS}:
        module = importlib.import_module(module_name)
        for name, value in vars(module).items():
            if not callable(value) or isinstance(value, type):
                continue
            if getattr(value, "__module__", None) != module_name:
                continue
            try:
                hints = inspect.get_annotations(value, eval_str=False)
            except (TypeError, ValueError):  # pragma: no cover - defensive
                continue
            if hints.get("return") in ("str", str):
                found.add((module_name, name))
    missed = found - listed
    assert not missed, f"unlisted routes to an index string: {sorted(missed)}"


def test_count_never_builds_an_index_fragment(
    parcel: Polygon, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Counting names no cell, proved against every route to a string.

    Patching the descent's own helpers would prove nothing: the property
    is that no *package* route to an index string is reached, and a
    descent that rendered through another module would satisfy a guard
    aimed at itself. So all thirteen enumerated doors are replaced with
    detonators at once, and the count has to survive with every one of
    them armed.

    Measured with the doors armed: 33,258 cells at resolution 7, none of
    the doors touched.
    """
    import importlib

    def _detonate(module_name: str, attribute: str):
        def fired(*args: object, **kwargs: object) -> str:
            raise AssertionError(f"the count reached {module_name}.{attribute}")

        return fired

    for module_name, attribute in INDEX_RENDERING_DOORS:
        module = importlib.import_module(module_name)
        monkeypatch.setattr(module, attribute, _detonate(module_name, attribute))

    assert itacart.count_internal_cells(parcel, 7) == 33_258


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


def test_border_children_vanish_below_a_fixed_area_threshold() -> None:
    """``get_children`` returns nothing for a border cell from resolution 8.

    ``P-5.2`` asked whether the overlap filter behind the validity gate
    is redundant. It is not, but what it rejects is not a non-overlap:
    ``hierarchy._OVERLAP_EPSILON_M2`` is a fixed 1 square metre, while a
    resolution-9 cell has a nominal area of exactly 1 square metre and
    everything finer is smaller. The comparison is strict, so a child
    wholly inside its parent is rejected for being the size it is
    supposed to be.

    Measured across twelve rows and every level: zero overlap rejections
    up to level 8, then 229 at level 9 and 39 at level 10 with none
    accepted. The F5 measurement saw zero because it went two levels
    deep.

    This test pins the defect rather than the intended behaviour, and
    fails the day it is fixed, which is the point. ``hierarchy.py`` is
    not this phase's to edit.
    """
    side = itacart.cell_size(1)
    row = 100
    cell = f"NE({itacart.last_lattice_column('NE', row, side):04d}/{row:04d})"
    assert itacart.absorbs_border(cell)
    counts = {}
    while itacart.get_resolution(cell) < 13:
        resolution = itacart.get_resolution(cell)
        children = list(itacart.get_children(cell))[0]
        counts[resolution] = len(children)
        absorbing = [k for k in children if itacart.absorbs_border(k)]
        if not absorbing:
            break
        cell = absorbing[-1]
    assert counts[8] == 0, counts
    assert all(counts[r] > 0 for r in counts if r < 8), counts
