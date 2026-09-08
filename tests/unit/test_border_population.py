"""The border population, enumerated in four quadrants and versioned here.

Two populations are described together because they are the same object
seen through two doors. One is the set of cells whose child count leaves
the refinement ratio; the other is the band of columns the filling
refuses. Both live on the outer edge of a row, and a phase that treats
them as separate problems fixes one and leaves the other.

Every number below came from a run over the whole resolution-1 grid --
5 100 221 cells, four quadrants, enumerated rather than sampled -- and
what this module re-runs is the part that fits in a suite. Where a claim
rests on the out-of-band sweep rather than on what runs here, the
docstring says so. The distinction matters: a number nobody can re-derive
is a number that rots, and this project has twice shipped a count that
was validated in one quadrant and wrong in the other three.
"""

from __future__ import annotations

import collections
import functools

import pytest
from shapely.geometry import Polygon

import itacart
from itacart.constants import CELL_SIZE_M, EXTENSION_ZONES
from itacart.exceptions import ITACaRTError

SIDE = CELL_SIZE_M[1]
assert SIDE is not None

QUADRANTS = ("NE", "NW", "SE", "SW")

#: Resolution 1 refines four-to-one. Any other count is the border.
CANONICAL_RATIO = 4

#: Column zero is the meridian column and the meridian belongs to the
#: east, so the western quadrants start at one. Counting from zero in all
#: four is how a previous enumeration reported one cell too many per
#: western row.
FIRST_COLUMN = {"NE": 0, "SE": 0, "NW": 1, "SW": 1}


def _first(quadrant: str) -> int:
    return FIRST_COLUMN[quadrant]


def _last_column_family() -> list[str]:
    """Every outermost cell of every row, in four quadrants.

    Rows whose last column falls below the first hold no cell at all,
    which is how the western quadrants come to have no polar cap: the cap
    sits in column zero, and they have no column zero.
    """
    family = []
    for quadrant in QUADRANTS:
        first = _first(quadrant)
        for row in range(0, 1001):
            last = itacart.last_lattice_column(quadrant, row, SIDE)
            if last >= first:
                family.append(f"{quadrant}({last:04d}/{row:04d})")
    return family


def _child_count(cell: str) -> int:
    # ``flatten`` matters. Without it the call yields one list per atom,
    # so counting what comes out counts the containers and reports one
    # child for every cell in the grid.
    return sum(1 for _ in itacart.get_children(cell, flatten=True))


def _fill_verdict(cell: str) -> str | None:
    """``None`` if the filling accepts the cell, else the exception name."""
    try:
        itacart.count_internal_cells(Polygon(itacart.cell_to_boundary(cell)), 2)
    except Exception as exc:
        return type(exc).__name__
    return None


@functools.lru_cache(maxsize=None)
def _refused_band(quadrant: str, row: int) -> tuple[tuple[int, str], ...]:
    """Walk inward from the last column until the filling accepts.

    The walk assumes the refusals are a contiguous suffix of the row.
    That assumption is not free and is checked below on whole rows,
    including both rows where an extension zone begins.
    """
    first = _first(quadrant)
    last = itacart.last_lattice_column(quadrant, row, SIDE)
    refused: list[tuple[int, str]] = []
    for column in range(last, first - 1, -1):
        verdict = _fill_verdict(f"{quadrant}({column:04d}/{row:04d})")
        if verdict is None:
            break
        refused.append((column, verdict))
    return tuple(refused)


@functools.lru_cache(maxsize=None)
def _band_profile() -> tuple[tuple[str, int, tuple[tuple[int, str], ...]], ...]:
    """The whole profile, walked once and shared by the tests that read it.

    Three tests ask different questions of the same walk. Walking it three
    times costs three times as much under coverage instrumentation, which
    traces every line of the package the walk touches, and the suite pays
    that cost on every run for the rest of the project.
    """
    return tuple(
        (quadrant, row, _refused_band(quadrant, row))
        for quadrant in QUADRANTS
        for row in range(0, 1001)
        if itacart.last_lattice_column(quadrant, row, SIDE) >= _first(quadrant)
    )


# --------------------------------------------------------------------------
# The grid the population is measured against
# --------------------------------------------------------------------------


def test_the_resolution_1_grid_holds_the_established_cell_count() -> None:
    """The enumerator answers the count the project already agreed on.

    An instrument is validated against a known answer before it measures
    an unknown one. This one is validated against 5 100 221, the count
    that replaced 5 102 223 once the enumerator was run in four quadrants
    instead of one, and the per-quadrant split is asserted too because
    the totals can agree while the halves do not.
    """
    counts = {}
    for quadrant in QUADRANTS:
        first = _first(quadrant)
        total = 0
        for row in range(0, 1001):
            last = itacart.last_lattice_column(quadrant, row, SIDE)
            if last >= first:
                total += last - first + 1
        counts[quadrant] = total

    assert counts == {
        "NE": 1279545,
        "NW": 1270568,
        "SE": 1276967,
        "SW": 1273141,
    }
    assert sum(counts.values()) == 5100221


def test_the_quadrants_are_not_congruent_and_the_extensions_say_why() -> None:
    """The eastern quadrants are larger, by exactly the two extensions.

    Read without this, the asymmetry above looks like the western
    off-by-one that has bitten this project twice. It is not. Each
    extension zone hands columns from a western quadrant to its eastern
    partner across the antemeridian, so every quadrant has two
    discontinuities and not one -- where the zone opens and where it
    closes -- and the pair moves in opposite directions on the two sides
    of the meridian. Looking only for the row where the count *drops*
    finds one of each pair and reads as one defect per quadrant.
    """
    lasts = {
        quadrant: [itacart.last_lattice_column(quadrant, r, SIDE) for r in range(1001)]
        for quadrant in QUADRANTS
    }
    jumps = {
        quadrant: [r for r in range(1000) if abs(row[r] - row[r + 1]) > 6]
        for quadrant, row in lasts.items()
    }
    assert jumps == {
        "NE": [708, 799],
        "NW": [708, 799],
        "SE": [170, 237],
        "SW": [170, 237],
    }

    # East gains where west loses, at the same row, in both zones.
    for east, west, opens, closes in (("NE", "NW", 708, 799), ("SE", "SW", 170, 237)):
        assert lasts[east][opens + 1] > lasts[east][opens]
        assert lasts[west][opens + 1] < lasts[west][opens]
        assert lasts[east][closes + 1] < lasts[east][closes]
        assert lasts[west][closes + 1] > lasts[west][closes]

    # And each discontinuity sits on a declared latitude limit of a zone.
    limits = {abs(zone.lat_min) for zone in EXTENSION_ZONES.values()} | {
        abs(zone.lat_max) for zone in EXTENSION_ZONES.values()
    }
    for quadrant, rows in jumps.items():
        for row in rows:
            ring = itacart.cell_to_boundary(
                f"{quadrant}"
                f"({itacart.last_lattice_column(quadrant, row, SIDE):04d}"
                f"/{row:04d})"
            )
            latitude = abs(sum(point[1] for point in ring) / len(ring))
            assert min(abs(latitude - limit) for limit in limits) < 0.3


# --------------------------------------------------------------------------
# The population that leaves the ratio
# --------------------------------------------------------------------------


@pytest.mark.slow
def test_the_out_of_ratio_population_is_inside_the_last_column_family() -> None:
    """Enumerated over the whole family, not over four rows of it.

    The family is every outermost cell of every row in four quadrants.
    Three of its properties are asserted rather than described: how many
    of its members leave the ratio, what counts they take, and what
    shapes they are. The two triangles are the polar caps, which the
    western quadrants do not have.
    """
    family = _last_column_family()
    assert len(family) == 4002

    counts = {cell: _child_count(cell) for cell in family}
    off = {cell: n for cell, n in counts.items() if n != CANONICAL_RATIO}

    assert len(off) == 2852
    assert len(family) - len(off) == 1150
    assert collections.Counter(n for n in off.values()) == {
        1: 2,
        2: 423,
        3: 752,
        5: 836,
        6: 839,
    }
    assert collections.Counter(cell[:2] for cell in off) == {
        "NE": 710,
        "NW": 720,
        "SE": 708,
        "SW": 714,
    }

    shapes = collections.Counter(str(itacart.cell_shape(cell)) for cell in off)
    assert shapes == {"trapezoid": 2850, "triangle": 2}
    assert {c for c in off if str(itacart.cell_shape(c)) == "triangle"} == {
        "NE(0000/1000)",
        "SE(0000/1000)",
    }
    assert all(itacart.absorbs_border(cell) for cell in off)


@pytest.mark.slow
def test_no_cell_away_from_the_last_column_leaves_the_ratio() -> None:
    """The converse half, which is what makes the claim a claim.

    The full sweep behind this ran over all 5 100 221 cells and found
    2 852 out of the ratio, every one of them the last column of its row.
    That sweep takes a quarter of an hour and does not belong in a suite,
    so what runs here is thirteen whole rows in four quadrants -- whole,
    because the last column of a row is one draw in two thousand and
    sampling is the surest way to miss it.
    """
    rows = (0, 1, 50, 100, 237, 300, 450, 500, 708, 799, 900, 990, 998)
    off = []
    seen = 0
    for quadrant in QUADRANTS:
        first = _first(quadrant)
        for row in rows:
            last = itacart.last_lattice_column(quadrant, row, SIDE)
            for column in range(first, last):  # last column excluded
                seen += 1
                if _child_count(f"{quadrant}({column:04d}/{row:04d})") != (
                    CANONICAL_RATIO
                ):
                    off.append(f"{quadrant}({column:04d}/{row:04d})")
    assert seen == 65627
    assert off == []


def test_the_ratio_pendency_measured_the_border_family_not_every_parent() -> None:
    """Why 291 in 1 000 and 99.9 per cent are both right.

    The pendency reports roughly three parents in ten staying on the
    ratio, and a sweep of whole rows reports one in a thousand leaving
    it. Neither is wrong and neither reproduces the other, because they
    count different denominators: the pendency's is the border family,
    where the fraction on the ratio is 1 150 of 4 002, and the sweep's is
    every parent in the row.

    The reconciliation is asserted rather than narrated so that it cannot
    be re-litigated from memory.
    """
    on_ratio_in_family = 1150 / 4002
    off_ratio_in_grid = 2852 / 5100221

    assert 0.28 < on_ratio_in_family < 0.30  # the pendency's 291/1000
    assert off_ratio_in_grid < 0.001  # the sweep's tenth of a per cent


# --------------------------------------------------------------------------
# The band the filling refuses
# --------------------------------------------------------------------------


@pytest.mark.slow
def test_the_refused_band_is_a_contiguous_suffix_of_its_row() -> None:
    """The assumption the inward walk rests on, checked on whole rows.

    Seven rows scanned column by column, chosen to include both rows
    where an extension zone begins, a polar row, and the shortest row
    that still holds cells. If a refusal ever appears in the middle of a
    row, the walk below stops measuring what it claims to measure, and
    this is where that shows up.
    """
    for quadrant, row in (
        ("NE", 0),
        ("NE", 500),
        ("NW", 708),
        ("NW", 900),
        ("SE", 237),
        ("SE", 999),
        ("SW", 998),
    ):
        first = _first(quadrant)
        last = itacart.last_lattice_column(quadrant, row, SIDE)
        refused = [
            column
            for column in range(first, last + 1)
            if _fill_verdict(f"{quadrant}({column:04d}/{row:04d})") is not None
        ]
        assert refused == list(
            range(last - len(refused) + 1, last + 1)
        ), f"{quadrant} row {row}: refusals are not a suffix"


@pytest.mark.slow
def test_the_refused_band_profile_over_every_row_of_every_quadrant() -> None:
    """The band, row by row, in four quadrants, pinned by equality.

    Equality rather than a ceiling, because a ceiling lets the number rot
    in both directions: widening the band would still pass, and narrowing
    it -- which is the point of this phase -- would pass without anyone
    having to lower the pin in the same commit.
    """
    widths: dict[str, dict[int, int]] = {quadrant: {} for quadrant in QUADRANTS}
    verdicts: collections.Counter[str] = collections.Counter()
    for quadrant, row, band in _band_profile():
        widths[quadrant][row] = len(band)
        verdicts.update(name for _, name in band)

    assert sum(sum(w.values()) for w in widths.values()) == 12259
    assert dict(verdicts) == {
        "NonExistentCellError": 12243,
        "DomainError": 13,
        "AntemeridianError": 2,
        "GeometryError": 1,
    }

    # Away from the four extension-zone rows the band is narrow, and it
    # is a share of the row rather than a constant: four columns is a
    # thousandth of the equatorial row and four fifths of row 998.
    assert widths["NE"][0] == 2
    assert widths["NE"][300] == 3
    assert widths["NE"][600] == 4
    assert widths["NE"][998] == 4

    # Eight rows are refused whole, all of them at the pole.
    whole = [
        (quadrant, row)
        for quadrant, rows in widths.items()
        for row, width in rows.items()
        if width == rows[row] and width == _row_length(quadrant, row)
    ]
    assert sorted(whole) == [
        ("NE", 999),
        ("NE", 1000),
        ("NW", 998),
        ("NW", 999),
        ("SE", 999),
        ("SE", 1000),
        ("SW", 998),
        ("SW", 999),
    ]


def _row_length(quadrant: str, row: int) -> int:
    return itacart.last_lattice_column(quadrant, row, SIDE) - _first(quadrant) + 1


@pytest.mark.slow
def test_the_band_width_is_the_rate_at_which_the_row_above_shortens() -> None:
    """Where the third and fourth refused columns come from.

    Two of them were derived already: a square of side ``s`` at column
    ``c`` spans the sheared coordinate over ``[(c - 1) s, (c + 1) s]``, so
    the outermost cell and its inner neighbour both reach the strip the
    screen protects. Three and four were not derived, and the reason they
    were not is that they do not come from this row at all.

    A cell's bounding box in lattice coordinates reaches the row above,
    and that row is shorter. Its own outermost cell is anomalous too, and
    it is the one that refuses. So the band is not a property of the row's
    last column; it is the difference between this row's last column and
    the next one's, plus the one column of the derivation above.

    Measured against every row of every quadrant, the derivation lands
    within a single column on all but four, and the residue is pinned by
    tally rather than by a bound so that it cannot drift in either
    direction: one column over where the neighbouring cell only touches
    the strip at a corner, one column under where the strip's outer bound
    does not quite reach.
    """
    residues: collections.Counter[int] = collections.Counter()
    rising: list[tuple[str, int]] = []
    for quadrant, row, band in _band_profile():
        last = itacart.last_lattice_column(quadrant, row, SIDE)
        above = itacart.last_lattice_column(quadrant, row + 1, SIDE)
        if above > last:
            # The row above is longer, which happens only where an
            # extension zone opens. The derivation assumes it is shorter.
            rising.append((quadrant, row))
            continue
        residues[len(band) - (last - above + 1)] += 1

    assert sorted(rising) == [("NE", 708), ("NW", 799), ("SE", 170), ("SW", 237)]
    assert dict(residues) == {0: 3382, 1: 368, -1: 248}
    assert sum(residues.values()) == 3998


@pytest.mark.slow
def test_the_widest_bands_sit_where_an_extension_zone_begins() -> None:
    """Four rows carry a band an order of magnitude wider than the rest.

    Elsewhere the band is one to five columns. On these four it is 23 to
    55, and they are exactly the rows where the last column jumps. Read
    without the zones this looks like a defect in the filling; it is the
    filling meeting a row whose neighbour has a different width.
    """
    widest = {
        ("NE", 799): 40,
        ("NW", 708): 55,
        ("SE", 237): 23,
        ("SW", 170): 23,
    }
    for (quadrant, row), expected in widest.items():
        assert len(_refused_band(quadrant, row)) == expected

    # And the row on either side is back to the ordinary width.
    for (quadrant, row), _ in widest.items():
        assert len(_refused_band(quadrant, row - 1)) <= 5
        assert len(_refused_band(quadrant, row + 1)) <= 5


@pytest.mark.slow
def test_the_refused_band_carries_no_exception_the_package_does_not_own() -> None:
    """Every refusal in the band is one of ours, by name.

    One was not. A single cell reached a caller carrying the geometry
    library's own topology exception, through a public name, and it did
    so because the intersection that splits a boundary by quadrant was
    handed a self-intersecting shape. The property is asserted over the
    whole band rather than for that cell, so the next one cannot appear
    quietly.
    """
    names = {name for _, _, band in _band_profile() for _, name in band}
    for name in names:
        assert issubclass(
            getattr(itacart.exceptions, name), ITACaRTError
        ), f"{name} is not part of this package's exception family"


@pytest.mark.slow
def test_the_band_gives_one_diagnosis_for_one_situation() -> None:
    """Fifty-seven cells used to be refused for the wrong reason.

    They are outermost cells like any other and the ordinary refusal is
    the true one, but a different exception arrived first and hid it.
    Cutting a densified boundary along an extension zone's own edge, when
    the boundary lies on that edge, left the union carrying the figure
    plus a train of zero-area pieces, and a collection is not something
    the projection accepts. Every one of the fifty-seven was the last
    column of its row in a western quadrant inside the zone latitudes,
    which is exactly where a cell's outer edge coincides with the cut.

    So the assertion is that the outermost family draws the outermost
    refusal, everywhere, rather than that a count went to zero.
    """
    for quadrant, row, band in _band_profile():
        last = itacart.last_lattice_column(quadrant, row, SIDE)
        for column, name in band:
            if column == last and name != "NonExistentCellError":
                assert (quadrant, row) in (
                    ("NE", 999),
                    ("NE", 1000),
                    ("NW", 998),
                    ("NW", 999),
                    ("SE", 999),
                    ("SE", 1000),
                    ("SW", 998),
                    ("SW", 999),
                ), f"{quadrant}({column:04d}/{row:04d}) refused as {name}"


def test_a_cell_whose_edge_lies_on_a_zone_edge_is_still_a_polygon() -> None:
    """The regression, at the seam where the debris was made.

    ``NW(0817/0714)`` has its outer edge on the Chukotka limit exactly.
    Before the fix the lift returned a collection of one polygon, one
    polygon of area under a ten-thousandth of a millionth, and forty-nine
    linestrings of no area at all. The property asserted is areal in,
    areal out -- a rule about type, so there is no threshold to choose
    and none to go stale.
    """
    from itacart.geometry import _auto_segment, _densify_any, _lift_extensions

    cell = "NW(0817/0714)"
    boundary = Polygon(itacart.cell_to_boundary(cell))
    assert min(x for x, _ in boundary.exterior.coords) == pytest.approx(-169.5)

    lifted = _lift_extensions(_densify_any(boundary, _auto_segment(2)))
    assert lifted.geom_type in {"Polygon", "MultiPolygon"}

    # And the cell now draws the refusal its position earns, which is the
    # one every other outermost cell draws.
    with pytest.raises(itacart.exceptions.NonExistentCellError):
        itacart.count_internal_cells(boundary, 2)


def test_the_polar_row_is_refused_everywhere_but_not_by_the_same_route() -> None:
    """Three quadrants refuse on the domain, one on the geometry.

    The polar row is refused by design: it is clipped by the pole and
    does not carry the nominal area. Three quadrants say exactly that.
    The fourth never reaches the check, because splitting its boundary by
    quadrant fails first.

    The cause is measured rather than inferred, and it is the same in all
    four. The ring of a polar-row cell stops short of the pole, and its
    densified form does not: an edge spanning half a turn of longitude is
    filled in along the geodesic joining its ends, and that geodesic runs
    through the pole. All four densified rings reach ninety degrees. Three
    of them come back without crossing themselves and one does, so the
    quadrant that behaves differently is not exhibiting a property of the
    grid -- it is the same defect landing on the wrong side of a validity
    test.
    """
    verdicts = {}
    for quadrant in QUADRANTS:
        cell = f"{quadrant}(0001/0999)"
        ring = itacart.cell_to_boundary(cell)
        polygon = Polygon(ring)
        assert polygon.is_valid
        assert max(abs(point[1]) for point in ring) < 89.99

        # The property is that densification leaves the cell, not that it
        # lands on a literal. It arrives a rounding short of ninety, and
        # asserting the printed value would pin the formatting instead.
        densified = itacart.densify_orthodromic(polygon, 1000.0)
        reached = max(abs(y) for _, y in densified.exterior.coords)
        assert reached > max(abs(point[1]) for point in ring)
        assert reached > 89.9999

        try:
            itacart.count_internal_cells(polygon, 2)
            verdicts[quadrant] = None
        except Exception as exc:
            verdicts[quadrant] = type(exc).__name__

    assert verdicts == {
        "NE": "DomainError",
        "NW": "DomainError",
        "SW": "DomainError",
        "SE": "GeometryError",
    }
