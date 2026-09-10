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

#: The whole refused band, in columns, over every row of every quadrant.
#:
#: It used to be twelve thousand and its accounting needed three derived
#: terms plus a per-machine residue, because the fill refused the entire
#: last-column family and the columns whose cells reached the strip it
#: protected. E3-L taught it to descend that family, which left the three
#: rows at the pole; E3-P found that two of those three were refused for
#: properties they do not have, and the band is empty.
#:
#: Empty is a stronger pin than any number, because every way of getting
#: it wrong -- a row miscounted as polar, a candidate position outside the
#: grid cancelling a valid cell, a screen widened by accident -- adds to
#: it rather than subtracting.
_REFUSED_COLUMNS = 0

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
    # Six was the top of this distribution while the child walk reached
    # one column east and stopped. It reaches as far as the absorbing
    # column now, and forty-five members gained a seventh child spelled
    # two columns out, which is where they came from: the six-count fell
    # from 839 to 794 by exactly that many.
    assert collections.Counter(n for n in off.values()) == {
        1: 2,
        2: 423,
        3: 752,
        5: 836,
        6: 794,
        7: 45,
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

    # Every refusal that remains is structural -- a polar row, or the
    # cap's own footprint against the antemeridian -- so none of them is
    # decided by a predicate on a marginal overlap and the whole profile
    # is pinned outright. ``NonExistentCellError`` has left the profile
    # entirely: it was the last lattice column, and the fill descends it.
    assert not verdicts, dict(verdicts)

    total = sum(sum(w.values()) for w in widths.values())
    assert total == _REFUSED_COLUMNS

    banded = sorted(
        (quadrant, row, width)
        for quadrant, rows in widths.items()
        for row, width in rows.items()
        if width
    )
    assert banded == []


def _row_length(quadrant: str, row: int) -> int:
    return itacart.last_lattice_column(quadrant, row, SIDE) - _first(quadrant) + 1


@pytest.mark.slow
def test_the_rows_where_a_zone_opens_are_the_rows_where_the_lattice_lengthens() -> None:
    """Four rows have a longer row above them, and they are the zones.

    This is what survives of a derivation that used to explain the width
    of the refused band. The band was the difference between a row's last
    column and the next one's, because a cell's bounding box in lattice
    coordinates reaches the row above and that row's own outermost cell
    was refused too. None of that is a fact about the grid; it was a fact
    about a screen that no longer refuses.

    What is a fact about the grid is where the lattice stops shortening
    and lengthens instead, which happens only where an extension zone
    carries the domain past the line. That is integer arithmetic over
    ``last_lattice_column`` and it moves on no machine, so it is pinned
    here on its own, without a fill in the assertion at all.

    The band is asserted empty on those rows: they were the widest in the
    profile, at 23 to 55 columns, and they are the rows a reader would
    check first if the descent had quietly stopped working.
    """
    rising: list[tuple[str, int]] = []
    for quadrant, row, band in _band_profile():
        last = itacart.last_lattice_column(quadrant, row, SIDE)
        above = itacart.last_lattice_column(quadrant, row + 1, SIDE)
        if above > last:
            rising.append((quadrant, row))
            assert not band, (quadrant, row, band)

    assert sorted(rising) == [("NE", 708), ("NW", 799), ("SE", 170), ("SW", 237)]


@pytest.mark.slow
def test_the_rows_that_carried_the_widest_bands_carry_none() -> None:
    """The four widest bands in the old profile, now empty.

    They ran 23 to 55 columns where the rest of the grid ran one to five,
    and they sit where the last column jumps because a zone begins. Read
    without the zones it looked like a defect in the filling; it was the
    filling meeting a row whose neighbour has a different width, and
    refusing rather than descending. Pinned by name so that a descent
    which regressed on the zone seams would fail here and not only in the
    aggregate.
    """
    for quadrant, row in (("NE", 799), ("NW", 708), ("SE", 237), ("SW", 170)):
        assert _refused_band(quadrant, row) == ()
        assert _refused_band(quadrant, row - 1) == ()
        assert _refused_band(quadrant, row + 1) == ()


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
    """Every refusal that remains names the polar family, and only it.

    Fifty-seven cells used to be refused for the wrong reason: cutting a
    densified boundary along an extension zone's own edge, when the
    boundary lies on that edge, left the union carrying the figure plus a
    train of zero-area pieces, and a collection is not something the
    projection accepts. That is repaired, and the family those cells
    belonged to -- the last lattice column -- is no longer refused at
    all, so the diagnosis it used to draw has left the profile.

    The assertion is now the complement: nothing outside the polar rows
    is refused for any reason, and inside them the reason is the domain
    or the antemeridian rather than a structural non-existence. A
    descent that regressed into refusing an ordinary column would fail
    here with the cell that did it named.
    """
    for quadrant, row, band in _band_profile():
        for column, name in band:
            assert row >= 998, (
                f"{quadrant}({column:04d}/{row:04d}) refused as {name} "
                f"outside the polar rows"
            )
            assert name in {
                "DomainError",
                "AntemeridianError",
            }, f"{quadrant}({column:04d}/{row:04d}) refused as {name}"


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

    # And the cell is filled, like every other outermost cell. It used to
    # be refused, and the refusal was the last thing standing between
    # this seam and the descent.
    assert itacart.count_internal_cells(boundary, 2) > 0


def test_row_999_densifies_over_the_pole_and_is_filled_anyway() -> None:
    """Four quadrants fill it, and the densification still crosses the pole.

    Row 999 used to be refused in all four quadrants as "the polar row".
    It is not the polar row: the pole falls in row 1000, and a cell of
    row 999 stops 1 966 metres short of it. The classification came from
    reading ``last_lattice_column`` for the row above and treating its
    answer of zero as "no cell", which is one cell in an eastern quadrant
    and none in a western one. These cells are ordinary border-absorbing
    trapezoids with six real children each, and they fill.

    What does not change is the geodesy, and it is measured rather than
    inferred, the same in all
    four. The ring of a polar-row cell stops short of the pole, and its
    densified form does not: an edge spanning half a turn of longitude is
    filled in along the geodesic joining its ends, and that geodesic runs
    through the pole. All four densified rings reach ninety degrees --
    that part is geodesy and stands. What did not stand was the longitude
    the walk gave the far side of the pole, which is what
    :func:`test_the_walk_over_the_pole_lands_on_the_branch_it_was_sent_to`
    pins.
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
        assert densified.is_valid, quadrant

        try:
            verdicts[quadrant] = itacart.count_internal_cells(polygon, 2)
        except Exception as exc:  # pragma: no cover - a regression would land here
            verdicts[quadrant] = type(exc).__name__

    # Filled, not refused, and by a count rather than by a name. The
    # count is small because the densified ring is not the cell -- an
    # edge spanning half a turn of longitude is filled in along a
    # geodesic that runs through the pole, so the figure the fill is
    # asked about at this latitude is not the figure the index names.
    # That loss is the subject of the interop measurements; what this
    # test pins is that the answer is a number in all four quadrants.
    assert all(isinstance(count, int) for count in verdicts.values()), verdicts
    assert all(count > 0 for count in verdicts.values()), verdicts


def test_the_walk_over_the_pole_lands_on_the_branch_it_was_sent_to() -> None:
    """A geodesic through the pole ends its walk where its segment ends.

    Two points half a turn of longitude apart lie on one meridian circle,
    so the geodesic joining them runs through the pole and the longitude
    changes by a hundred and eighty degrees at the crossing. Both
    ``+180`` and ``-180`` name that change, and the direct solution
    returns whichever its normalisation picks. Placing every interior
    point on the branch of the *start* left the far half free to land on
    the opposite branch from the vertex it was walking towards, which
    folds the ring across the globe.

    Three of the four polar-row cells survived that because their far
    vertex sits on the prime meridian, where both branches agree. The
    fourth walks towards a hundred and eighty and did not. The property
    asserted here is over all four, since one passing example is what the
    grid had before.
    """
    for quadrant in QUADRANTS:
        ring = itacart.cell_to_boundary(f"{quadrant}(0001/0999)", close=True)
        edges = [
            (ring[index], ring[index + 1])
            for index in range(len(ring) - 1)
            if abs(abs(ring[index + 1][0] - ring[index][0]) - 180.0) < 1e-9
        ]
        assert len(edges) == 1, f"{quadrant}: expected one edge spanning half a turn"

        start, end = edges[0]
        walk = itacart.densify_segment(start, end, 1000.0)
        assert len(walk) > 2, quadrant

        # The pole is reached, and it is reached at the midpoint because
        # the edge is symmetric about it. Everything after that midpoint
        # belongs to the destination's half turn, and the test asks for
        # the whole half rather than for the last point alone.
        # Distance to the pole, not equality with it. The walk arrives by
        # integration and stops a fraction of a micrometre short, so
        # asserting the printed ninety would pin the formatting rather
        # than the geodesy. One metre is three orders below the step this
        # walk takes, which makes the claim "it reaches the pole" and not
        # "it rounds well".
        shortfall = 90.0 - max(abs(latitude) for _, latitude in walk)
        assert shortfall * 111_319.0 < 1.0, (quadrant, shortfall)

        far = walk[len(walk) // 2 + 1 :]
        for longitude, _ in far:
            assert abs(longitude - end[0]) < 90.0, (quadrant, longitude, end[0])
