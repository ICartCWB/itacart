"""The engine facade, and the declarations it makes.

Three kinds of test live here, and they are different in kind.

The construction and delegation tests are ordinary: they check that the
facade fixes defaults, refuses arguments it cannot honour, and hands work
to the functional API rather than doing any of its own.

The description tests are transcription checks. ``describe`` and ``crs``
report what the paper says the system is, so each assertion here names the
section it comes from and compares against the package's own constants
rather than against a literal written twice.

The conformance tests are the load-bearing ones. ``conformance`` is what a
consumer reads to decide whether ITACaRT serves, so what is asserted is
not that it returns something well shaped -- though that too -- but that
every requirement of the standard has a record, that every record that is
not "met" says what diverges, and that the statuses agree with what the
conformance suite measures. A declaration that drifted from the
measurements would be worse than no declaration, because it would be
believed.

Provenance, carried forward from the placeholder this file replaces:
nothing here is portable from ``itacart_core``, which has no coverage for
this module. The 406 figure that once stood in that note was the size of
its whole suite rather than a count of anything reusable, measured in an
earlier phase with a grep that returned nothing.
"""

from __future__ import annotations

import inspect
import re

import pytest

import itacart
from itacart import constants, engine

#: Every requirement of OGC 20-040r3 that describes a reference system's
#: behaviour. 1 to 5 are the Common Classes module and describe conceptual
#: schemas, which this package cannot assert about itself.
REPORTED_REQUIREMENTS = tuple(str(number) for number in range(6, 30))


class TestConstruction:
    """Defaults are fixed once, and bad arguments are refused at the door."""

    def test_the_defaults_are_the_documented_ones(self) -> None:
        grid = engine.ITACaRT()
        assert grid.default_resolution == constants.MAX_RESOLUTION
        assert grid.edge_model == "WGS84_GEODESIC"
        assert grid.max_segment_m == engine.DEFAULT_MAX_SEGMENT_M

    def test_the_defaults_can_be_overridden(self) -> None:
        grid = engine.ITACaRT(5, "WGS84_SMALL_CIRCLE", 250.0)
        assert (grid.default_resolution, grid.edge_model) == (5, "WGS84_SMALL_CIRCLE")
        assert grid.max_segment_m == 250.0
        assert "WGS84_SMALL_CIRCLE" in repr(grid)

    @pytest.mark.parametrize("resolution", [-1, constants.MAX_RESOLUTION + 1])
    def test_a_resolution_off_the_table_is_refused(self, resolution: int) -> None:
        with pytest.raises(itacart.ResolutionError):
            engine.ITACaRT(resolution)

    def test_an_edge_model_the_system_does_not_define_is_refused(self) -> None:
        with pytest.raises(ValueError, match="edge model"):
            engine.ITACaRT(edge_model="GREAT_CIRCLE")

    @pytest.mark.parametrize("bound", [0.0, -1.0])
    def test_a_non_positive_densification_bound_is_refused(self, bound: float) -> None:
        with pytest.raises(ValueError, match="positive"):
            engine.ITACaRT(max_segment_m=bound)


class TestDelegation:
    """The facade computes nothing; it holds defaults and forwards."""

    def test_geo_to_cell_uses_the_default_resolution_when_none_is_given(self) -> None:
        grid = engine.ITACaRT(default_resolution=1)
        assert grid.geo_to_cell(50.0, 27.0) == itacart.geo_to_cell(50.0, 27.0, 1)

    def test_geo_to_cell_honours_an_explicit_resolution(self) -> None:
        grid = engine.ITACaRT(default_resolution=1)
        assert grid.geo_to_cell(50.0, 27.0, 3) == itacart.geo_to_cell(50.0, 27.0, 3)

    def test_the_two_geometry_delegations_match_the_functional_api(self) -> None:
        grid, cell = engine.ITACaRT(), "NE(0500/0300)"
        assert grid.cell_to_centroid(cell) == itacart.cell_to_centroid(cell)
        assert grid.cell_to_boundary(cell) == itacart.cell_to_boundary(cell)

    @pytest.mark.parametrize("method", ["cell_to_centroid", "cell_to_boundary"])
    def test_a_compositional_index_is_refused_by_the_single_cell_methods(
        self, method: str
    ) -> None:
        """The signatures promise one cell, so many cells are refused.

        The functional API answers a compositional index with a positionally
        aligned list. These two methods are typed as returning one answer,
        and silently returning the first of several would be worse than
        refusing: it would be wrong in a way the type checker endorses.
        """
        index = itacart.compose(["NE(0500/0300)", "NE(0501/0300)"])
        with pytest.raises(ValueError, match="compositional index"):
            getattr(engine.ITACaRT(), method)(index)


class TestDescription:
    """``describe`` and ``crs`` transcribe the paper, and are checked against it."""

    def test_the_crs_is_wgs84_read_from_the_constants(self) -> None:
        """Section 3 of the paper, and the GNSS-compatibility criterion."""
        reported = itacart.crs()
        assert reported["datum"] == "WGS84"
        assert reported["code"] == 4326
        ellipsoid = reported["ellipsoid"]
        assert ellipsoid["semi_major_axis_m"] == constants.WGS84_A
        assert ellipsoid["inverse_flattening"] == constants.WGS84_INV_F
        assert ellipsoid["semi_minor_axis_m"] == constants.WGS84_B
        assert ellipsoid["first_eccentricity_squared"] == constants.WGS84_E2

    def test_the_projection_is_named_and_proj_is_not_a_dependency(self) -> None:
        projection = itacart.crs()["projection"]
        assert projection["equal_area"] is True
        assert projection["proj_string"] == constants.SINUSOIDAL_PROJ
        assert "not called at runtime" in projection["note"]

    def test_the_resolution_table_is_read_and_not_restated(self) -> None:
        """Table 1 of the paper, through the package's own transcription."""
        resolutions = itacart.describe()["resolutions"]
        assert resolutions["table"] == list(itacart.resolution_table())
        assert resolutions["count"] == constants.RESOLUTION_COUNT == 14
        assert (resolutions["min"], resolutions["max"]) == (0, 13)
        assert resolutions["refinement"]["even_levels"] == 4
        assert resolutions["refinement"]["odd_levels"] == 25

    def test_the_description_declares_no_polyhedron(self) -> None:
        """Section 3: base and height are measured on the ellipsoid itself."""
        tessellation = itacart.describe()["tessellation"]
        assert tessellation["base_unit_polyhedron"] is None
        assert tessellation["method"] == "direct surface tessellation"

    def test_the_description_declares_the_equal_size_constraint_and_its_exceptions(
        self,
    ) -> None:
        """Requirement 21 asks for the declaration; honesty asks for the rest."""
        constraints = itacart.describe()["constraints"]
        assert constraints["cellEqualSized"] is True
        assert constraints["exceptions"]

    def test_the_boundary_treatments_name_the_zones_the_package_defines(self) -> None:
        treatments = itacart.describe()["boundary_treatments"]
        assert set(treatments["extension_zones"]) == set(constants.EXTENSION_ZONES)
        for name, zone in constants.EXTENSION_ZONES.items():
            assert treatments["extension_zones"][name]["lon_limit"] == zone.lon_limit

    def test_the_identity_is_read_from_the_package(self) -> None:
        described = itacart.describe()
        assert described["version"] == itacart.__version__
        assert described["paper_doi"] == itacart.__paper_doi__
        assert described["standard"]["document"] == "OGC 20-040r3"


class TestConformance:
    """The declaration, and whether it says what the measurements say."""

    def test_every_reported_requirement_has_exactly_one_record(self) -> None:
        records = itacart.conformance()["requirements"]
        identifiers = [record["requirement"] for record in records]
        assert identifiers == list(REPORTED_REQUIREMENTS)

    def test_every_record_carries_a_uri_a_verdict_and_a_justification(self) -> None:
        for record in itacart.conformance()["requirements"]:
            assert record["uri"].startswith("http://www.opengis.net/spec/DGGS/2.0/req/")
            assert isinstance(record["conformant"], bool)
            assert record["coverage"] in engine._COVERAGE
            assert len(record["justification"]) > 40, record["requirement"]

    def test_the_verdict_is_binary_because_the_requirements_say_shall(self) -> None:
        """No record may hedge, and none may claim to conform in part.

        Every requirement of the standard is written with *shall*, so the
        answer to it is yes or no: a requirement that holds for every cell
        but three is not met. Coverage sits beside the verdict to size an
        exception, and the record helper refuses the one combination that
        would let it soften the verdict instead.
        """
        for record in itacart.conformance()["requirements"]:
            if record["conformant"]:
                assert record["coverage"] == "full", record["requirement"]

        with pytest.raises(ValueError, match="not a claim"):
            engine._record("11", "core/rs/cell/simple", True, "partial", "...")

    def test_the_summary_counts_the_records_it_summarises(self) -> None:
        reported = itacart.conformance()
        summary = reported["summary"]
        records = reported["requirements"]
        assert summary["conformant"] + summary["non_conformant"] == len(records)
        assert summary["conformant"] == sum(r["conformant"] for r in records)
        assert summary["core_conformant"] is True
        assert summary["eaers_conformant"] is False
        assert sum(summary["coverage"].values()) == len(records)

    def test_an_unknown_coverage_cannot_be_recorded(self) -> None:
        """The vocabulary is closed, so a typo cannot invent a coverage."""
        with pytest.raises(ValueError, match="unknown coverage"):
            engine._record("6", "core/rs/harmonized_model", False, "mostly", "...")

    def test_the_divergences_are_the_ones_the_suite_measured(self) -> None:
        """The declaration and the conformance suite have to agree.

        Five requirements are not met and nineteen are, and every verdict
        is measured in ``tests/conformance``. Pinning them by number is what
        stops the declaration from drifting into optimism one edit at a
        time -- and, as Requirements 11, 25 and 26 show, from drifting into
        pessimism either. Three of the four the paper refuses as a block
        are met, and the paper's own reason for refusing them was the
        figures rather than the text.
        """
        by_id = {
            record["requirement"]: record
            for record in itacart.conformance()["requirements"]
        }
        assert [r for r in by_id if not by_id[r]["conformant"]] == [
            "20",
            "22",
            "27",
            "28",
            "29",
        ]
        for requirement in ("11", "12", "21", "23", "24", "25", "26"):
            assert by_id[requirement]["conformant"] is True, requirement

    @pytest.mark.parametrize("requirement", ["27", "28", "29"])
    def test_a_failure_over_a_named_family_says_how_large_it_is(
        self, requirement: str
    ) -> None:
        """An exception must be sized, or it will be read as the whole.

        Not met over three cells in five million and not met over the
        entire grid are the same two words, and a consumer reading only the
        words would take the second. Each of these three fails over a named
        family, so each justification carries the figure that bounds it.
        """
        by_id = {
            record["requirement"]: record
            for record in itacart.conformance()["requirements"]
        }
        record = by_id[requirement]
        assert record["conformant"] is False
        assert record["coverage"] == "partial"
        assert any(
            figure in record["justification"]
            for figure in ("0.999915", "0.7 m", "0.032")
        ), record["requirement"]

    def test_the_geometry_records_name_cells_and_never_functions(self) -> None:
        """The distinction the project has had to re-learn several times.

        Requirements 11, 26, 28 and 29 are about the tessellation itself:
        what shape a cell is and what area it holds. Nothing an operation
        refuses belongs in their justifications. ``polyfill`` declines
        trapezoids, the polar row and the caps, and ``get_children`` raises
        descending the cap from resolution 6, but those are limitations of
        algorithms, and writing one of them here would assert something
        false about the grid to a reader who never opens the code.

        The rule is enforced rather than remembered: no operation that
        takes a cell may be named in those four justifications.
        Requirements that really are about an operation -- 16, 17, 18, 19,
        27 -- are free to name one, and do.

        The forbidden set is derived rather than listed, from the first
        parameter of every public callable, plus the two quantizers that
        take a geometry and refuse cell families anyway. Deriving it means
        an operation added later is covered without anyone remembering to
        add it here. It also means the set holds operations and not English
        words: ``describe`` and ``conformance`` take no cell and are not in
        it, which is why a justification may say that a figure describes a
        family.
        """
        functions = set()
        for name in itacart.__all__:
            value = getattr(itacart, name)
            if not callable(value) or isinstance(value, type):
                continue
            try:
                first = next(iter(inspect.signature(value).parameters))
            except (StopIteration, ValueError):  # pragma: no cover - no parameters
                continue
            if first in ("cell", "cells", "index", "region"):
                functions.add(name)
        functions |= {"polyfill", "vertex_to_cell"}
        assert {"get_children", "polyfill", "uncompact_cells"} <= functions
        assert not {"describe", "conformance"} & functions

        by_id = {
            record["requirement"]: record
            for record in itacart.conformance()["requirements"]
        }
        for requirement in ("11", "26", "28", "29"):
            words = set(re.findall(r"[A-Za-z_]+", by_id[requirement]["justification"]))
            named = sorted(functions & words)
            assert named == [], f"requirement {requirement} names {named}"

    def test_the_class_and_the_shortcuts_report_the_same_thing(self) -> None:
        grid = engine.ITACaRT()
        assert grid.describe() == itacart.describe()
        assert grid.crs() == itacart.crs()
        assert grid.conformance() == itacart.conformance()
