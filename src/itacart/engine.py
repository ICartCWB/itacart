"""The OGC-facing engine object.

Carries the requirements of OGC 20-040r3 that are properties of the
system rather than of any one cell -- 6 (harmonized model), 7 (defined
CRS) and 21's declaration clause -- and the self-assessment against the
whole requirement set, which have nowhere else to live.

The governing standard is OGC 20-040r3, *Topic 21 -- Discrete Global Grid
Systems -- Part 1: Core Reference System and Operations and Equal Area
Earth Reference System* (2021), the edition ISO 19170-1 follows. Its
predecessor OGC 15-104r5 (2017) stops at Requirement 18 and has no EAERS
class, so the requirement numbers used here exist only in 20-040r3.

The functional API is complete on its own; this class exists to expose
ITACaRT as a described, introspectable DGGS to standards-aware consumers
such as an OGC API-DGGS server.

Provenance: ``itacart_core/engine.py`` (``IDGGSEngine``), widened here to
cover the declarative requirements.
"""

from __future__ import annotations

from typing import Any

from . import constants
from .cells import cell_to_boundary, cell_to_centroid, geo_to_cell
from .exceptions import NonAtomicIndexError, ResolutionError
from .resolutions import resolution_table

__all__ = ["ITACaRT", "describe", "crs", "conformance"]

#: Default densification bound, in metres, shared by the class and the
#: geometry module's own default so that a caller who fixes nothing gets
#: the same edges from either route.
DEFAULT_MAX_SEGMENT_M = 1000.0

#: How much of the governed domain a requirement holds over. Descriptive
#: only: it never softens ``conformant``. What it adds is the size of an
#: exception, so that a requirement failing over three cells in five
#: million is not read as failing over the grid.
_COVERAGE = ("full", "partial", "none")


def _record(
    identifier: str, uri: str, conformant: bool, coverage: str, justification: str
) -> dict[str, object]:
    """One conformance record.

    ``conformant`` is a boolean and not a three-valued status, because
    every requirement of the standard is written with *shall*. A
    requirement either holds across the domain it governs or it does not,
    and an implementation that holds it almost everywhere has not met it.
    "Partially met" is a useful description of an implementation and an
    invalid answer to a *shall*, so it survives here as ``coverage`` and
    never as a verdict.
    """
    if coverage not in _COVERAGE:
        raise ValueError(f"unknown coverage: {coverage!r}")
    if conformant and coverage != "full":
        raise ValueError(
            f"requirement {identifier}: conformant with {coverage} coverage is "
            "not a claim this vocabulary can make"
        )
    return {
        "requirement": identifier,
        "uri": f"http://www.opengis.net/spec/DGGS/2.0/req/{uri}",
        "conformant": conformant,
        "coverage": coverage,
        "justification": justification,
    }


class ITACaRT:
    """Stateful facade over the ITACaRT reference system.

    Wraps the functional API so a caller can fix defaults once and reuse
    them, and provides the descriptive metadata the OGC model asks for.

    Args:
        default_resolution: Resolution assumed when a call omits one.
        edge_model: Default edge interpretation for geometry operations.
        max_segment_m: Default densification bound in metres.

    Raises:
        ResolutionError: If ``default_resolution`` is outside the table.
        ValueError: If ``edge_model`` is not an edge type ITACaRT defines,
            or ``max_segment_m`` is not positive.
    """

    #: Edge types this reference system declares, in the vocabulary of
    #: clause 9.1.4.2 of the standard. ITACaRT states two: its cell edges
    #: are geodesics on the WGS84 ellipsoid, and the rows and the polar cap
    #: close on small circles of constant latitude.
    EDGE_MODELS = ("WGS84_GEODESIC", "WGS84_SMALL_CIRCLE")

    def __init__(
        self,
        default_resolution: int = 13,
        edge_model: str = "WGS84_GEODESIC",
        max_segment_m: float = DEFAULT_MAX_SEGMENT_M,
    ) -> None:
        if (
            not constants.MIN_RESOLUTION
            <= default_resolution
            <= constants.MAX_RESOLUTION
        ):
            raise ResolutionError(
                f"resolution {default_resolution} is outside "
                f"{constants.MIN_RESOLUTION}..{constants.MAX_RESOLUTION}"
            )
        if edge_model not in self.EDGE_MODELS:
            raise ValueError(f"unknown edge model: {edge_model!r}")
        if not max_segment_m > 0.0:
            raise ValueError(f"max_segment_m must be positive, got {max_segment_m!r}")
        self.default_resolution = default_resolution
        self.edge_model = edge_model
        self.max_segment_m = max_segment_m

    def __repr__(self) -> str:
        return (
            f"ITACaRT(default_resolution={self.default_resolution}, "
            f"edge_model={self.edge_model!r}, max_segment_m={self.max_segment_m!r})"
        )

    # -- Descriptive (OGC Req 6, 7) -------------------------------------

    def describe(self) -> dict[str, Any]:
        """Machine-readable description of the reference system.

        Covers identity and DOI, the CRS and datum, the tessellation
        method, cell geometry, the full resolution table, refinement
        ratios, and the boundary treatments. Shaped to feed an OGC
        API-DGGS ``/dggs/{dggrsId}`` response.

        The content is read from the package rather than restated: the
        resolution table comes from :func:`itacart.resolution_table`, which
        transcribes Table 1 of the paper, and the ellipsoid from
        :mod:`itacart.constants`, which transcribes section 3. Anything
        here that duplicated a value instead of reading it would be a
        second place for that value to be wrong.

        ``constraints`` carries ``cellEqualSized``, which Requirement 21
        asks the reference system to declare, alongside the exception that
        makes the declaration honest.

        Returns:
            The description as a plain mapping.
        """
        from . import __paper_doi__ as paper_doi
        from . import __version__ as version

        return {
            "identifier": "ITACaRT",
            "title": "ITA Cadastral Ellipsoidal Reference Tessellation",
            "version": version,
            "paper_doi": paper_doi,
            "standard": {
                "document": "OGC 20-040r3",
                "title": (
                    "Topic 21 - Discrete Global Grid Systems - Part 1: Core "
                    "Reference System and Operations and Equal Area Earth "
                    "Reference System"
                ),
                "also_published_as": "ISO 19170-1",
            },
            "crs": self.crs(),
            "tessellation": {
                "method": "direct surface tessellation",
                "projection": "ellipsoidal sinusoidal (parallels plane)",
                "base_unit_polyhedron": None,
                "note": (
                    "Cells are tessellated directly on the WGS84 ellipsoid from "
                    "base and height measured on it, per section 3 of the paper. "
                    "No polyhedron is mapped to the surface, which is why "
                    "Requirement 22 is not met."
                ),
            },
            "domain": {
                "extent": "whole surface of the WGS84 ellipsoid",
                "dimensionality": {"spatial": 2, "temporal": 0, "topological": 2},
                "level_zero": list(constants.QUADRANTS),
            },
            "constraints": {
                "cellEqualSized": True,
                "exceptions": ["trapezoid", "polar cap", "polar row"],
            },
            "cell_geometry": {
                "shapes": ["parallelogram", "triangle", "trapezoid"],
                "edge_models": list(self.EDGE_MODELS),
                "default_edge_model": self.edge_model,
                "representative_position": "centroid",
            },
            "resolutions": {
                "count": constants.RESOLUTION_COUNT,
                "min": constants.MIN_RESOLUTION,
                "max": constants.MAX_RESOLUTION,
                "default": self.default_resolution,
                "table": list(resolution_table()),
                "refinement": {
                    "even_levels": 4,
                    "odd_levels": 25,
                    "note": (
                        "Resolution 0 is a global quadrant and resolution 1 is "
                        "the Cartesian base grid; neither is produced by "
                        "refining anything, so neither carries a ratio."
                    ),
                },
            },
            "index": {
                "method": "hierarchy-based, compositional",
                "grammar": "QQ(XXXX/YYYY(refinement(refinement...)))",
                "quadrants": list(constants.QUADRANTS),
                "example": "SE(1400/0374(3(C2(3))))",
                "note": (
                    "Section 3.1 of the paper. Column 0 does not exist in the "
                    "western quadrants, so the canonical spelling of a cell on "
                    "the prime meridian is the eastern one."
                ),
            },
            "boundary_treatments": {
                "prime_meridian": (
                    "The first column of each quadrant is a triangle, bounded "
                    "by the prime meridian."
                ),
                "antemeridian": (
                    "Not crossed. The last column of a row is a trapezoid "
                    "clipped on the antemeridian, and the two extension zones "
                    "carry the inhabited land that would otherwise straddle it."
                ),
                "extension_zones": {
                    name: {
                        "quadrant": zone.quadrant,
                        "lon_limit": zone.lon_limit,
                        "lat_min": zone.lat_min,
                        "lat_max": zone.lat_max,
                        "description": zone.description,
                    }
                    for name, zone in constants.EXTENSION_ZONES.items()
                },
                "poles": (
                    "Each pole is a single cap cell, in the eastern quadrant "
                    "only: NE(0000/1000) and SE(0000/1000). There is no cap in "
                    "NW or SW."
                ),
            },
            "geometry_defaults": {"max_segment_m": self.max_segment_m},
        }

    def crs(self) -> dict[str, Any]:
        """The coordinate reference system this grid is defined on.

        Reports WGS84 as the datum, satisfying requirement 7 and the
        GNSS-compatibility design criterion, alongside the ellipsoidal
        sinusoidal projection used internally.

        The ellipsoid parameters are read from :mod:`itacart.constants`,
        which transcribes section 3 of the paper and derives everything
        from the two defining parameters rather than transcribing decimals.

        Returns:
            A mapping with datum, ellipsoid parameters, the PROJ string
            and the projection's defining equations.
        """
        return {
            "datum": "WGS84",
            "authority": "EPSG",
            "code": 4326,
            "coordinate_epoch": None,
            "ellipsoid": {
                "name": "WGS 84",
                "semi_major_axis_m": constants.WGS84_A,
                "inverse_flattening": constants.WGS84_INV_F,
                "semi_minor_axis_m": constants.WGS84_B,
                "flattening": constants.WGS84_F,
                "first_eccentricity_squared": constants.WGS84_E2,
            },
            "projection": {
                "name": "ellipsoidal sinusoidal (parallels plane)",
                "proj_string": constants.SINUSOIDAL_PROJ,
                "equations": {
                    "x": "lambda * nu(phi) * cos(phi)",
                    "y": "M(phi), the meridian arc from the equator",
                },
                "equal_area": True,
                "note": (
                    "Equations (1) and (2) of section 3 of the paper, computed "
                    "directly. The PROJ string is given for interoperability; "
                    "PROJ is not called at runtime and is not a dependency."
                ),
            },
        }

    def conformance(self) -> dict[str, Any]:
        """Self-reported conformance against DGGS Core and EAERS.

        Reports one record per requirement of OGC 20-040r3, each with its
        URI, a boolean verdict, a coverage description and a justification.

        The verdict is boolean because the requirements are written with
        *shall*, which admits no middle answer: a requirement that holds
        for every cell but three is not met. Coverage is reported beside it
        and never instead of it, because the size of an exception is
        information a consumer needs and is not a softening of the verdict.
        Where a requirement is not met, the justification names what
        diverges and, wherever the package has measured it, by how much and
        over how many cells.

        Three things this deliberately does not do. It does not repeat the
        paper's grouping of Requirements 22 to 25 under one heading and one
        verdict: the figures of the standard put a polyhedron in the path
        of all four, but the normative text puts it in the path of only the
        first, and the other three are met. It does not write any
        limitation of an algorithm as a property of the grid: where a
        refusal belongs to an operation rather than to the tessellation,
        the record says which. And it does not answer a requirement about
        the figure on the surface of the Earth by measuring the chart the
        package happens to draw it in, which is how the polar caps were
        once read out of Requirements 11 and 26 that they in fact meet.

        Returns:
            A mapping with the two conformance classes and, under
            ``requirements``, one record per requirement, each with an
            identifier, a status and a justification.
        """
        core = [
            _record(
                "6",
                "core/rs/harmonized_model",
                True,
                "full",
                "describe() reports the reference system against the DGGS Core "
                "RS data model of Figure 13 and Tables 40 to 46.",
            ),
            _record(
                "7",
                "core/rs/crs",
                True,
                "full",
                "crs() reports WGS84 as the datum, which is also the "
                "GNSS-compatibility design criterion of the paper. No coordinate "
                "epoch is provided, the datum being static.",
            ),
            _record(
                "8",
                "core/rs/global_domain",
                True,
                "full",
                "The domain is the whole surface of the WGS84 ellipsoid, "
                "two-dimensional in space, with no temporal dimension.",
            ),
            _record(
                "9",
                "core/rs/global_domain/complete",
                True,
                "full",
                "Level zero is the four quadrants, which cover the ellipsoid. "
                "Below it, an enumerated lattice of 2376 positions has every "
                "position inside the plane figure of the cell that addresses it.",
            ),
            _record(
                "10",
                "core/rs/global_domain/unique",
                True,
                "full",
                "geo_to_cell answers exactly once, and four enumerated rows of "
                "resolution 1 -- 2004, 1785, 312 and 2 cells -- have no adjacent "
                "pair overlapping in area.",
            ),
            _record(
                "11",
                "core/rs/cell/simple",
                True,
                "full",
                "Clause 8.2.4.1 asks that a cell not self-intersect, be "
                "topologically a circle, and enclose a region measurable in "
                "its own dimension. Every cell family satisfies all three as a "
                "region of the ellipsoid surface, the polar caps included: a "
                "cap is bounded by one parallel, which is a simple closed "
                "curve, is homeomorphic to a disc, and has a surface area "
                "agreeing with the closed form to two parts in a hundred "
                "million.",
            ),
            _record(
                "12",
                "core/rs/cell/direct_position",
                True,
                "full",
                "The direct position is the centroid, which lies within the "
                "boundary in every family measured. The anchor, which the paper "
                "offers for this requirement, is a vertex and so lies on the "
                "boundary rather than within it.",
            ),
            _record(
                "13",
                "core/rs/cell/address",
                True,
                "full",
                "The compositional index is a globally unique identifier per "
                "zone, enumerated over whole rows without collision. A second "
                "spelling of a prime-meridian cell is accepted and rewritten "
                "into the canonical eastern one; that is one zone with two "
                "spellings, not two zones with one name.",
            ),
            _record(
                "14",
                "core/rs/discrete_global_grid",
                True,
                "full",
                "Each resolution is the aggregation of all cells at that level.",
            ),
            _record(
                "15",
                "core/rs/discrete_global_grid/sequence",
                True,
                "full",
                "Fourteen levels in increasing refinement order, 10 km to 1 cm.",
            ),
            _record(
                "16",
                "core/functions/quantization",
                True,
                "full",
                "polyfill and vertex_to_cell quantize areal, linear and point "
                "geometry; decompose and cells_to_geometry read it back.",
            ),
            _record(
                "17",
                "core/functions/query/zonequery",
                True,
                "full",
                "Parent, child, neighbour, ring, disk and distance resolve from "
                "the index alone. grid_distance refuses three named cases: a "
                "trapezoid endpoint, an antemeridian crossing, and a "
                "non-atomic index.",
            ),
            _record(
                "18",
                "core/functions/interoperation/query",
                True,
                "full",
                "from_geojson, recover_from_geojson and from_geodataframe read "
                "an external query back into cell sets.",
            ),
            _record(
                "19",
                "core/functions/interoperation/broadcast",
                True,
                "full",
                "cells_to_geojson, cell_to_wkt, cells_to_wkt and to_geodataframe "
                "emit standard formats, with the index written into the GeoJSON "
                "feature id that section 3.2 of RFC 7946 reserves for it.",
            ),
        ]

        eaers = [
            _record(
                "20",
                "ea/ers/harmonized_model",
                False,
                "partial",
                "The data models of Figures 20 and 22 are satisfied. Figure 21 "
                "is not: it describes a polyhedral interface, and ITACaRT "
                "tessellates the ellipsoid directly, so there is nothing for "
                "that interface to describe. Two models out of three, which a "
                "shall does not accept.",
            ),
            _record(
                "21",
                "ea/ers/global_domain",
                True,
                "full",
                "The domain is the whole surface of the Earth model, and "
                "describe() declares cellEqualSized among its constraints "
                "together with the three cell families that are its exceptions.",
            ),
            _record(
                "22",
                "ea/ers/tessellation/initial",
                False,
                "none",
                "No base unit polyhedron is mapped to the surface. The initial "
                "tessellation is direct on the ellipsoid, which is the design "
                "choice the paper makes for absolute geodetic fidelity, and it "
                "is a property of the reference system rather than a gap.",
            ),
            _record(
                "23",
                "ea/ers/tessellation/sequence",
                True,
                "full",
                "Refinement generates grids with strictly smaller cells at every "
                "level, 10 km down to 1 cm.",
            ),
            _record(
                "24",
                "ea/ers/tessellation/sequence/max",
                True,
                "full",
                "The sequence stops at resolution 13; a finer level is refused.",
            ),
            _record(
                "25",
                "ea/ers/tessellation/global_area_preservation",
                True,
                "full",
                "The requirement offers two branches and says one shall apply. "
                "This CRS has a static datum, so the branch that governs asks "
                "that domain completeness and position uniqueness hold "
                "throughout the sequence for all cells within their respective "
                "grids. Measured one grid down, across parents and at the "
                "antemeridian: no overlapping pair, and every parent covered by "
                "its own children. The children of a trapezoid cover 1.000248 "
                "of their parent's plane figure, and every point of that "
                "surplus is addressed to the same trapezoid at resolution 1 and "
                "to one of its children at resolution 2, so the grid owns it at "
                "both levels and owns it consistently.",
            ),
            _record(
                "26",
                "ea/ers/cell/simple/2d_polygon",
                True,
                "full",
                "Clause 9.1.5.1 asks that edges meet only at vertices, that "
                "exactly two edges meet at each vertex, that edges and vertices "
                "be equal in number, and that the region have measurable area. "
                "The criteria are about the figure on the surface model of the "
                "Earth and say nothing about straight segments, geodesic edges "
                "or a change of tangent at a vertex; a parallel is a small "
                "circle, which is an allowed edge type. Parallelograms, "
                "triangles and trapezoids satisfy the four criteria directly. "
                "The polar cap satisfies them with its bounding parallel "
                "partitioned into arcs: for any k, k arcs meet only at k "
                "vertices, two at each, and the region has area.",
            ),
            _record(
                "27",
                "ea/ers/cell/direct_position/centroid",
                False,
                "partial",
                "cell_to_centroid is the direct position, it lies on the "
                "surface, and it is within the cell everywhere. Over the "
                "ordinary families it agrees with the geodesic centre of "
                "surface area to within ten metres, which on a cell 10 km "
                "across is a thousandth of the cell; over the parallelogram "
                "families the agreement is 0.7 m and stable under refinement "
                "of the search. It is nonetheless the area centroid of the "
                "equal-area projection plane inverted to geodetic, and the "
                "requirement names the geodesic centre. The two part company "
                "at the polar row and the two caps, where the cap's geodesic "
                "centre is the pole itself and the plane answer sits about "
                "1.3 km short of it: three orders of magnitude, and three "
                "families out of ten. A shall is not met by seven.",
            ),
            _record(
                "28",
                "ea/ers/cell/equal_area/error_budget",
                False,
                "partial",
                "The requirement fixes a budget of 1 per cent or less against "
                "the theoretical average cell area, which at resolution 1 is "
                "100 008 533 m squared over 5 100 221 cells. Formal conformance "
                "is no, because a budget covering the whole domain cannot be "
                "stated. Substantive coverage is partial, and the exception is "
                "small and named: 3941 antemeridian trapezoids and the 2 polar "
                "caps, 3943 cells or 0.077 per cent of the grid. Every other "
                "cell, parallelogram and meridian triangle alike, sits at "
                "0.999915 of the theoretical average, which is 0.0085 per cent "
                "off and well inside the budget. The trapezoid family runs from "
                "0.032 to 2.067 and the caps sit at 0.121; those figures "
                "describe that family, not the grid.",
            ),
            _record(
                "29",
                "ea/ers/cell/equal_area",
                False,
                "partial",
                "The projection is equal-area by construction, so every "
                "ordinary cell of a resolution has the same area as every "
                "other, measured uniformly at 0.999915 of the theoretical "
                "average. The requirement is per cell geometry, and the "
                "trapezoid geometry is not internally equal-area: clipping on "
                "the antemeridian leaves a family running from 0.032 to 2.067. "
                "The polar caps carry the same failure into the triangle "
                "geometry, at 0.121 against the meridian triangles' 0.999915.",
            ),
        ]

        return {
            "standard": "OGC 20-040r3 (ISO 19170-1)",
            "classes": {
                "core": "http://www.opengis.net/spec/DGGS/2.0/conf/core",
                "eaers": "http://www.opengis.net/spec/DGGS/2.0/conf/ea/ers",
            },
            "summary": {
                "conformant": sum(bool(r["conformant"]) for r in core + eaers),
                "non_conformant": sum(not r["conformant"] for r in core + eaers),
                "core_conformant": all(r["conformant"] for r in core),
                "eaers_conformant": all(r["conformant"] for r in eaers),
                "coverage": {
                    value: sum(r["coverage"] == value for r in core + eaers)
                    for value in _COVERAGE
                },
            },
            "requirements": core + eaers,
            "note": (
                "Requirements 1 to 5 belong to the Common Classes module and "
                "describe conceptual schemas rather than a reference system's "
                "behaviour; they are outside what this package can assert about "
                "itself and are not reported."
            ),
        }

    # -- Delegating operations ------------------------------------------

    def geo_to_cell(self, lon: float, lat: float, resolution: int | None = None) -> str:
        """Address the cell containing a position.

        See :func:`itacart.cells.geo_to_cell`.
        """
        return geo_to_cell(
            lon, lat, self.default_resolution if resolution is None else resolution
        )

    def cell_to_centroid(self, cell: str) -> tuple[float, float]:
        """Geodetic centroid of a cell.

        See :func:`itacart.cells.cell_to_centroid`.

        Raises:
            NonAtomicIndexError: If the index names more than one cell.
                The functional API answers a compositional index with a
                positionally aligned list; this signature promises one
                answer, and returning the first of several would be wrong
                in a way the type checker endorses. The refusal uses the
                package's own exception rather than a bare one, because
                ``are_neighbor_cells`` and ``grid_distance`` refuse the
                identical argument shape with this exception and a caller
                should not have to know which of the three it called.
        """
        centroid = cell_to_centroid(cell)
        if not isinstance(centroid, tuple):
            raise NonAtomicIndexError(
                f"expected one cell, got a compositional index: {cell!r}"
            )
        return centroid

    def cell_to_boundary(self, cell: str) -> list[tuple[float, float]]:
        """Geodetic vertices bounding a cell.

        See :func:`itacart.cells.cell_to_boundary`.

        Raises:
            NonAtomicIndexError: If the index names more than one cell, for
                the reason given on :meth:`cell_to_centroid`.
        """
        ring = cell_to_boundary(cell)
        if isinstance(ring[0], list):
            raise NonAtomicIndexError(
                f"expected one cell, got a compositional index: {cell!r}"
            )
        return ring


def describe() -> dict[str, Any]:
    """Module-level shortcut to :meth:`ITACaRT.describe` with defaults."""
    return ITACaRT().describe()


def crs() -> dict[str, Any]:
    """Module-level shortcut to :meth:`ITACaRT.crs` with defaults."""
    return ITACaRT().crs()


def conformance() -> dict[str, Any]:
    """Module-level shortcut to :meth:`ITACaRT.conformance` with defaults."""
    return ITACaRT().conformance()
