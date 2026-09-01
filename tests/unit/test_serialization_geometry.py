"""Tests for :mod:`itacart.serialization.geometry_blob`."""

from __future__ import annotations

import pytest

import itacart
from itacart.exceptions import (
    GeometryError,
    IncompatibleProfileError,
    MalformedBlobError,
    UnsupportedGeometryTypeError,
)
from itacart.serialization import geometry_blob as gb
from itacart.serialization import tree_blob as tb

RING = [
    "NE(0625/0451(2(A1)))",
    "NE(0625/0451(1(A2)))",
    "NE(0625/0451(1(A1)))",
]
LINE = ["NE(0625/0451(1(A1)))", "NE(0625/0451(1(A2)))"]
POINT = ["NE(0625/0451(1(A1)))"]


def _polygon(rings: list[list[str]] | None = None, **kwargs: object) -> bytes:
    return gb.encode_geometry(
        rings if rings is not None else [RING], "POLYGON", resolution=3, **kwargs
    )


# --------------------------------------------------------------------------
# Header, nine bytes, as the specification lays it out
# --------------------------------------------------------------------------


def test_the_header_is_nine_bytes_and_says_what_the_spec_says() -> None:
    blob = _polygon()
    assert blob[0] == gb.MAGIC_GEOMETRY == 0xC8
    assert blob[1] == gb.FORMAT_VERSION == 0x01
    assert (blob[2] >> 4) == 0  # MIN_LEX_CYCLIC_ROTATION
    assert (blob[2] & 0xF) == gb.EDGE_MODELS.index("WGS84_GEODESIC")
    assert (blob[3] >> 4) == gb.DENSIFICATION_MODELS.index("ORTHODROMIC_VINCENTY")
    assert (blob[3] & 0xF) == gb.RESOLUTION_MODES.index("UNIFORM")
    assert ((blob[4] << 8) | blob[5]) == 1000
    assert (blob[6] >> 4) == 3
    assert (blob[6] & 0xF) == gb.GEOMETRY_TYPES.index("POLYGON")
    assert ((blob[7] << 8) | blob[8]) == 1
    assert gb.HEADER_SIZE_BYTES == 9


def test_geometry_type_codes_start_at_one() -> None:
    """Code zero is reserved, so POINT is 1 and MULTIPOLYGON is 6.

    Numbering the enum from zero shifts every type by one. The shift is
    silent: the blob still decodes, into the wrong type, and its identity
    hash still verifies against itself.
    """
    assert gb.GEOMETRY_TYPES[0] == "RESERVED"
    assert gb.GEOMETRY_TYPES.index("POINT") == 1
    assert gb.GEOMETRY_TYPES.index("MULTIPOLYGON") == 6
    for name in ("POINT", "LINESTRING", "POLYGON"):
        rings = {"POINT": [POINT], "LINESTRING": [LINE], "POLYGON": [RING]}[name]
        blob = gb.encode_geometry(rings, name, resolution=3)
        assert blob[6] & 0xF == gb.GEOMETRY_TYPES.index(name)
        assert gb.read_geometry_type(blob) == name


def test_the_reserved_type_code_is_refused() -> None:
    with pytest.raises(UnsupportedGeometryTypeError, match="reserved"):
        gb.encode_geometry([RING], "RESERVED", resolution=3)


def test_an_unknown_type_is_refused() -> None:
    with pytest.raises(UnsupportedGeometryTypeError, match="unknown geometry type"):
        gb.encode_geometry([RING], "TESSERACT", resolution=3)


# --------------------------------------------------------------------------
# Round trip
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("rings", "geometry_type"),
    [
        ([POINT], "POINT"),
        ([LINE], "LINESTRING"),
        ([RING], "POLYGON"),
        ([POINT, ["NE(0625/0451(2(A1)))"]], "MULTIPOINT"),
        ([LINE, list(reversed(LINE))], "MULTILINESTRING"),
    ],
)
def test_round_trip_preserves_rings_and_profile(
    rings: list[list[str]], geometry_type: str
) -> None:
    blob = gb.encode_geometry(rings, geometry_type, resolution=3, canonicalize=False)
    decoded = gb.decode_geometry(blob)
    assert decoded["rings"] == rings
    assert decoded["geometry_type"] == geometry_type
    assert decoded["edge_model"] == "WGS84_GEODESIC"
    assert decoded["densification_model"] == "ORTHODROMIC_VINCENTY"
    assert decoded["max_segment_m"] == 1000
    assert decoded["resolution"] == 3


def test_round_trip_preserves_vertex_order() -> None:
    """Order is the thing GeometryBlob keeps and TreeBlob throws away."""
    forward = gb.encode_geometry([RING], "POLYGON", resolution=3, canonicalize=False)
    backward = gb.encode_geometry(
        [list(reversed(RING))], "POLYGON", resolution=3, canonicalize=False
    )
    assert forward != backward
    assert gb.decode_geometry(forward)["rings"][0] == RING
    assert gb.decode_geometry(backward)["rings"][0] == list(reversed(RING))


def test_hole_order_is_preserved() -> None:
    hole = ["NE(0625/0451(3(A1)))", "NE(0625/0451(3(A2)))", "NE(0625/0451(3(A3)))"]
    blob = gb.encode_geometry([RING, hole], "POLYGON", resolution=3, canonicalize=False)
    assert gb.decode_geometry(blob)["rings"] == [RING, hole]


def test_the_edge_model_survives_the_round_trip() -> None:
    blob = gb.encode_geometry(
        [RING],
        "POLYGON",
        edge_model="PLANAR_SINUSOIDAL_STRAIGHT",
        densification_model="ASSUMED_BY_PRODUCER",
        max_segment_m=250,
        resolution=3,
    )
    decoded = gb.decode_geometry(blob)
    assert decoded["edge_model"] == "PLANAR_SINUSOIDAL_STRAIGHT"
    assert decoded["densification_model"] == "ASSUMED_BY_PRODUCER"
    assert decoded["max_segment_m"] == 250


# --------------------------------------------------------------------------
# The compatibility matrix
# --------------------------------------------------------------------------


def test_a_geodesic_edge_needs_densification() -> None:
    with pytest.raises(IncompatibleProfileError, match="other than NONE"):
        gb.encode_geometry([RING], "POLYGON", densification_model="NONE", resolution=3)


def test_a_straight_edge_refuses_vincenty_densification() -> None:
    with pytest.raises(IncompatibleProfileError, match="contradicts"):
        gb.encode_geometry(
            [RING],
            "POLYGON",
            edge_model="PLANAR_SINUSOIDAL_STRAIGHT",
            densification_model="ORTHODROMIC_VINCENTY",
            resolution=3,
        )


# --------------------------------------------------------------------------
# Structural preconditions the specification states
# --------------------------------------------------------------------------


def test_ring_closure_is_implicit_and_explicit_closure_is_refused() -> None:
    """The last vertex does not repeat the first.

    Accepting a closed ring and then rotating it to canonical form moves
    the repeated vertex into the middle, which turns a closed ring into
    one that is neither closed nor valid.
    """
    with pytest.raises(GeometryError, match="closure is implicit"):
        gb.encode_geometry([RING + [RING[0]]], "POLYGON", resolution=3)


def test_a_ring_needs_three_vertices() -> None:
    with pytest.raises(GeometryError, match="at least 3 vertices"):
        gb.encode_geometry([RING[:2]], "POLYGON", resolution=3)


def test_consecutive_duplicates_are_refused() -> None:
    with pytest.raises(GeometryError, match="consecutive duplicate"):
        gb.encode_geometry([[RING[0], RING[0], RING[1]]], "POLYGON", resolution=3)


def test_a_linestring_needs_two_vertices() -> None:
    with pytest.raises(GeometryError, match="at least 2 vertices"):
        gb.encode_geometry([POINT], "LINESTRING", resolution=3)


def test_a_point_holds_exactly_one_vertex() -> None:
    with pytest.raises(GeometryError, match="exactly one vertex"):
        gb.encode_geometry([LINE], "POINT", resolution=3)


def test_an_empty_geometry_is_refused() -> None:
    with pytest.raises(GeometryError, match="no rings"):
        gb.encode_geometry([], "POLYGON", resolution=3)


def test_an_empty_ring_is_refused() -> None:
    with pytest.raises(GeometryError, match="empty ring"):
        gb.encode_geometry([[]], "POLYGON", resolution=3)


def test_a_single_part_type_refuses_extra_parts() -> None:
    with pytest.raises(GeometryError, match="single part"):
        gb.encode_geometry([LINE, LINE], "LINESTRING", resolution=3)


HOLE = ["NE(0625/0451(3(A1)))", "NE(0625/0451(3(A2)))", "NE(0625/0451(3(A3)))"]
OTHER = ["NE(0626/0451(2(A1)))", "NE(0626/0451(1(A2)))", "NE(0626/0451(1(A1)))"]


def test_multipolygon_nests_one_level_deeper() -> None:
    """Each entry is a polygon, and each polygon is a list of rings."""
    blob = gb.encode_geometry(
        [[RING, HOLE], [OTHER]], "MULTIPOLYGON", resolution=3, canonicalize=False
    )
    assert gb.read_geometry_type(blob) == "MULTIPOLYGON"
    assert ((blob[7] << 8) | blob[8]) == 2
    assert gb.decode_geometry(blob)["rings"] == [RING, HOLE, OTHER]


def test_a_flat_list_is_refused_for_multipolygon() -> None:
    """Attaching every hole to the first polygon would decode cleanly.

    It would also be the wrong geometry, and nothing downstream could
    tell. Refusing the ambiguous spelling is the only reading that does
    not invent an answer.
    """
    with pytest.raises(GeometryError, match="which hole belongs to which"):
        gb.encode_geometry([RING, HOLE], "MULTIPOLYGON", resolution=3)


def test_multipolygon_hole_grouping_survives_the_round_trip() -> None:
    """Two spellings of the same rings are two different geometries."""
    together = gb.encode_geometry(
        [[RING, HOLE], [OTHER]], "MULTIPOLYGON", resolution=3, canonicalize=False
    )
    apart = gb.encode_geometry(
        [[RING], [HOLE], [OTHER]], "MULTIPOLYGON", resolution=3, canonicalize=False
    )
    assert together != apart
    assert gb.geometry_hash(together) != gb.geometry_hash(apart)
    assert gb.geometry_to_tree(together) == gb.geometry_to_tree(apart)


def test_vertices_must_all_sit_at_the_declared_resolution() -> None:
    mixed = [RING[0], RING[1], "NE(0625/0451(1))"]
    with pytest.raises(GeometryError, match="not 3"):
        gb.encode_geometry([mixed], "POLYGON", resolution=3)


@pytest.mark.parametrize(
    "branching",
    ["NE(0625/0451(1(A1,A2)))", "NE(0625/0451(1(A1))),NW(0625/0451(1(A1)))"],
)
def test_a_branching_index_is_not_a_vertex(branching: str) -> None:
    """One vertex is one cell, whether the branch is inside or across."""
    with pytest.raises(GeometryError, match="one cell, not a tree"):
        gb.encode_geometry([[branching, RING[0], RING[1]]], "POLYGON", resolution=3)


@pytest.mark.parametrize("bad", ["", "   ", "not an index", 7])
def test_a_malformed_vertex_is_refused(bad: object) -> None:
    with pytest.raises(GeometryError):
        gb.encode_geometry(
            [[bad, RING[0], RING[1]]],  # type: ignore[list-item]
            "POLYGON",
            resolution=3,
        )


def test_an_out_of_range_resolution_is_refused() -> None:
    with pytest.raises(GeometryError, match="resolution 14"):
        gb.encode_geometry([RING], "POLYGON", resolution=14)


def test_an_out_of_range_segment_bound_is_refused() -> None:
    with pytest.raises(GeometryError, match="max_segment_m"):
        gb.encode_geometry([RING], "POLYGON", max_segment_m=70000, resolution=3)


# --------------------------------------------------------------------------
# Canonicalisation and identity
# --------------------------------------------------------------------------


def test_the_hash_is_invariant_under_ring_rotation() -> None:
    """What makes the blob content-addressable."""
    digests = set()
    for shift in range(len(RING)):
        rotated = RING[shift:] + RING[:shift]
        digests.add(gb.geometry_hash(_polygon([rotated])))
    assert len(digests) == 1


def test_without_canonicalisation_rotation_changes_the_hash() -> None:
    a = _polygon([RING], canonicalize=False)
    b = _polygon([RING[1:] + RING[:1]], canonicalize=False)
    assert gb.geometry_hash(a) != gb.geometry_hash(b)


def test_the_hash_is_keccak_and_not_sha3() -> None:
    """The two differ in one padding byte, and the identity is on-chain.

    Substituting :mod:`hashlib` here changes every digest silently.
    """
    import hashlib

    blob = _polygon()
    assert len(gb.geometry_hash(blob)) == 32
    assert gb.geometry_hash(blob) != hashlib.sha3_256(blob).digest()
    assert gb._keccak256(b"") == bytes.fromhex(
        "c5d2460186f7233c927e7db2dcc703c0e500b653ca82273b7bfad8045d85a470"
    )


def test_hashing_refuses_a_blob_that_is_not_one() -> None:
    with pytest.raises(MalformedBlobError):
        gb.geometry_hash(b"\x00" * 9)


# --------------------------------------------------------------------------
# The bridge to TreeBlob, one-way by construction
# --------------------------------------------------------------------------


def test_geometry_to_tree_loses_order() -> None:
    """The derived TreeBlob is a function of the vertex set alone.

    Order, ring topology and edge model all vanish. This is the property
    the format is built on, not a defect: coverage survives, identity
    does not.
    """
    forward = _polygon([RING], canonicalize=False)
    backward = _polygon([list(reversed(RING))], canonicalize=False)
    assert forward != backward
    assert gb.geometry_to_tree(forward) == gb.geometry_to_tree(backward)


def test_geometry_to_tree_loses_the_edge_model() -> None:
    geodesic = _polygon([RING])
    straight = gb.encode_geometry(
        [RING],
        "POLYGON",
        edge_model="PLANAR_SINUSOIDAL_STRAIGHT",
        densification_model="ASSUMED_BY_PRODUCER",
        resolution=3,
    )
    assert geodesic != straight
    assert gb.geometry_to_tree(geodesic) == gb.geometry_to_tree(straight)


def test_geometry_to_tree_loses_the_geometry_type() -> None:
    as_polygon = _polygon([RING])
    as_line = gb.encode_geometry([RING], "LINESTRING", resolution=3)
    assert gb.geometry_to_tree(as_polygon) == gb.geometry_to_tree(as_line)


def test_geometry_to_tree_agrees_with_encoding_the_vertex_set() -> None:
    blob = _polygon([RING])
    cells = sorted(set(gb.decode_geometry(blob)["rings"][0]))
    assert gb.geometry_to_tree(blob) == tb.encode_tree(itacart.compose(cells))


def test_geometry_to_tree_deduplicates_repeated_vertices() -> None:
    hole = [RING[0], "NE(0625/0451(3(A2)))", "NE(0625/0451(3(A3)))"]
    blob = gb.encode_geometry([RING, hole], "POLYGON", resolution=3)
    tree = gb.geometry_to_tree(blob)
    assert tb.count_vertices(tree) == len(set(RING) | set(hole))


def test_geometry_to_tree_carries_vertices_from_two_quadrants() -> None:
    """A geometry may cross the equator, and its tree has to follow.

    A TreeBlob holds one group per quadrant, so a region straddling a
    quadrant boundary is an ordinary blob rather than a refusal.
    """
    mixed = [RING[0], RING[1], "NW(0625/0451(1(A1)))"]
    blob = gb.encode_geometry([mixed], "POLYGON", resolution=3, canonicalize=False)
    tree = gb.geometry_to_tree(blob)
    assert tb.count_vertices(tree) == 3
    assert tb.decode_tree(tree).count("(") > 0
    assert sorted({cell[:2] for cell in itacart.decompose(tb.decode_tree(tree))}) == [
        "NE",
        "NW",
    ]


# --------------------------------------------------------------------------
# Malformed blobs
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        ((0, 0x00), "bad magic"),
        ((1, 0x02), "unsupported format version"),
        ((2, 0x10), "unknown canonical profile"),
        ((2, 0x0F), "unknown edge model"),
        ((3, 0xF0), "unknown densification model"),
        ((3, 0x0F), "unknown resolution mode"),
        ((6, 0xDF), "unknown geometry type"),
        ((6, 0x00), "resolution 0 out of range"),
    ],
)
def test_a_corrupt_header_is_refused(mutate: tuple[int, int], message: str) -> None:
    blob = bytearray(_polygon())
    index, value = mutate
    blob[index] = value
    with pytest.raises(MalformedBlobError, match=message):
        gb.validate_geometry(bytes(blob))


def test_a_short_blob_is_refused() -> None:
    with pytest.raises(MalformedBlobError, match="shorter than"):
        gb.validate_geometry(_polygon()[:5])


def test_a_non_bytes_blob_is_refused() -> None:
    with pytest.raises(MalformedBlobError, match="must be bytes"):
        gb.validate_geometry("not bytes")  # type: ignore[arg-type]


def test_a_multi_resolution_blob_is_refused() -> None:
    blob = bytearray(_polygon())
    blob[3] = (blob[3] & 0xF0) | gb.RESOLUTION_MODES.index("MULTI")
    with pytest.raises(MalformedBlobError, match="not written here"):
        gb.validate_geometry(bytes(blob))


def test_reserved_part_bits_must_be_zero() -> None:
    blob = bytearray(_polygon())
    blob[11] |= 0x0F
    with pytest.raises(MalformedBlobError, match="reserved part bits"):
        gb.validate_geometry(bytes(blob))


def test_trailing_bytes_after_the_stream_are_refused() -> None:
    with pytest.raises(MalformedBlobError, match="unread bits"):
        gb.validate_geometry(_polygon() + b"\x00\x00")


def test_trailing_padding_must_be_zero() -> None:
    blob = bytearray(_polygon())
    blob[-1] |= 0x01
    with pytest.raises(MalformedBlobError, match="padding is not zero"):
        gb.validate_geometry(bytes(blob))


def test_a_truncated_vertex_stream_is_refused() -> None:
    with pytest.raises(MalformedBlobError, match="ends inside a field"):
        gb.validate_geometry(_polygon()[:12])


def test_a_vertex_that_is_not_a_string_is_refused() -> None:
    with pytest.raises(GeometryError, match="non-empty index string"):
        gb.encode_geometry(
            [[7, RING[0], RING[1]]],  # type: ignore[list-item]
            "POLYGON",
            resolution=3,
        )


def test_a_quadrant_without_a_cell_is_not_a_vertex() -> None:
    with pytest.raises(GeometryError, match="malformed vertex index"):
        gb.encode_geometry([["NE()", RING[0], RING[1]]], "POLYGON", resolution=3)


def test_a_corrupt_column_in_a_vertex_stream_is_refused() -> None:
    """Eleven bits hold 2047; the grammar admits 2004."""
    blob = bytearray(_polygon())
    blob[13] |= 0b0011_1111
    blob[14] |= 0b1111_1111
    with pytest.raises(MalformedBlobError, match="not addressable"):
        gb.validate_geometry(bytes(blob))


def test_an_impossible_shared_level_count_is_refused() -> None:
    """A vertex cannot share more levels than the resolution has.

    The field is four bits wide and holds fifteen; a resolution-3 blob
    has three levels to share.
    """
    blob = bytearray(
        gb.encode_geometry([RING], "POLYGON", resolution=3, canonicalize=False)
    )
    blob[17] |= 0b0000_0011  # first two bits of the shared field, at bit 142
    blob[18] |= 0b1100_0000  # its remaining two bits
    with pytest.raises(MalformedBlobError, match="claims 15 shared levels"):
        gb.validate_geometry(bytes(blob))


def test_a_reserved_refinement_code_in_a_vertex_is_refused() -> None:
    """Five bits hold 31; a resolution-3 level addresses 25 tokens.

    Codes 25 to 31 name nothing, so a blob carrying one is malformed
    however well the rest of it reads.
    """
    blob = bytearray(
        gb.encode_geometry([RING], "POLYGON", resolution=3, canonicalize=False)
    )
    blob[17] |= 0b0111_1100  # the five-bit resolution-3 field of vertex 0
    with pytest.raises(MalformedBlobError, match="reserved at resolution 3"):
        gb.validate_geometry(bytes(blob))
