"""GeometryBlob: binary encoding of an ordered ITACaRT vertex sequence.

Preserves what TreeBlob discards: vertex order, ring topology, edge model,
densification parameters and OGC SFA type. That makes it the form to hash
when a geometry must be identified rather than merely covered.

Header layout, 9 bytes:

=================  ====  =======================  ============================
Byte(s)            Bits  Field                    Values
=================  ====  =======================  ============================
``0x00``           8     Magic                    ``0xC8``
``0x01``           8     Format version           ``0x01``
``0x02`` hi        4     Canonical profile        ``0x0`` MIN_LEX_CYCLIC_ROTATION
``0x02`` lo        4     Edge model               ``0x0`` PLANAR_SINUSOIDAL_STRAIGHT,
                                                  ``0x1`` WGS84_GEODESIC
``0x03`` hi        4     Densification model      ``0x0`` NONE,
                                                  ``0x1`` ORTHODROMIC_VINCENTY,
                                                  ``0x2`` ASSUMED_BY_PRODUCER
``0x03`` lo        4     Resolution mode          ``0x0`` UNIFORM, ``0x1`` MULTI
``0x04``-``0x05``  16    ``max_segment_m``        ``0x03E8`` = 1000 m
``0x06`` hi        4     Uniform resolution       ``0xD`` = res 13, ``0x0`` if MULTI
``0x06`` lo        4     Geometry type (OGC SFA)  ``0x0`` reserved, ``0x1`` POINT
                                                  .. ``0x6`` MULTIPOLYGON
``0x07``-``0x08``  16    ``component_count``      1 for single, N for MULTI*
=================  ====  =======================  ============================

Components, parts and rings follow, then the vertex stream. Vertices are
delta-encoded by ``shared_levels``, the count of hierarchy levels a vertex
shares with its predecessor; vertex 0 is always written in full.

Density is a curve in the depth of the shared prefix, not a constant. A
vertex sharing nothing with its predecessor costs 65 bits at resolution 13
plus a 4-bit ``shared_levels`` field; one sharing all thirteen levels costs
the 4-bit field alone. Measured on a synthetic resolution-13 polygon whose
vertices share levels 1 to 7, the amortised cost is about 25 bits per vertex
for a uniform blob and about 29 for a multi-resolution one, easing towards
those figures as vertex count rises.
"""

from __future__ import annotations

from typing import Any, Sequence

from ..exceptions import (
    GeometryError,
    IncompatibleProfileError,
    MalformedBlobError,
    UnsupportedGeometryTypeError,
)
from .tree_blob import (
    _QUADRANT_TO_BITS,
    MAX_RESOLUTION,
    MIN_RESOLUTION,
    X_BITS,
    Y_BITS,
    _alphabet,
    _BitReader,
    _BitWriter,
    _column_is_addressable,
    _component_width,
    encode_tree,
)

__all__ = [
    "encode_geometry",
    "decode_geometry",
    "validate_geometry",
    "geometry_hash",
    "geometry_to_tree",
    "read_geometry_type",
]

MAGIC_GEOMETRY = 0xC8
FORMAT_VERSION = 0x01
HEADER_SIZE_BYTES = 9

_BITS_TO_QUADRANT = {bits: name for name, bits in _QUADRANT_TO_BITS.items()}

#: OGC Simple Feature Access type codes. Code 0 is reserved, so POINT is 1
#: and MULTIPOLYGON is 6. Numbering these from zero shifts every type by
#: one and silently changes the identity hash of every blob.
GEOMETRY_TYPES = (
    "RESERVED",
    "POINT",
    "LINESTRING",
    "POLYGON",
    "MULTIPOINT",
    "MULTILINESTRING",
    "MULTIPOLYGON",
)
EDGE_MODELS = ("PLANAR_SINUSOIDAL_STRAIGHT", "WGS84_GEODESIC")
DENSIFICATION_MODELS = ("NONE", "ORTHODROMIC_VINCENTY", "ASSUMED_BY_PRODUCER")
RESOLUTION_MODES = ("UNIFORM", "MULTI")
CANONICAL_PROFILES = ("MIN_LEX_CYCLIC_ROTATION",)

#: Ring roles, written in the high nibble of each part header.
RING_ROLES = ("EXTERIOR", "INTERIOR", "LINE", "POINT")

_MULTI_TYPES = {"MULTIPOINT", "MULTILINESTRING", "MULTIPOLYGON"}
_SINGLE_VERTEX_TYPES = {"POINT", "MULTIPOINT"}
_CLOSED_TYPES = {"POLYGON", "MULTIPOLYGON"}


# ---------------------------------------------------------------------------
# Vertices
# ---------------------------------------------------------------------------


class _Vertex:
    """A parsed atomic index, decomposed for bit-level work."""

    __slots__ = ("quadrant", "column", "row", "refinements")

    def __init__(
        self, quadrant: str, column: int, row: int, refinements: "tuple[str, ...]"
    ) -> None:
        self.quadrant = quadrant
        self.column = column
        self.row = row
        self.refinements = refinements

    @property
    def resolution(self) -> int:
        return 1 + len(self.refinements)

    def to_index(self) -> str:
        inner = ""
        for token in reversed(self.refinements):
            inner = f"({token}{inner})"
        return f"{self.quadrant}({self.column:04d}/{self.row:04d}{inner})"


def _parse_vertex(index: str) -> _Vertex:
    """Parse one atomic index. Branching indices belong to TreeBlob."""
    from .tree_blob import _merge, _pair, _parse

    try:
        globe = _merge(_parse(index.strip()))
    except Exception as exc:
        raise GeometryError(f"malformed vertex index {index!r}: {exc}") from exc
    if len(globe.children) > 1:
        raise GeometryError(f"a vertex is one cell, not a tree: {index!r}")
    root = globe.children[0]
    quadrant = str(root.component)
    node = root
    refinements: list[str] = []
    column = row = -1
    while node.children:
        if len(node.children) > 1:
            raise GeometryError(f"a vertex is one cell, not a tree: {index!r}")
        node = node.children[0]
        if node.resolution == 1:
            column, row = _pair(node.component)
        else:
            refinements.append(str(node.component))
    return _Vertex(quadrant, column, row, tuple(refinements))


def _write_vertex(writer: _BitWriter, vertex: _Vertex) -> None:
    """Write a vertex in full: quadrant, column, row, then refinements.

    The column, the row and every refinement were settled when the index
    was parsed. Checking them again here would put a second copy of each
    rule downstream of the one that decides, and a copy that cannot fire
    is a copy that cannot be kept honest.
    """
    writer.write(_QUADRANT_TO_BITS[vertex.quadrant], 2)
    writer.write(vertex.column, X_BITS)
    writer.write(vertex.row, Y_BITS)
    _write_refinements(writer, vertex, 0)


def _write_refinements(writer: _BitWriter, vertex: _Vertex, start: int) -> None:
    for level, token in enumerate(vertex.refinements):
        if level < start:
            continue
        resolution = level + 2
        writer.write(_alphabet(resolution).index(token), _component_width(resolution))


def _read_vertex(reader: _BitReader, resolution: int) -> _Vertex:
    """Read a vertex written in full. The header settled the resolution."""
    quadrant = _BITS_TO_QUADRANT[reader.read(2)]
    column = reader.read(X_BITS)
    row = reader.read(Y_BITS)
    refinements = _read_refinements(reader, resolution, [], 0)
    if not _column_is_addressable(column, bool(refinements), quadrant):
        raise MalformedBlobError(f"column {column} is not addressable in {quadrant}")
    return _Vertex(quadrant, column, row, tuple(refinements))


def _read_refinements(
    reader: _BitReader, resolution: int, prefix: "Sequence[str]", start: int
) -> "list[str]":
    tokens = list(prefix[:start])
    for level in range(start, resolution - 1):
        target = level + 2
        alphabet = _alphabet(target)
        code = reader.read(_component_width(target))
        if code >= len(alphabet):
            raise MalformedBlobError(
                f"refinement code {code} is reserved at resolution {target}"
            )
        tokens.append(alphabet[code])
    return tokens


def _shared_levels(previous: _Vertex, current: _Vertex) -> int:
    """Hierarchy levels the two vertices have in common."""
    if (
        previous.quadrant != current.quadrant
        or previous.column != current.column
        or previous.row != current.row
    ):
        return 0
    shared = 1
    # Falling out of this loop without breaking would mean two identical
    # consecutive vertices, which `_validate_part` refuses, and every vertex
    # here sits at the blob's uniform resolution, so the two refinement
    # tuples are always the same length. The loop therefore always breaks.
    for a, b in zip(  # pragma: no cover - the exhausting arc cannot be taken
        previous.refinements, current.refinements
    ):
        if a != b:
            break
        shared += 1
    return shared


# ---------------------------------------------------------------------------
# Enumerations
# ---------------------------------------------------------------------------


def _code_of(name: str, table: "Sequence[str]", what: str) -> int:
    try:
        return table.index(name.upper())
    except (AttributeError, ValueError):
        raise UnsupportedGeometryTypeError(
            f"unknown {what}: {name!r}; expected one of {', '.join(table)}"
        ) from None


def _geometry_type_code(name: str) -> int:
    code = _code_of(name, GEOMETRY_TYPES, "geometry type")
    if code == 0:
        raise UnsupportedGeometryTypeError("geometry type code 0 is reserved")
    return code


def _check_profile(edge_model: str, densification_model: str) -> None:
    """Refuse combinations the format cannot mean.

    A geodesic edge with no densification claims a curve was preserved by
    a producer that ran nothing, and a straight sinusoidal edge with
    Vincenty densification claims the opposite.
    """
    if edge_model == "WGS84_GEODESIC" and densification_model == "NONE":
        raise IncompatibleProfileError(
            "WGS84_GEODESIC edges need a densification model other than NONE"
        )
    if (
        edge_model == "PLANAR_SINUSOIDAL_STRAIGHT"
        and densification_model == "ORTHODROMIC_VINCENTY"
    ):
        raise IncompatibleProfileError(
            "PLANAR_SINUSOIDAL_STRAIGHT edges are straight; "
            "ORTHODROMIC_VINCENTY densification contradicts them"
        )


# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------


def _components_of(
    rings: "Sequence[Any]", geometry_type: str
) -> "list[list[list[str]]]":
    """Group the caller's rings into the components and parts written out.

    A POLYGON is one component whose parts are its exterior and its
    holes; MULTIPOINT and MULTILINESTRING are one component per entry.
    MULTIPOLYGON is the one type a flat list cannot express, because a
    flat list cannot say which hole belongs to which polygon, so it is
    nested one level deeper: a list of polygons, each a list of rings.
    """
    if not rings:
        raise GeometryError("geometry has no rings")
    if geometry_type == "MULTIPOLYGON":
        components = [[list(ring) for ring in polygon] for polygon in _polygons(rings)]
    elif geometry_type == "POLYGON":
        components = [[list(ring) for ring in rings]]
    else:
        if geometry_type in ("POINT", "LINESTRING") and len(rings) != 1:
            raise GeometryError(f"{geometry_type} is a single part, got {len(rings)}")
        components = [[list(ring)] for ring in rings]
    for component in components:
        for ring in component:
            if not ring:
                raise GeometryError("geometry has an empty ring")
            _validate_part(ring, geometry_type)
    return components


def _polygons(rings: "Sequence[Any]") -> "list[Sequence[Sequence[str]]]":
    """Read a MULTIPOLYGON argument as a list of polygons.

    Each entry must itself be a list of rings. A caller who passes the
    flat form a POLYGON takes gets told so, rather than having its rings
    silently attached to one polygon.
    """
    polygons: "list[Sequence[Sequence[str]]]" = []
    for polygon in rings:
        if isinstance(polygon, str) or not all(
            not isinstance(ring, str) for ring in polygon
        ):
            raise GeometryError(
                "MULTIPOLYGON takes a list of polygons, each a list of rings; "
                "a flat list of rings cannot say which hole belongs to which"
            )
        polygons.append(polygon)
    return polygons


def _validate_part(ring: "Sequence[str]", geometry_type: str) -> None:
    """Enforce the structural preconditions the format states.

    Ring closure is implicit: the last vertex does not repeat the first,
    and an explicitly closed ring is refused rather than silently
    accepted. Accepting one and then rotating it to canonical form moves
    the repeated vertex into the middle of the sequence, which is how a
    closed ring becomes a ring that is neither closed nor valid.
    """
    for vertex in ring:
        if not isinstance(vertex, str) or not vertex.strip():
            raise GeometryError(f"vertex must be a non-empty index string: {vertex!r}")
    for previous, current in zip(ring, ring[1:]):
        if previous == current:
            raise GeometryError(
                f"{geometry_type} has consecutive duplicate vertices: {current}"
            )
    if geometry_type in _SINGLE_VERTEX_TYPES:
        if len(ring) != 1:
            raise GeometryError(
                f"{geometry_type} parts hold exactly one vertex, got {len(ring)}"
            )
        return
    if geometry_type in _CLOSED_TYPES:
        if len(ring) < 3:
            raise GeometryError(
                f"{geometry_type} rings need at least 3 vertices, got {len(ring)}"
            )
        if ring[0] == ring[-1]:
            raise GeometryError(
                "ring closure is implicit; drop the repeated closing vertex"
            )
        return
    if len(ring) < 2:
        raise GeometryError(
            f"{geometry_type} needs at least 2 vertices, got {len(ring)}"
        )


def _ring_role(part_index: int, geometry_type: str) -> int:
    if geometry_type in _SINGLE_VERTEX_TYPES:
        return RING_ROLES.index("POINT")
    if geometry_type in _CLOSED_TYPES:
        return RING_ROLES.index("EXTERIOR" if part_index == 0 else "INTERIOR")
    return RING_ROLES.index("LINE")


# ---------------------------------------------------------------------------
# Encoding
# ---------------------------------------------------------------------------


def encode_geometry(
    rings: "Sequence[Any]",
    geometry_type: str = "POLYGON",
    edge_model: str = "WGS84_GEODESIC",
    densification_model: str = "ORTHODROMIC_VINCENTY",
    max_segment_m: int = 1000,
    resolution: int = 13,
    canonicalize: bool = True,
) -> bytes:
    """Encode ordered vertex rings as a GeometryBlob.

    Ring orientation is significant and is preserved exactly.
    Canonicalisation may rotate a closed ring to its least starting
    vertex; it never reverses one, because rotating is a change of
    spelling and reversing is a change of meaning.

    Which ring is the exterior and which are holes is carried here by
    **part order** — the first ring of a component is its exterior — and
    not by reading the winding. Winding is preserved to the byte all the
    same, because many producers and consumers use direction to carry
    that same distinction, and a ring handed on reversed may be read
    downstream as the opposite kind of ring. The common convention is
    counter-clockwise for an exterior boundary and clockwise for a hole::

        exterior = [A, B, C, D]   # counter-clockwise
        hole     = [E, H, G, F]   # clockwise: E, F, G, H reversed

    Reversing either sequence changes its identity: the bytes differ and
    so does the hash. Nothing here reverses a ring to make it agree with
    a convention, and nothing here checks that it does.

    Args:
        rings: Rings as lists of atomic index strings, exterior first.
            MULTIPOLYGON nests one level deeper, as a list of polygons,
            because a flat list cannot say which hole belongs to which
            polygon.
        geometry_type: One of POINT, LINESTRING, POLYGON, MULTIPOINT,
            MULTILINESTRING, MULTIPOLYGON.
        edge_model: How edges between vertices are interpreted.
        densification_model: Whether the producer ran Vincenty
            densification or merely guarantees the segment bound.
        max_segment_m: Densification bound in metres, 0 to 65 535.
        resolution: Uniform resolution of every vertex cell.
        canonicalize: Apply minimum lexicographic cyclic rotation to
            closed rings. Directional types are never rotated.

    Returns:
        The binary blob.

    Raises:
        UnsupportedGeometryTypeError: On unsupported types.
        GeometryError: If ring structure contradicts the declared type.
        IncompatibleProfileError: On an invalid edge model and
            densification model combination.

    Example:
        >>> blob = encode_geometry([["NE(0625/0451)"]], "POINT", resolution=1)
        >>> read_geometry_type(blob)
        'POINT'

        Rotation is a spelling; reversal is not:

        >>> ring = ["NE(0625/0451(1))", "NE(0625/0451(2))", "NE(0625/0451(3))"]
        >>> encode_geometry([ring], resolution=2) == encode_geometry(
        ...     [ring[1:] + ring[:1]], resolution=2
        ... )
        True
        >>> encode_geometry([ring], resolution=2) == encode_geometry(
        ...     [list(reversed(ring))], resolution=2
        ... )
        False
    """
    from ..geometry import canonicalize_rings

    geometry_type = geometry_type.upper()
    type_code = _geometry_type_code(geometry_type)
    edge_code = _code_of(edge_model, EDGE_MODELS, "edge model")
    densification_code = _code_of(
        densification_model, DENSIFICATION_MODELS, "densification model"
    )
    _check_profile(edge_model.upper(), densification_model.upper())
    if not MIN_RESOLUTION <= resolution <= MAX_RESOLUTION:
        raise GeometryError(f"resolution {resolution} out of range")
    if not 0 <= max_segment_m <= 0xFFFF:
        raise GeometryError(f"max_segment_m {max_segment_m} outside [0, 65535]")

    components = _components_of(rings, geometry_type)
    if canonicalize and geometry_type in _CLOSED_TYPES:
        components = [canonicalize_rings(part) for part in components]
    if len(components) > 0xFFFF:  # pragma: no cover - 65 536 components
        raise GeometryError(f"too many components: {len(components)}")

    parsed: "list[list[list[_Vertex]]]" = []
    for component in components:
        parsed_parts: "list[list[_Vertex]]" = []
        for part in component:
            vertices = [_parse_vertex(index) for index in part]
            for vertex in vertices:
                if vertex.resolution != resolution:
                    raise GeometryError(
                        f"vertex {vertex.to_index()} is at resolution "
                        f"{vertex.resolution}, not {resolution}"
                    )
            parsed_parts.append(vertices)
        parsed.append(parsed_parts)

    writer = _BitWriter()
    writer.write(MAGIC_GEOMETRY, 8)
    writer.write(FORMAT_VERSION, 8)
    writer.write(0, 4)  # canonical profile: MIN_LEX_CYCLIC_ROTATION
    writer.write(edge_code, 4)
    writer.write(densification_code, 4)
    writer.write(0, 4)  # resolution mode: UNIFORM
    writer.write(max_segment_m, 16)
    writer.write(resolution, 4)
    writer.write(type_code, 4)
    writer.write(len(components), 16)

    for parsed_component in parsed:
        if len(parsed_component) > 0xFFFF:  # pragma: no cover - 65 536 parts
            raise GeometryError(f"too many parts: {len(parsed_component)}")
        writer.write(len(parsed_component), 16)
        for part_index, vertices in enumerate(parsed_component):
            writer.write(_ring_role(part_index, geometry_type), 4)
            writer.write(0, 4)  # reserved per-part nibble
            if len(vertices) > 0xFFFF:  # pragma: no cover - 65 536 vertices
                raise GeometryError(f"too many vertices: {len(vertices)}")
            writer.write(len(vertices), 16)
            _write_stream(writer, vertices)
    return writer.to_bytes()


def _write_stream(writer: _BitWriter, vertices: "Sequence[_Vertex]") -> None:
    previous: "_Vertex | None" = None
    for vertex in vertices:
        if previous is None:
            _write_vertex(writer, vertex)
        else:
            shared = _shared_levels(previous, vertex)
            writer.write(shared, 4)
            if shared == 0:
                _write_vertex(writer, vertex)
            else:
                _write_refinements(writer, vertex, shared - 1)
        previous = vertex


# ---------------------------------------------------------------------------
# Decoding
# ---------------------------------------------------------------------------


def _read_header(blob: bytes) -> "dict[str, Any]":
    if not isinstance(blob, (bytes, bytearray)):
        raise MalformedBlobError(f"blob must be bytes, got {type(blob).__name__}")
    if len(blob) < HEADER_SIZE_BYTES:
        raise MalformedBlobError("blob is shorter than the geometry header")
    if blob[0] != MAGIC_GEOMETRY:
        raise MalformedBlobError(f"bad magic byte 0x{blob[0]:02X}")
    if blob[1] != FORMAT_VERSION:
        raise MalformedBlobError(f"unsupported format version {blob[1]}")
    profile = (blob[2] >> 4) & 0xF
    if profile >= len(CANONICAL_PROFILES):
        raise MalformedBlobError(f"unknown canonical profile {profile}")
    edge = blob[2] & 0xF
    if edge >= len(EDGE_MODELS):
        raise MalformedBlobError(f"unknown edge model {edge}")
    densification = (blob[3] >> 4) & 0xF
    if densification >= len(DENSIFICATION_MODELS):
        raise MalformedBlobError(f"unknown densification model {densification}")
    mode = blob[3] & 0xF
    if mode >= len(RESOLUTION_MODES):
        raise MalformedBlobError(f"unknown resolution mode {mode}")
    if mode != 0:
        raise MalformedBlobError("multi-resolution blobs are not written here")
    resolution = (blob[6] >> 4) & 0xF
    if not MIN_RESOLUTION <= resolution <= MAX_RESOLUTION:
        raise MalformedBlobError(f"resolution {resolution} out of range")
    type_code = blob[6] & 0xF
    if type_code == 0 or type_code >= len(GEOMETRY_TYPES):
        raise MalformedBlobError(f"unknown geometry type code {type_code}")
    return {
        "canonical_profile": CANONICAL_PROFILES[profile],
        "edge_model": EDGE_MODELS[edge],
        "densification_model": DENSIFICATION_MODELS[densification],
        "resolution_mode": RESOLUTION_MODES[mode],
        "max_segment_m": (blob[4] << 8) | blob[5],
        "resolution": resolution,
        "geometry_type": GEOMETRY_TYPES[type_code],
        "component_count": (blob[7] << 8) | blob[8],
    }


def decode_geometry(blob: bytes) -> "dict[str, Any]":
    """Decode a GeometryBlob into rings and profile metadata.

    Rings come back in the order they were written and with the winding
    they were given. The exterior of a component is its first ring.

    Args:
        blob: The binary blob.

    Returns:
        A mapping with ``rings``, ``geometry_type``, ``edge_model``,
        ``densification_model``, ``max_segment_m`` and ``resolution``.

    Raises:
        MalformedBlobError: If the blob fails structural validation.

    Example:
        >>> blob = encode_geometry([["NE(0625/0451)"]], "POINT", resolution=1)
        >>> decode_geometry(blob)["rings"]
        [['NE(0625/0451)']]
    """
    header = _read_header(blob)
    resolution = int(header["resolution"])
    reader = _BitReader(bytes(blob), start_bit=HEADER_SIZE_BYTES * 8)
    rings: "list[list[str]]" = []
    for _ in range(int(header["component_count"])):
        for _ in range(reader.read(16)):
            reader.read(4)  # ring role, recoverable from the geometry type
            if reader.read(4):
                raise MalformedBlobError("reserved part bits are set")
            count = reader.read(16)
            rings.append(
                [
                    vertex.to_index()
                    for vertex in _read_stream(reader, count, resolution)
                ]
            )
    left = reader.remaining_bits()
    if left >= 8:
        raise MalformedBlobError(f"{left} unread bits after the vertex stream")
    if left and reader.read(left) != 0:
        raise MalformedBlobError("trailing padding is not zero")
    header.pop("component_count")
    header["rings"] = rings
    return header


def _read_stream(reader: _BitReader, count: int, resolution: int) -> "list[_Vertex]":
    vertices: "list[_Vertex]" = []
    previous: "_Vertex | None" = None
    for position in range(count):
        if previous is None:
            vertex = _read_vertex(reader, resolution)
        else:
            shared = reader.read(4)
            if shared > resolution:
                raise MalformedBlobError(
                    f"vertex {position} claims {shared} shared levels "
                    f"at resolution {resolution}"
                )
            if shared == 0:
                vertex = _read_vertex(reader, resolution)
            else:
                tokens = _read_refinements(
                    reader, resolution, previous.refinements, shared - 1
                )
                vertex = _Vertex(
                    previous.quadrant,
                    previous.column,
                    previous.row,
                    tuple(tokens),
                )
        vertices.append(vertex)
        previous = vertex
    return vertices


def validate_geometry(blob: bytes) -> None:
    """Validate blob structure, raising on the first problem found.

    Args:
        blob: The binary blob.

    Raises:
        MalformedBlobError: With a message identifying the failure.

    Example:
        >>> blob = encode_geometry([["NE(0625/0451)"]], "POINT", resolution=1)
        >>> validate_geometry(blob) is None
        True
    """
    decode_geometry(blob)


def read_geometry_type(blob: bytes) -> str:
    """Read the OGC SFA type from the header without full decoding.

    Args:
        blob: The binary blob.

    Returns:
        The geometry type name.

    Raises:
        MalformedBlobError: If the header is unreadable.

    Example:
        >>> blob = encode_geometry([["NE(0625/0451)"]], "POINT", resolution=1)
        >>> read_geometry_type(blob)
        'POINT'
    """
    return str(_read_header(blob)["geometry_type"])


# ---------------------------------------------------------------------------
# Identity
# ---------------------------------------------------------------------------


def _keccak256(data: bytes) -> bytes:
    """Keccak-256, the digest the on-chain identity is keyed on.

    Not SHA3-256: the two differ in their padding byte, so substituting
    :mod:`hashlib` here would silently change every identity.
    """
    mask = 0xFFFFFFFFFFFFFFFF
    rounds = [
        0x0000000000000001,
        0x0000000000008082,
        0x800000000000808A,
        0x8000000080008000,
        0x000000000000808B,
        0x0000000080000001,
        0x8000000080008081,
        0x8000000000008009,
        0x000000000000008A,
        0x0000000000000088,
        0x0000000080008009,
        0x000000008000000A,
        0x000000008000808B,
        0x800000000000008B,
        0x8000000000008089,
        0x8000000000008003,
        0x8000000000008002,
        0x8000000000000080,
        0x000000000000800A,
        0x800000008000000A,
        0x8000000080008081,
        0x8000000000008080,
        0x0000000080000001,
        0x8000000080008008,
    ]
    offsets = [
        [0, 36, 3, 41, 18],
        [1, 44, 10, 45, 2],
        [62, 6, 43, 15, 61],
        [28, 55, 25, 21, 56],
        [27, 20, 39, 8, 14],
    ]

    def rotate(value: int, amount: int) -> int:
        return ((value << amount) | (value >> (64 - amount))) & mask

    def permute(state: "list[list[int]]") -> "list[list[int]]":
        for constant in rounds:
            parity = [
                state[x][0] ^ state[x][1] ^ state[x][2] ^ state[x][3] ^ state[x][4]
                for x in range(5)
            ]
            theta = [
                parity[(x - 1) % 5] ^ rotate(parity[(x + 1) % 5], 1) for x in range(5)
            ]
            for x in range(5):
                for y in range(5):
                    state[x][y] ^= theta[x]
            scratch = [[0] * 5 for _ in range(5)]
            for x in range(5):
                for y in range(5):
                    scratch[y][(2 * x + 3 * y) % 5] = rotate(state[x][y], offsets[x][y])
            for x in range(5):
                for y in range(5):
                    state[x][y] = scratch[x][y] ^ (
                        (~scratch[(x + 1) % 5][y]) & scratch[(x + 2) % 5][y] & mask
                    )
            state[0][0] ^= constant
        return state

    rate = 136
    padded = bytearray(data)
    padded.append(0x01)
    while len(padded) % rate != rate - 1:
        padded.append(0x00)
    padded.append(0x80)

    state = [[0] * 5 for _ in range(5)]
    for start in range(0, len(padded), rate):
        block = padded[start : start + rate]
        for lane in range(rate // 8):
            value = int.from_bytes(block[lane * 8 : lane * 8 + 8], "little")
            state[lane % 5][lane // 5] ^= value
        state = permute(state)

    digest = bytearray()
    for lane in range(4):
        digest.extend(state[lane % 5][lane // 5].to_bytes(8, "little"))
    return bytes(digest[:32])


def geometry_hash(blob: bytes) -> bytes:
    """Content hash of a canonical GeometryBlob.

    Meaningful only for canonicalised blobs: without canonicalisation the
    same geometry written from a different starting vertex hashes
    differently.

    Canonicalisation covers the starting vertex and nothing else, so a
    ring and its reverse are two geometries and hash apart. That is the
    intended reading: direction is meaning, not spelling.

    Args:
        blob: The binary blob.

    Returns:
        A 32-byte digest.

    Example:
        >>> blob = encode_geometry([["NE(0625/0451)"]], "POINT", resolution=1)
        >>> len(geometry_hash(blob))
        32
    """
    _read_header(blob)
    return _keccak256(bytes(blob))


def geometry_to_tree(blob: bytes) -> bytes:
    """Derive the TreeBlob of a GeometryBlob's vertex set.

    One-way by construction. The TreeBlob is determined by the *set* of
    vertices alone, so distinct sequences, ring topologies, windings and
    edge models over the same vertices all yield the same TreeBlob.
    Coverage is preserved, identity is not: a ring and its reverse cover
    the same ground and arrive here as one tree.

    Args:
        blob: A GeometryBlob.

    Returns:
        The corresponding TreeBlob bytes.

    Raises:
        MalformedBlobError: If the blob fails structural validation.

    Example:
        >>> blob = encode_geometry([["NE(0625/0451)"]], "POINT", resolution=1)
        >>> geometry_to_tree(blob) == encode_tree("NE(0625/0451)")
        True
    """
    decoded = decode_geometry(blob)
    cells: "list[str]" = []
    seen: "set[str]" = set()
    for ring in decoded["rings"]:
        for index in ring:
            if index not in seen:
                seen.add(index)
                cells.append(index)
    groups: "dict[str, list[str]]" = {}
    for index in cells:
        groups.setdefault(index[:2], []).append(index[3:-1])
    return encode_tree(
        ",".join(
            f"{quadrant}({','.join(sorted(bodies))})"
            for quadrant, bodies in sorted(groups.items())
        )
    )
