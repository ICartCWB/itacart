"""GeometryBlob: binary encoding of an ordered ITACaRT vertex sequence.

Preserves what TreeBlob discards: vertex order, ring topology, edge model,
densification parameters and OGC SFA type. That makes it the form to hash
when a geometry must be identified rather than merely covered.

Header layout, 9 bytes:

===========  ====  ==============================  ==========================
Byte(s)      Bits  Field                           Values
===========  ====  ==============================  ==========================
``0x00``        8  Magic                           ``0xC8``
``0x01``        8  Format version                  ``0x01``
``0x02`` hi     4  Canonical profile               ``0x0`` MIN_LEX_CYCLIC_ROTATION
``0x02`` lo     4  Edge model                      ``0x1`` WGS84_GEODESIC
``0x03`` hi     4  Densification model             ``0x1`` ORTHODROMIC_VINCENTY,
                                                   ``0x2`` ASSUMED_BY_PRODUCER
``0x03`` lo     4  Resolution mode                 ``0x0`` UNIFORM
``0x04``-``05``16  ``max_segment_m``               ``0x03E8`` = 1000 m
``0x06`` hi     4  Uniform resolution              ``0xD`` = res 13
``0x06`` lo     4  Geometry type (OGC SFA)         ``0x0`` POINT .. ``0x5`` MULTIPOLYGON
``0x07``-``08``16  ``n_components``                1 for single, N for MULTI*
===========  ====  ==============================  ==========================

Components, parts and rings follow, then the vertex stream. Vertices are
delta-encoded by ``shared_levels``, the count of hierarchy levels a vertex
shares with its predecessor; vertex 0 is always written in full.

Measured density on real cadastral datasets ranges from 27 to 38
bits/vertex, rising with geographic extent as shared prefixes shorten.

Origem: itacart_core/geometry_blob.py (F2.2 v2),
spec em docs/geometry_encoding_spec.md.
"""

from __future__ import annotations

from typing import Any

__all__ = [
    "encode_geometry",
    "decode_geometry",
    "validate_geometry",
    "geometry_hash",
    "geometry_to_tree",
    "read_geometry_type",
]


def encode_geometry(
    rings: list[list[str]],
    geometry_type: str = "POLYGON",
    edge_model: str = "WGS84_GEODESIC",
    densification_model: str = "ORTHODROMIC_VINCENTY",
    max_segment_m: int = 1000,
    resolution: int = 13,
    canonicalize: bool = True,
) -> bytes:
    """Encode ordered vertex rings as a GeometryBlob.

    Args:
        rings: Rings as lists of atomic index strings, exterior first.
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
    """
    raise NotImplementedError


def decode_geometry(blob: bytes) -> dict[str, Any]:
    """Decode a GeometryBlob into rings and profile metadata.

    Args:
        blob: The binary blob.

    Returns:
        A mapping with ``rings``, ``geometry_type``, ``edge_model``,
        ``densification_model``, ``max_segment_m`` and ``resolution``.

    Raises:
        MalformedBlobError: If the blob fails structural validation.
    """
    raise NotImplementedError


def validate_geometry(blob: bytes) -> None:
    """Validate blob structure, raising on the first problem found.

    Args:
        blob: The binary blob.

    Raises:
        MalformedBlobError: With a message identifying the failure.
    """
    raise NotImplementedError


def geometry_hash(blob: bytes) -> bytes:
    """Content hash of a canonical GeometryBlob.

    Meaningful only for canonicalised blobs: without canonicalisation the
    same geometry written from a different starting vertex hashes
    differently.

    Args:
        blob: The binary blob.

    Returns:
        A 32-byte digest.
    """
    raise NotImplementedError


def geometry_to_tree(blob: bytes) -> bytes:
    """Derive the TreeBlob of a GeometryBlob's vertex set.

    One-way by construction. The TreeBlob is determined by the *set* of
    vertices alone, so distinct sequences, ring topologies and edge
    models over the same vertices all yield the same TreeBlob. Coverage
    is preserved, identity is not.

    Args:
        blob: A GeometryBlob.

    Returns:
        The corresponding TreeBlob bytes.
    """
    raise NotImplementedError


def read_geometry_type(blob: bytes) -> str:
    """Read the OGC SFA type from the header without full decoding.

    Args:
        blob: The binary blob.

    Returns:
        The geometry type name.

    Raises:
        MalformedBlobError: If the header is unreadable.
    """
    raise NotImplementedError
