# itacart

An equal-area parallelogram Discrete Global Grid System for terrestrial
cadastral mapping, tessellated directly on the WGS84 ellipsoid.

This is the reference implementation of **ITACaRT** — *ITA Cadastral
Ellipsoidal Reference Tessellation* — described in Silva, Dietzsch &
Shiguemori (2025), *Revista Brasileira de Cartografia*, v. 77,
[10.14393/rbcv77n0a-79281](https://doi.org/10.14393/rbcv77n0a-79281).
When this documentation and the paper disagree, the paper is right and the
disagreement is a bug.

```{code-block} python
import itacart

cell = itacart.geo_to_cell(-46.6328862, -23.5508962, resolution=13)
lon, lat = itacart.cell_to_centroid(cell)
ring = itacart.cell_to_boundary(cell, close=True)
```

## What is implemented

The package is being delivered phase by phase. The API reference below is
generated from the source, so a function that appears there exists; a function
listed in the paper but absent from the reference has not been written yet.

```{toctree}
:maxdepth: 2
:caption: Contents

api/index
_generated/figures/index
building
```

## The compositional index

Every public operation accepts a **compositional index**, a single string that
may address one cell or a whole region:

```text
SE(1400/0374(3(C2(3))))                one cell, resolution 4
NW(0625/0451(1(E1(3(B2(4(A2,B2)))))))  two sibling cells under one parent
```

Parentheses descend one resolution level; commas separate siblings at the same
level. Operations that return one value per cell return a list aligned with
{func}`itacart.index.decompose` order.

Uniqueness — OGC DGGS Core requirement 13 — holds against the canonical form
produced by {func}`itacart.index.normalize`, not against arbitrary spellings:
`4(1,2,3,4)` and `4` denote the same space.
