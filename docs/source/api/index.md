# API reference

Everything below is generated from the source at build time. The page for a
module is its module docstring plus one entry per public name, taken from the
docstring of that name — there is no second copy of the API to keep in sync.

```{eval-rst}
.. autosummary::
   :toctree: generated
   :template: autosummary/module.rst
   :recursive:

   itacart
```

## Reading order

The package is layered, and the layers correspond to delivery phases:

`itacart.constants`
: Immutable values transcribed from the paper — the WGS84 ellipsoid, Table 1,
  the refinement alphabets, the antemeridian extension zones. Data, not
  behaviour, with a single exception noted on the page.

`itacart.exceptions`
: Every error the package raises, rooted at `ITACaRTError`, so a caller can
  guard a whole pipeline with one `except` and never catch a stray
  `ValueError`.

`itacart.geodesy`
: The ellipsoidal sinusoidal projection and the geodesic solutions, computed
  from equations (1) and (2) of the paper. PROJ is not called at runtime.

`itacart.index`
: Parsing, composition and canonical form of the compositional index. The
  bridge between one string and the cells it addresses.

The remaining modules are declared in the package facade and delivered by later
phases; until then their pages show the contract without an implementation
behind it.
