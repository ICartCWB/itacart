# Contributing

## Environment

```bash
git clone https://github.com/ICartCWB/itacart
cd itacart
pip install -e ".[dev,docs,geo]"
```

## Before opening a PR

```bash
black src tests
isort src tests
flake8 src tests --max-line-length=88 --extend-ignore=E203,W503
mypy src
pytest
```

## Code discipline

- Docstrings, identifiers and comments in English.
- Type hints required on every public signature.
- Coverage target >= 85% per module.
- A fixed bug earns a permanent regression test.
- An architectural decision is recorded with its rationale **and the
  alternatives that were rejected**. Recording the rejection is what keeps a
  later change from relitigating a question that was already settled.

## Suggested order of implementation

The dependency graph between modules sets the order:

1. `constants`, `exceptions` — no dependencies
2. `geodesy` — projection and Vincenty
3. `index` — parser, compose, normalize
4. `resolutions` — table and scales
5. `cells` — quantization and inverse geometry
6. `boundary` — border families (depends on cells)
7. `hierarchy` — navigation
8. `topology` — neighbourhood (depends on boundary for deflection)
9. `geometry` — polyfill, densification, canonicalization
10. `serialization` — TreeBlob and GeometryBlob
11. `interop`, `engine`

## How the work is sequenced

Development proceeds in numbered phases, each one closed by a handoff that
records what was delivered, which decisions were taken, and what was left
open. The phase documents are coordination artefacts and are not versioned
here; the CHANGELOG is the record this repository keeps.
