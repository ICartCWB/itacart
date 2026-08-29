# Changelog

Format based on [Keep a Changelog](https://keepachangelog.com/1.1.0/);
versioning follows [SemVer](https://semver.org/).

## [Unreleased]

### Fixed
- `build-system` required `setuptools>=61`, incompatible with the SPDX form
  `license = "MIT"` (PEP 639). Raised to `>=77`.
- `LICENSE` was missing, breaking `license-files` in `pyproject.toml`.

### Added
- Complete v1 public surface defined in stubs with type hints and docstrings.
- OGC DGGS Core / EAERS conformance suite mapping Frames 4 and 5 of the paper.
- `boundary` module covering the prime meridian, the antemeridian and the
  extension zones.
- `nominal_cell_area` and `effective_cell_area` kept distinct.
- `cell_to_anchor` (Req 12) and `cell_to_centroid` kept distinct.
- H3 v4 compatible aliases.

## [0.1.0a4] - 2026-08

### Changed
- `requires-python` raised from 3.8 to 3.10.
- Package layout defined under `src/itacart/`.

## [0.1.0a3]

- Placeholder release on PyPI.
