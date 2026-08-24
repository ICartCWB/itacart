# Changelog

Formato baseado em [Keep a Changelog](https://keepachangelog.com/1.1.0/);
versionamento conforme [SemVer](https://semver.org/lang/pt-BR/).

## [Unreleased]

### Fixed
- `build-system` exigia `setuptools>=61`, incompativel com a forma SPDX
  `license = "MIT"` (PEP 639). Elevado para `>=77`.
- `LICENSE` ausente, quebrando `license-files` do `pyproject.toml`.

### Added
- Superficie publica completa da v1 definida em stubs com type hints e docstrings.
- Suite de conformidade OGC DGGS Core / EAERS mapeando os Frames 4 e 5 do artigo.
- Modulo `boundary` cobrindo meridiano principal, antemeridiano e zonas de extensao.
- Distincao entre `nominal_cell_area` e `effective_cell_area`.
- Distincao entre `cell_to_anchor` (Req 12) e `cell_to_centroid`.
- Aliases compativeis com H3 v4.

## [0.1.0a4] - 2026-08

### Changed
- `requires-python` elevado de 3.8 para 3.10.
- Estrutura de pacote definida sob `src/itacart/`.

## [0.1.0a3]

- Placeholder no PyPI.
