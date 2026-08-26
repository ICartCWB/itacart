"""Sphinx configuration for the itacart documentation.

The build reads the **working tree**, never an installed copy — see the
``sys.path`` bootstrap below. On a development branch the version pip has
installed is almost never the version under test, and documentation generated
from the wrong source is worse than no documentation at all.

The API reference is not written by hand. ``sphinx.ext.autosummary`` walks the
package recursively and generates one page per module and per public name out
of the docstrings themselves, so a module delivered by a later phase appears in
the documentation without editing anything here.
"""

from __future__ import annotations

import sys
from pathlib import Path

# --------------------------------------------------------------------------
# Import the package from the working tree
# --------------------------------------------------------------------------

SOURCE_DIR = Path(__file__).resolve().parent
ROOT = SOURCE_DIR.parents[1]

sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(SOURCE_DIR / "_ext"))

import itacart  # noqa: E402  (must follow the sys.path bootstrap)

_loaded_from = Path(itacart.__file__).resolve().parent
if _loaded_from != (ROOT / "src" / "itacart").resolve():
    raise RuntimeError(
        f"itacart was imported from {_loaded_from}, not from the working tree "
        f"at {ROOT / 'src' / 'itacart'}. An installed copy is shadowing it; "
        "uninstall it or build in a clean environment."
    )

# --------------------------------------------------------------------------
# Project metadata, read from the package rather than restated
# --------------------------------------------------------------------------

project = "ITACaRT"
copyright = "2026, Israel Nunes da Silva, Gabriel Dietzsch, Elcio Hideiti Shiguemori"
author = "Israel Nunes da Silva"
release = itacart.__version__
version = release

# --------------------------------------------------------------------------
# Extensions
# --------------------------------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx.ext.coverage",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "myst_parser",
    "phase_figures",  # local, in _ext/
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

# --------------------------------------------------------------------------
# autosummary and autodoc — what makes the reference self-generating
# --------------------------------------------------------------------------

autosummary_generate = True

# The package facade re-exports 109 names from the submodules. Documenting them
# at package level as well would duplicate every entry and produce a page of
# cross-references to itself, so the facade shows only what it defines and the
# submodule pages carry the rest.
autosummary_imported_members = False

# Deliberately no "members" here. A module page is an index: its own docstring
# plus the summary tables the template builds. The detail lives on one page per
# name, generated into api/generated/. Turning members on as well documents
# every object twice — once inline and once on its own page — and Sphinx says
# so, loudly, as "duplicate object description".
autodoc_default_options = {
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
autodoc_preserve_defaults = True
# Only the optional extras are mocked. shapely is a hard runtime dependency,
# so mocking it would render every shapely type in a signature as a mock.
autodoc_mock_imports = ["geopandas", "joblib"]

# Google style throughout, as the repository already had it. The one file that
# used NumPy style, itacart_core/compositional_index.py, is an origin document
# and not part of the package. Flip numpy_docstring back on if a module ported
# from itacart_core ever lands with "Parameters / ----------" headings.
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_use_admonition_for_notes = True
napoleon_use_rtype = False

always_document_param_types = True
typehints_defaults = "comma"

coverage_show_missing_items = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "shapely": ("https://shapely.readthedocs.io/en/stable", None),
    "geopandas": ("https://geopandas.org/en/stable", None),
}

# --------------------------------------------------------------------------
# MyST
# --------------------------------------------------------------------------

myst_enable_extensions = ["colon_fence", "deflist", "dollarmath"]
myst_heading_anchors = 3

# --------------------------------------------------------------------------
# HTML output
# --------------------------------------------------------------------------

html_theme = "furo"
html_title = f"itacart {release}"
# O orientador recomenda maximizar graficos e imagens: figuras do artigo
# ficam em docs/_static e sao referenciadas nas paginas de conceito.
#
# Two entries, on purpose. "../_static" is the repository's own directory and
# is kept so that any page already pointing at _static/f1/... keeps resolving.
# "_static" is source-local, for assets that belong to the documentation
# rather than to a phase. Drop the first if no page references it directly:
# the phase galleries do not need it, they copy what they use.
html_static_path = ["_static", "../_static"]
html_copy_source = False
html_show_sphinx = False

html_theme_options = {
    "source_repository": "https://github.com/ICartCWB/itacart/",
    "source_branch": "main",
    "source_directory": "docs/source/",
    "light_css_variables": {"color-brand-primary": "#2f6f9f"},
    "dark_css_variables": {"color-brand-primary": "#7fb3d5"},
}

# --------------------------------------------------------------------------
# Phase figure galleries
# --------------------------------------------------------------------------
#
# The verification notebooks write their figures to docs/_static/fN/, which is
# outside this source directory: the notebooks themselves are not versioned
# (P-1.10) but the figures they produce are. The phase_figures extension copies
# each directory into _generated/ at build time and writes one gallery page per
# phase, so a new phase's figures need no edit here either.

phase_figures_root = ROOT / "docs" / "_static"
phase_figures_titles = {
    "f0": "F0 — Bootstrap",
    "f1": "F1 — Geodesy",
    "f2": "F2 — Compositional index",
    "f3": "F3 — Resolutions and cells",
    "f4": "F4 — Boundary behaviour",
    "f5": "F5 — Hierarchy",
    "f6": "F6 — Topology",
    "f7": "F7 — Vector geometry",
    "f8": "F8 — Serialization",
    "f9": "F9 — Interoperability and conformance",
}
