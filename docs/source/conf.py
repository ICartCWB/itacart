"""Sphinx configuration for the itacart documentation."""

from __future__ import annotations

import itacart

project = "ITACaRT"
copyright = "2026, Israel Nunes da Silva, Gabriel Dietzsch, Elcio Hideiti Shiguemori"
author = "Israel Nunes da Silva"
release = itacart.__version__
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx_autodoc_typehints",
    "sphinx_copybutton",
    "myst_parser",
]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
napoleon_google_docstring = True
napoleon_numpy_docstring = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "shapely": ("https://shapely.readthedocs.io/en/stable", None),
    "geopandas": ("https://geopandas.org/en/stable", None),
}

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "furo"
html_static_path = ["../_static"]
html_title = f"ITACaRT {release}"

# O orientador recomenda maximizar graficos e imagens: figuras do artigo
# ficam em docs/_static e sao referenciadas nas paginas de conceito.
html_theme_options = {
    "source_repository": "https://github.com/ICartCWB/itacart/",
    "source_branch": "main",
    "source_directory": "docs/source/",
}
