# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

sys.path.insert(0, os.path.abspath("../../src"))

project = "itacart"
copyright = "2025, ITACaRT Project"
author = "ITACaRT Project"
release = "0.1.0a3"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",  # Imports docstrings from code
    "sphinx.ext.napoleon",  # Allows Sphinx to understand Google-style docstrings
    "sphinx.ext.viewcode",  # Adds links to source code in your docs
    "myst_parser",  # Allows contributing.rst to copy from CONTRIBUTING.md directly
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
