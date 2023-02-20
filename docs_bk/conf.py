# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("."))

# rebuilding with sphinx-autobuild:
# sphinx-autobuild docs docs/_build/html --watch pyhdx

# -- Project information -----------------------------------------------------

project = "PyHDX"
copyright = "2022, Jochem Smit"
author = "Jochem Smit"


# -- General configuration ---------------------------------------------------

# Custom param gui formatting
def setup(app):
    from paramdoc import param_format_basic

    app.connect("autodoc-process-docstring", param_format_basic, priority=-100)


# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",  # autodoc docstrings
    "sphinxcontrib.bibtex",  # allow bibtex references
    "sphinx.ext.intersphinx",
    "nbsphinx",  # jupyter notebooks in docs
    "sphinx.ext.mathjax",  # render latex style math
    "sphinx.ext.viewcode",  # view code links
    "sphinx_autodoc_typehints",  # autodoc type hints
    "sphinx.ext.napoleon",  # google style docstrings
]

# bibtex configuration
bibtex_bibfiles = ["refs.bib"]
bibtex_reference_style = "label"

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("http://docs.scipy.org/doc/scipy/reference/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "torch": ("https://pytorch.org/docs/stable/", None),
    "symfit": ("https://symfit.readthedocs.io/en/stable/", None),
}
