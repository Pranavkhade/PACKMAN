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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import packman

# -- Sphinx information  -----------------------------------------------------
master_doc = 'index'

# -- Project information -----------------------------------------------------

project = 'py-PACKMAN'
copyright = '2020, Pranav Khade(https://github.com/Pranavkhade)'
author = 'Pranav Khade(https://github.com/Pranavkhade)'

# The full version, including alpha/beta/rc tags
release = packman.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.napoleon','sphinx.ext.todo', 'sphinx.ext.autodoc', 'sphinx_rtd_theme']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['../_static']

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {"collapse_navigation" : False}

html_favicon = '../_static/gallary/logo.gif'

def run_apidoc(_):

    from sphinx.ext import apidoc

    ignore_paths = []

    argv = [
        "-f",
        "-T",
        "-e",
        "-M",
        "-o", "source",
        "../../packman"
    ] + ignore_paths

    apidoc.main(argv)

def setup(app):
    app.connect('builder-inited', run_apidoc)