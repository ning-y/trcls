# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import os, sys
import pytest

sys.path.append(os.path.join(os.path.abspath('../src')))
print(sys.path)

# -- Project information -----------------------------------------------------

project = 'trcls'
copyright = ''
author = 'LEE Ning Yuan'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinxarg.ext', 'sphinx.ext.autodoc', 'numpydoc']

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
html_theme = 'press'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# So, I want to customise the CSS, but since I'm using an external theme I don't
# think I can edit the HTML templates to include this. I'll have to use a raw
# HTML directive in my RST files to hack this in.
html_css_files = ['custom.css']

# Explicitly specify the sidebar contents to avoid the ugly quick search box
html_sidebars = {
    '**': [
        'util/sidetoc.html'
    ]
}

