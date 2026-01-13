# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
import sys
import pathlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "py/"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PyORBIT3'
copyright = '2025, PyORBIT Collaboration'
author = 'PyORBIT Collaboration'
release = 'v3.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.coverage',
              'sphinx.ext.todo',
              'sphinx.ext.viewcode',
              'sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.napoleon',
              'sphinx_copybutton',
              'myst_parser',
              'breathe',
              'exhale']

napoleon_numpy_docstring = True
autosummary_imported_members = True

templates_path = ['_templates']
exclude_patterns = []

# -- Breathe and Exhale Options ---------------------------------------------------------
breathe_projects = {
        "PyORBIT3": "./_doxygen/xml",
    }
breathe_default_project = "PyORBIT3"

doxyfile = '\n'.join([
    "INPUT = ../../src",
    "EXCLUDE_PATTERNS = *wrap* *_init.cc",
    "PREDEFINED += DOXYGEN_SHOULD_SKIP_THIS",
    "PREDEFINED += PyObject_HEAD=\"PyObject ob_base;\""
    ])

exhale_args = {
        "containmentFolder": "./api",
        "rootFileName": "pyorbit_root.rst",
        "doxygenStripFromPath": "..",
        "rootFileTitle": "PyORBIT3 C++ API",
        "createTreeView": True,
        "exhaleExecutesDoxygen": True,
        "exhaleDoxygenStdin": doxyfile,
        }

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
html_show_source_link = True
html_theme_options = {
    'github_url': 'https://github.com/PyORBIT-Collaboration/PyORBIT3',
    'logo': {
        'text': "PyORBIT3",
    }
}
