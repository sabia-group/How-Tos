# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('../../'))
html_extra_path = ['../../tools']
html_static_path = ['_static']


# -- Project information

project = 'SabIA-How-To'
copyright = '2024, SabIA Team'
author = 'SabIA Team'

release = '0.0'
version = '0.0.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'nbsphinx',
    'sphinx.ext.viewcode',  # This allows linking or embedding source code
    'sphinx_copybutton', # pip install sphinx-copybutton
    'sphinx_new_tab_link'
]

# testing web code rendering
autosummary_generate = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
