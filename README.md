Template for the Read the Docs tutorial
=======================================

This GitHub template includes fictional Python library
with some basic Sphinx docs.

The website is served at

https://how-tos.readthedocs.io/en/latest/


--

To serve locally run:

`pip install -r docs/requirements.txt`

followed by

`sphinx-autobuild docs/source docs/_build/html`

Then just open `http://127.0.0.1:8000` on your browser.


## Guide to Contributing

We suggest the following workflow for adding tutorials:

  1. In `docs/source` create a sub-directory dedicated to your tutorial, e.g. `cheese`

  1. Within `docs/source/cheese`, create an RST (reStructuredText) file where you will typeset the howto. 

  1. Register the new entry in the table of contents by adding the line `cheese/cheese` to `docs/source/index.rst` under the `.. toctree::` block.
