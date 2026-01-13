# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ACT'
copyright = '2026, David van der Spoel, Paul J. van Maaren and Mohammad M. Ghahremanpour'
author = 'David van der Spoel, Paul J. van Maaren and Mohammad M. Ghahremanpour'
release = '1.2 :math:`\\beta`'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']

extensions = ['sphinxcontrib.bibtex']
bibtex_bibfiles = ['library.bib']
bibtex_default_style = 'unsrt'
