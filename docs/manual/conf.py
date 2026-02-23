# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import datetime
import os
import sys

project = u'ACT'

#execfile('conf-vars.py')

# The master toctree document.
master_doc = 'index'

copyright = str(datetime.datetime.now().year) + u', ACT development team'
thisyear_string = str(datetime.datetime.now().year)
github_user = "AlexandriaChemistry"
github_repo_name = "ACT"
github_version = "main"
conf_py_path = "/docs/"
# The full version, including alpha/beta/rc tags
release = '1.2b'

extlinks = {'issue': ('https://github.com/AlexandriaChemistry/ACT/issues/%s',
                      'Issue ')}
author = 'David van der Spoel, Paul J. van Maaren and Mohammad M. Ghahremanpour'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ['_templates']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
#html_static_path = ['_static']

extensions = ['sphinxcontrib.bibtex']
bibtex_bibfiles = ['manual-refs.bib']
bibtex_default_style = 'unsrt'

# The master toctree document.
master_doc = 'index'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["sphinx_lesson",
    "sphinx.ext.githubpages",
    "sphinx_rtd_theme_ext_color_contrast",
    'sphinx_copybutton',
    'sphinxcontrib.bibtex'
]

mermaid_output_format = 'raq'
mermaid_output_format = "png"
mermaid_params = [
    "--theme",
    "forest",
    "--backgroundColor",
    "transparent",
    '-p' 'docs/puppeteer-config.json'
]

#jupyter_execute_notebooks = "cache"

myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist",
    "colon_fence",
]
copybutton_exclude = '.linenos, .gp, .go'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'classic'

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
#html_logo = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']
#def setup(app):
#    app.add_css_file('custom_theme.css')

# HTML context:
from os.path import basename, dirname, realpath

html_context = {
    "display_github": True,
    "github_user": github_user,
    # Auto-detect directory name.  This can break, but
    # useful as a default.
    # "github_repo": github_repo_name or basename(dirname(realpath(__file__))),
    "github_repo": github_repo_name or basename(dirname(realpath(__file__))),
    "github_version": github_version,
    "conf_py_path": conf_py_path,
}
