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
# import sphinx_rtd_theme
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'Using Python in an HPC environment'
copyright = '2025, UPPMAX/HPC2N/LUNARC/InfraVis'
author = 'UPPMAX/HPC2N/LUNARC/InfraVis'
github_user = "UPPMAX"
github_repo_name = "HPC-python"  # auto-detected from dirname if blank
github_version = "main"
conf_py_path = "/docs/"
# The full version, including alpha/beta/rc tags
release = '2.0'

# The master toctree document.
master_doc = 'manual/index'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ["sphinx_lesson",
    "sphinx.ext.githubpages",
    "sphinx_rtd_theme_ext_color_contrast",
    'sphinx_copybutton'

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

jupyter_execute_notebooks = "cache"

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
# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# exclude_patterns = [
#     "README*",
#     "_build",
#     "Thumbs.db",
#     ".DS_Store",
#     "jupyter_execute",
#     "*venv*",
# ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_logo = "img/hpc2n-lunarc-uppmax-hpc-course.png"
#html_logo = "img/logo-python-hpc-course.svg"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
def setup(app):
    app.add_css_file('custom_theme.css')

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

