# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# OS path commands should be uncommented in conf.py file
# so that html utility could access the right project files to
# generate documentation.
import os
import sys

sys.path.insert(0, os.path.abspath("../src/"))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "X-FP"
copyright = "2022-2023, Marcel Baltruschat and Aishvarya Tandon"
author = "Marcel Baltruschat and Aishvarya Tandon"
release = "1.1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- MyST configuration ------------------------------------------------------
# https://myst-parser.readthedocs.io/en/latest/configuration.html

myst_all_links_external = True


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "nature"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_logo = "../src/xfp/resources/X-FP_logo.png"
