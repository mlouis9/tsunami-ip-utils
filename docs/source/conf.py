import os, sys
sys.path.insert(0, os.path.abspath('../../src'))

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'tsunami_ip_utils'
copyright = '2024, Matthew Louis'
author = 'Matthew Louis'
release = '0.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_rtd_theme'
]

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']


# -- Custom Configuration For Generating Public and Private API Docs --------
from docutils.parsers.rst import Directive

# Suppress warnings from documenting the same module twice (public and private)
suppress_warnings = [
    'automodule.duplicate',
]

class ApiTypeDirective(Directive):
    has_content = False
    required_arguments = 1

    def run(self):
        env = self.state.document.settings.env
        env.config.api_type = self.arguments[0]
        return []

def setup(app):
    app.add_config_value('api_type', 'public', 'env')
    app.add_directive('api_type', ApiTypeDirective)
    app.connect('autodoc-skip-member', skip_member)

def skip_member(app, what, name, obj, skip, options):
    api_type = app.config.api_type
    if api_type == 'public' and name.startswith('_'):
        return True
    elif api_type == 'private' and not name.startswith('_'):
        return True
    return skip