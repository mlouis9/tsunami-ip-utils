import os, sys
import logging
from docutils.parsers.rst import Directive
import inspect
from sphinx_gallery.scrapers import figure_rst

src_path = os.path.abspath('../../src')
ext_path = os.path.abspath('./_ext')
sys.path.insert(0, src_path)
sys.path.insert(0, ext_path)

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

# -- Plotly Image Scraper -----------------------------------------------------
def plotly_scraper(block, block_vars, gallery_conf):
    from IPython.display import HTML

    # Check if 'fig' is in the example_globals and if it can generate HTML
    fig = block_vars['example_globals'].get('fig', None)
    if fig and hasattr(fig, 'to_html'):
        html_output = fig.to_html(include_plotlyjs='cdn', full_html=False)
        # Wrap the HTML output in necessary RST for Sphinx to handle it correctly
        html_rst = f"""
.. raw:: html

    {html_output}

        """
        return html_rst
    else:
        # Handle cases where 'fig' is not available or does not have 'to_html'
        print("No 'fig' found or 'fig' lacks 'to_html' method")
        return ''

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_rtd_theme',
    'sphinx_autodoc_typehints',
    'sphinx.ext.viewcode',         # Add links to source code
    'sphinx.ext.mathjax',          # Render math equations using MathJax
    'globalparam',                 # Custom extension for global parameters
    'sphinx_gallery.gen_gallery',  # Generate gallery of examples
]

sphinx_gallery_conf = {
    'examples_dirs': '../../examples',  # Path to example scripts
    'gallery_dirs': 'auto_examples',  # Path to save gallery generated output
    'filename_pattern': r'plot_.*\.py$',  # Adjusted regex to ensure it captures all intended files
    'example_extensions': ['.py'],
    'image_scrapers': ('matplotlib', plotly_scraper),  # If using matplotlib for plots
    'doc_module': ('tsunami_ip_utils',)
}

autodoc_default_options = {
    'member-order': 'bysource',
    'special-members': '__init__',
    'inherited-members': True,
}

templates_path = ['_templates']
exclude_patterns = []
show_headings = False

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Custom Class for Filtering Duplicate Object Warnings --------------------
"""NOTE: These warnings could not be suppressed using the sphinx suppress_warnings option,
so this manual approach was necessary. NOTE: also that although the .. automodule:: my_module directive is
used twice for each module (once for private API and once for public API), there is no overlap in the index because
the methods/classes in the private and public APIs are mutually exclusive."""

class FilterDuplicateObjectWarnings(logging.Filter):
    def filter(self, record):
        is_duplicate_warning = 'duplicate object description of %s, other instance in %s, use :no-index: for one of them' in record.msg
        return not is_duplicate_warning

# -- Custom Configuration For Generating Public and Private API Docs --------

# Suppress warnings from documenting the same module twice (public and private)

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

    # Add filter to Sphinx logger to suppress duplicate object warnings
    logger = logging.getLogger('sphinx')
    logger.addFilter(FilterDuplicateObjectWarnings())

def has_private_members(obj):
    """Check if the given object has any private members."""
    for name, member in inspect.getmembers(obj):
        if name.startswith('_') and not name.startswith('__'):
            return True
    return False

def skip_by_api_type(api_type, name, obj):
    if inspect.isclass(obj):
        # Decide based on whether the class has private members
        if api_type == 'public':
            return False  # Do not skip, include in documentation for public API
        elif api_type == 'private' and not has_private_members(obj):
            return True  # Skip if it's meant to be private but has no private members
    else:
        # For non-class objects, use the existing name-based logic
        private_method = name.startswith('_') and ( not name.startswith('__') )
        if api_type == 'public' and private_method:
            return True
        elif api_type == 'private' and not private_method:
            return True

def skip_member(app, what, name, obj, skip, options):
    api_type = app.config.api_type
    skip = skip_by_api_type(api_type, name, obj)

    # Exclude inherited members for specific classes
    classes_to_exclude_inherited_members = ['EnhancedPlotlyFigure']
    if isinstance(obj, type) and obj.__name__ in classes_to_exclude_inherited_members:
        options['inherited-members'] = {}

    return skip
