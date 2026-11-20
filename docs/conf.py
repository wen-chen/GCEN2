# docs/conf.py
"""Sphinx configuration for GCEN2 documentation."""
import os
import sys

# Make the gcen2 package importable for autodoc
sys.path.insert(0, os.path.abspath('..'))

project = 'GCEN2'
copyright = '2026, Wen Harold Chen'
author = 'Wen Harold Chen'
release = '2.0.0'
version = '2.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',       # NumPy / Google style docstrings
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'pandas': ('https://pandas.pydata.org/docs', None),
    'scipy': ('https://docs.scipy.org/doc/scipy', None),
}

autosummary_generate = True
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
    'member-order': 'bysource',
}
autodoc_mock_imports = [
    'anndata', 'scanpy', 'matplotlib', 'seaborn',
    'pyvis', 'networkx', 'Jinja2', 'jinja2',
    'weasyprint', 'typer', 'magic',
]

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'navigation_depth': 4,
    'collapse_navigation': False,
}
html_static_path = []

templates_path = []
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'plans', 'design']

nitpicky = False
