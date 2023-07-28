# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'Kestrel'
copyright = 'University of Bristol'
author = 'Jake Langham & Mark Woodhouse'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'alabaster'
html_logo = 'logo.png'

# -- Options for EPUB output
epub_show_urls = 'footnote'

# -- Options for latex
latex_elements = {
        'preamble': r'''\usepackage{amssymb}
        \usepackage{bm}
        '''
}
