# Configuration file for the Sphinx documentation builder.

# -- Project information

import os
import sys
sys.path.insert(0, os.path.abspath('../../synthpop/'))
sys.path.insert(0, os.path.abspath('../'))
sys.path.insert(0, os.path.abspath('../../'))
sys.path.insert(0, os.path.abspath('../../synthpop/modules/'))
sys.path.insert(0, os.path.abspath('../../synthpop/modules/age/'))
sys.path.insert(0, os.path.abspath('../../synthpop/modules/metallicity/'))

#import migrate_interactive_part
#migrate_interactive_part.migrate('../../synthpop')

project = 'SynthPop'
copyright = 'tbd'
author = 'M Huston'
master_doc = 'index'

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

extensions.append('autoapi.extension')
extensions.append('sphinx.ext.autosectionlabel')
autoapi_dirs = ['../../synthpop','../../synthpop/modules/']
autoapi_ignore = ['conf.py']

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'
html_css_files = ["custom.css"]

# -- Options for EPUB output
epub_show_urls = 'footnote'
