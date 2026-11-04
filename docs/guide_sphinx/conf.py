"""Minimal Sphinx frontend for the handwritten chapters of the Flint site."""
from pathlib import Path
from sphinxcontrib.serializinghtml import JSONHTMLBuilder

project = 'MOM6'
root_doc = 'guide-index'
extensions = ['sphinx.ext.intersphinx', 'sphinxcontrib.bibtex']
bibtex_bibfiles = ['ocean.bib', 'references.bib', 'zotero.bib']
intersphinx_mapping = {
    'manual': ('https://mom6.readthedocs.io/en/main/',
               str(Path(__file__).resolve().parent / 'guide-inventory.inv')),
}
html_permalinks = False
nitpicky = True


class GuideJSONBuilder(JSONHTMLBuilder):
    """HTML fragments whose links match the site's sibling .html pages."""
    name = 'json'

    def get_target_uri(self, docname, typ=None):
        return docname + '.html'


def setup(app):
    app.add_builder(GuideJSONBuilder, override=True)
    return {'parallel_read_safe': True, 'parallel_write_safe': True}
