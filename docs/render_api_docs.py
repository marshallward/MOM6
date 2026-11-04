"""MOM6's HTML presentation layer for independently generated reST pages."""
import argparse
import contextlib
from html import escape
from html.parser import HTMLParser
from pathlib import Path
import shutil
from string import Template
import xml.etree.ElementTree as ET
from urllib.parse import unquote

from docutils import nodes
from docutils.core import publish_parts
from docutils.readers.standalone import Reader
from docutils.transforms import Transform
from render_guide_docs import GUIDE_PAGES, SCIENCE_PAGES, render_guides


THEME = Path(__file__).parent / 'api_theme'


def graph_text(image):
    """Build linked text edges from the same SVG used by the visual graph."""
    tree = ET.parse(unquote(image['uri']))
    ns = {'svg': 'http://www.w3.org/2000/svg'}
    labels = {}
    edges = []
    for group in tree.findall('.//svg:g', ns):
        title = group.find('svg:title', ns)
        if title is None:
            continue
        if group.get('class') == 'node':
            link = group.find('.//svg:a', ns)
            url = (link.get('{http://www.w3.org/1999/xlink}href', link.get('href'))
                   if link is not None else None)
            # Qualified identities disambiguate same-named procedures/modules.
            labels[title.text] = url
        elif group.get('class') == 'edge':
            edges.append(title.text.split('->', 1))
    listing = nodes.bullet_list(classes=['call-graph-links'])
    for caller, callee in edges:
        paragraph = nodes.paragraph()
        for index, name in enumerate((caller, callee)):
            if index:
                paragraph += nodes.Text(' → ')
            url = labels.get(name)
            paragraph += (nodes.reference('', name, refuri=url) if url else nodes.Text(name))
        listing += nodes.list_item('', paragraph)
    return listing


class GraphPresentation(Transform):
    """Replace graph-backed call lists with a diagram and linked text view."""
    default_priority = 850

    def apply(self):
        titles = {'Callees': 'Calls into', 'Callers': 'Called by'}
        for container in list(self.document.findall(nodes.container)):
            if 'api-procedure' not in container['classes']:
                continue
            # Only direct children belong to this procedure, not nested ones.
            for rubric in list(container.children):
                if not isinstance(rubric, nodes.rubric) or rubric.astext() not in titles:
                    continue
                position = container.index(rubric)
                if position + 1 >= len(container):
                    continue
                image = container[position + 1]
                if not isinstance(image, nodes.image) or not image['uri'].endswith('.svg'):
                    continue
                text = graph_text(image)
                title = rubric.astext()
                for listing in list(container.children):
                    if not isinstance(listing, nodes.definition_list):
                        continue
                    for item in list(listing.children):
                        term = next((child for child in item if isinstance(child, nodes.term)), None)
                        if term is not None and term.astext() == titles[title]:
                            listing.remove(item)
                    if not listing.children:
                        container.remove(listing)
                position = container.index(rubric)
                opening = nodes.raw('',
                    '<details class="call-graph-panel"><summary>' + title + ' call graph</summary>'
                    '<div class="call-graph" tabindex="0" role="region" '
                    'aria-label="Call graph (scroll to explore)">', format='html')
                text_open = nodes.raw('', '</div><details class="call-graph-text">'
                                      '<summary>Text view</summary>', format='html')
                closing = nodes.raw('', '</details></details>', format='html')
                container[position:position + 2] = [opening, image, text_open, text, closing]


class GraphReader(Reader):
    def get_transforms(self):
        return super().get_transforms() + [GraphPresentation]


class Outline(HTMLParser):
    """Read navigation labels from rendered headings and definition targets."""
    def __init__(self):
        super().__init__()
        self.sections = []
        self.definitions = []
        self.items = []
        self.capture = None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'section':
            self.sections.append(attrs.get('id'))
        elif tag == 'dl':
            self.definitions.append(attrs.get('id'))
        elif tag == 'h2' and self.sections and self.sections[-1]:
            self.capture = (tag, self.sections[-1], 'section', [])
        elif tag == 'dt' and self.definitions and self.definitions[-1]:
            self.capture = (tag, self.definitions[-1], 'definition', [])
            self.definitions[-1] = None

    def handle_data(self, text):
        if self.capture is not None:
            self.capture[3].append(text)

    def handle_endtag(self, tag):
        if self.capture is not None and tag == self.capture[0]:
            _, target, kind, parts = self.capture
            self.items.append((target, kind, ' '.join(''.join(parts).split())))
            self.capture = None
        if tag == 'section' and self.sections:
            self.sections.pop()
        elif tag == 'dl' and self.definitions:
            self.definitions.pop()


def render_reference(rst_dir, html_dir, warning_stream):
    """Render standard reST with local theme assets; source documents are untouched."""
    html_dir.mkdir(parents=True, exist_ok=True)
    assets = html_dir / 'assets'
    assets.mkdir(exist_ok=True)
    for name in ('api.css', 'api.js', 'mathjax-3.2.2.js', 'mathjax-LICENSE'):
        shutil.copyfile(THEME / name, assets / name)
    template = Template((THEME / 'page.html').read_text(encoding='utf-8'))
    files = sorted(rst_dir.glob('*.rst'), key=lambda path: path.stem.lower())
    guide_bodies = render_guides(rst_dir, warning_stream, html_dir)
    modules = [path.stem for path in files
               if path.stem not in ('index', 'api-index') and path.stem not in guide_bodies]
    labels = {}
    for path in files:
        if path.stem in guide_bodies:
            labels[path.stem] = {**GUIDE_PAGES, **SCIENCE_PAGES}[path.stem]
            continue
        # Flint emits an overlined document title. Read only that header rather
        # than retaining all rendered pages just to build the sidebar.
        with path.open(encoding='utf-8') as stream:
            next(stream, '')
            title = next(stream, '').strip()
        suffix = ' module reference'
        labels[path.stem] = title[:-len(suffix)] if title.endswith(suffix) else path.stem
    for path in files:
        # publish_parts resolves embedded image paths against CWD. Inline SVG
        # preserves clickable graph nodes without separate deployed assets.
        source = path.read_text(encoding='utf-8')
        source_path = str(path.resolve())
        if path.stem in guide_bodies:
            parts = {'html_body': guide_bodies[path.stem]}
        else:
            with contextlib.chdir(path.parent):
                parts = publish_parts(source=source, source_path=source_path, reader=GraphReader(),
                                      writer='html5', settings_overrides={
                                          'warning_stream': warning_stream, 'halt_level': 2,
                                          'stylesheet_path': [], 'embed_stylesheet': False,
                                          'image_loading': 'embed',
                                      })
        outline = Outline()
        outline.feed(parts['html_body'])
        navigation = []
        for name in modules:
            current = ' aria-current="page"' if name == path.stem else ''
            navigation.append(f'<li><a href="{escape(name)}.html"{current}>{escape(labels[name])}</a></li>')
        contents = ''.join(
            f'<li class="{kind}"><a href="#{escape(target)}" title="{escape(label)}">'
            f'{escape(label)}</a></li>' for target, kind, label in outline.items)
        home = path.stem == 'index'
        module_index = path.stem == 'api-index'
        guide = path.stem in guide_bodies
        guides = []
        for name, label in GUIDE_PAGES.items():
            if name in labels:
                current = ' aria-current="page"' if name == path.stem else ''
                guides.append(f'<a href="{name}.html"{current}>{escape(label)}</a>')
        title = ('MOM6 documentation' if home else
                 'MOM6 API reference' if module_index else labels[path.stem])
        page = template.substitute(
            title=escape(title), body=parts['html_body'], modules='\n'.join(navigation),
            guides='\n'.join(guides),
            contents=contents,
            math_script=('<script defer src="assets/mathjax-3.2.2.js"></script>'
                         if 'class="math ' in parts['html_body'] else ''),
            page_class=('home-page' if home else 'module-index' if module_index else
                        'guide-page' if guide else 'module-page'),
            breadcrumb='Home' if home else 'API reference' if module_index else escape(labels[path.stem]),
            module_count=str(len(modules)),
            contents_hidden=' hidden' if not contents else '',
        )
        (html_dir / (path.stem + '.html')).write_text(page, encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rst_directory', type=Path)
    parser.add_argument('html_directory', type=Path)
    args = parser.parse_args()
    import sys
    render_reference(args.rst_directory, args.html_directory, sys.stderr)


if __name__ == '__main__':
    main()
