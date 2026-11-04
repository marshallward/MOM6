"""MOM6's HTML presentation layer for independently generated reST pages."""
import argparse
import contextlib
import re
from html import escape
from html.parser import HTMLParser
from pathlib import Path
import shutil
from string import Template

from docutils.core import publish_parts


THEME = Path(__file__).parent / 'api_theme'


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
    for name in ('api.css', 'api.js'):
        shutil.copyfile(THEME / name, assets / name)
    template = Template((THEME / 'page.html').read_text(encoding='utf-8'))
    files = sorted(rst_dir.glob('*.rst'), key=lambda path: path.stem.lower())
    modules = [path.stem for path in files if path.stem != 'index']
    labels = {}
    for path in files:
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
        with contextlib.chdir(path.parent):
            parts = publish_parts(source=source, source_path=source_path,
                                  writer='html5', settings_overrides={
                                      'warning_stream': warning_stream, 'halt_level': 2,
                                      'stylesheet_path': [], 'embed_stylesheet': False,
                                      'image_loading': 'embed',
                                  })
        # Preserve natural graph label sizes on narrow screens. The wrapper
        # provides keyboard-accessible scrolling instead of shrinking the SVG.
        parts['html_body'] = re.sub(
            r'(<svg\b.*?</svg>)',
            r'<div class="call-graph" tabindex="0" role="region" '
            r'aria-label="Call graph (scroll to explore)">\1</div>',
            parts['html_body'], flags=re.DOTALL)
        outline = Outline()
        outline.feed(parts['html_body'])
        navigation = []
        for name in modules:
            current = ' aria-current="page"' if name == path.stem else ''
            navigation.append(f'<li><a href="{escape(name)}.html"{current}>{escape(labels[name])}</a></li>')
        contents = ''.join(
            f'<li class="{kind}"><a href="#{escape(target)}" title="{escape(label)}">'
            f'{escape(label)}</a></li>' for target, kind, label in outline.items)
        landing = path.stem == 'index'
        title = 'MOM6 API reference' if landing else labels[path.stem]
        page = template.substitute(
            title=escape(title), body=parts['html_body'], modules='\n'.join(navigation),
            contents=contents, page_class='module-index' if landing else 'module-page',
            breadcrumb='Module index' if landing else escape(labels[path.stem]),
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
