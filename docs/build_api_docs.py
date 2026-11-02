#!/usr/bin/env python3
"""Build the experimental Flint API reference, then validate local HTML links."""
import argparse
import contextlib
from dataclasses import replace
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import unquote, urlsplit

from docutils.core import publish_file
import flint
from flint.docgen.model import iter_module_documentation
from flint.docgen.references import HtmlDestinations, ReferenceIndex
from flint.docgen.rst import write_module


class HtmlLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.links = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'id' in attrs:
            self.ids.add(attrs['id'])
        if tag == 'a' and 'href' in attrs:
            self.links.append(attrs['href'])


def check_links(directory):
    """Check relative page links and fragments in the rendered site."""
    pages = {}
    for path in sorted(directory.glob('*.html')):
        page = HtmlLinks()
        page.feed(path.read_text(encoding='utf-8'))
        pages[path.name] = page
    count = cross_file = 0
    errors = []
    for filename, page in pages.items():
        for href in page.links:
            url = urlsplit(href)
            if url.scheme or url.netloc:
                continue
            target = unquote(url.path) or filename
            fragment = unquote(url.fragment)
            count += 1
            cross_file += target != filename
            if target not in pages or (fragment and fragment not in pages[target].ids):
                errors.append(f'{filename}: {href}')
    if errors:
        raise ValueError('Broken HTML links:\n' + '\n'.join(errors))
    return count, cross_file


def generate_reference(source, rst_dir, reports):
    """Parse one file at a time; retain only documentation snapshots for linking."""
    missing = sorted(path for path in source.rglob('*') if path.is_symlink() and not path.exists())
    if missing:
        raise FileNotFoundError(
            'Source symlinks have missing targets. Initialize the MOM6 submodules with '
            '`git submodule update --init --recursive`:\n' + '\n'.join(map(str, missing)))
    files = sorted(path for path in source.rglob('*') if path.suffix in ('.F90', '.f90'))
    if not files:
        raise ValueError(f'No Fortran files found in {source}')
    includes = sorted({str(path.parent) for path in source.rglob('*') if path.is_file()})
    modules, externals = [], []
    external = None
    with (reports / 'parse.log').open('w', encoding='utf-8') as log:
        for path in files:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                project = flint.parse(str(path), includes=includes)
                if len(project.sources) != 1:
                    raise ValueError(f'Unable to parse {path}; see {reports / "parse.log"}')
                for module in iter_module_documentation(project):
                    if module.name == 'External' and not module.kind:
                        external = module
                        externals.extend(module.procedures)
                    else:
                        modules.append(module)
                del project
    if externals:
        modules.append(replace(external, procedures=tuple(externals)))
    index = ReferenceIndex(modules)
    (reports / 'unresolved.txt').write_text(
        ''.join(message + '\n' for message in index.diagnostics), encoding='utf-8')
    destinations = HtmlDestinations(index)
    for module in modules:
        with (rst_dir / (module.name + '.rst')).open('w', encoding='utf-8') as stream:
            write_module(stream, module, references=index, destinations=destinations)
    with (rst_dir / 'index.rst').open('w', encoding='utf-8') as stream:
        stream.write('MOM6 API reference\n==================\n\n')
        stream.write('Experimental reference generated from the MOM6 ``src`` tree.\n'
                     'Unresolved procedure references are displayed as plain text.\n\n')
        for module in sorted(modules, key=lambda entry: entry.name.lower()):
            stream.write(f'* `{module.name} <{module.name}.html>`__\n')
    return {'source_files': len(files), 'module_pages': len(modules),
            'call_references': sum(len(calls) for calls in index.calls.values()),
            'unresolved_references': len(index.diagnostics)}


def build(source, output):
    # A fresh directory prevents old module pages from masking broken links.
    output.mkdir(parents=True, exist_ok=False)
    rst_dir, html_dir, reports = (output / name for name in ('rst', 'html', 'reports'))
    for directory in (rst_dir, html_dir, reports):
        directory.mkdir()
    summary = generate_reference(source, rst_dir, reports)
    with (reports / 'render.log').open('w', encoding='utf-8') as warnings:
        for rst in sorted(rst_dir.glob('*.rst')):
            publish_file(source_path=str(rst),
                         destination_path=str(html_dir / (rst.stem + '.html')),
                         writer='html5', settings_overrides={
                             'warning_stream': warnings, 'halt_level': 2,
                             'input_encoding': 'utf-8', 'output_encoding': 'utf-8',
                         })
    links, cross_file = check_links(html_dir)
    summary.update(html_pages=len(list(html_dir.glob('*.html'))),
                   checked_links=links, cross_file_links=cross_file)
    (html_dir / '.nojekyll').touch()
    (reports / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('src'))
    parser.add_argument('--output', type=Path, default=Path('docs/_build/flint-api'))
    args = parser.parse_args()
    build(args.source.resolve(), args.output)


if __name__ == '__main__':
    main()
