#!/usr/bin/env python3
"""Build the experimental Flint API reference, then validate local HTML links."""
import argparse
import contextlib
from dataclasses import replace
from html.parser import HTMLParser
import json
import tomllib
from pathlib import Path
from urllib.parse import unquote, urlsplit

import flint
from flint.docgen.model import iter_module_documentation
from flint.docgen.references import HtmlDestinations, ReferenceIndex
from flint.docgen.rst import write_module
from render_api_docs import render_reference


DEFAULT_SOURCES = (
    'src',
    'config_src/drivers/solo_driver',
    'config_src/external',
    'config_src/infra/FMS2',
    'config_src/memory/dynamic_symmetric',
)


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
    sources = [source] if isinstance(source, Path) else list(source)
    for root in sources:
        if not root.is_dir():
            raise FileNotFoundError(f'Source directory does not exist: {root}')
    paths = sorted({path for root in sources for path in root.rglob('*')})
    missing = [path for path in paths if path.is_symlink() and not path.exists()]
    if missing:
        raise FileNotFoundError(
            'Source symlinks have missing targets. Initialize the MOM6 submodules with '
            '`git submodule update --init --recursive`:\n' + '\n'.join(map(str, missing)))
    files = [path for path in paths if path.suffix in ('.F90', '.f90')]
    if not files:
        raise ValueError(f'No Fortran files found in {source}')
    includes = sorted({str(path.parent) for path in paths if path.is_file()})
    modules, externals = [], []
    owners = {}
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
                        name = module.name.lower()
                        if name in owners:
                            raise ValueError(f'Duplicate module {module.name}: {owners[name]} and {path}')
                        owners[name] = path
                        modules.append(module)
                del project
    if externals:
        modules.append(replace(external, procedures=tuple(externals)))
    dependency_file = Path(__file__).with_name('api-dependencies.toml')
    with dependency_file.open('rb') as stream:
        dependencies = tomllib.load(stream)['libraries']
    external_modules = {name: library for library, names in dependencies.items() for name in names}
    index = ReferenceIndex(modules, external_modules=external_modules)
    (reports / 'unresolved.txt').write_text(
        ''.join(message + '\n' for message in index.diagnostics), encoding='utf-8')
    external_counts = {}
    with (reports / 'external.txt').open('w', encoding='utf-8') as stream:
        for scope, calls in index.calls.items():
            for ref in calls:
                if ref.external is not None:
                    symbol = ref.external
                    external_counts[symbol.library] = external_counts.get(symbol.library, 0) + 1
                    stream.write(f'{scope}: {ref.label}: {symbol.library}: {symbol.module}::{symbol.name}\n')
    destinations = HtmlDestinations(index)
    for module in modules:
        with (rst_dir / (module.name + '.rst')).open('w', encoding='utf-8') as stream:
            write_module(stream, module, references=index, destinations=destinations)
    with (rst_dir / 'index.rst').open('w', encoding='utf-8') as stream:
        stream.write('MOM6 API reference\n==================\n\n')
        stream.write('Experimental reference generated from these source directories:\n\n')
        for root in sources:
            stream.write(f'* ``{root}``\n')
        stream.write('\nExternal calls are grouped by dependency. Unknown references are retained '
                     'in build diagnostics, not published as calls.\n\n'
                     'Modules\n-------\n\n')
        for module in sorted(modules, key=lambda entry: entry.name.lower()):
            stream.write(f'* `{module.name} <{module.name}.html>`__\n')
    return {'source_files': len(files), 'module_pages': len(modules),
            'call_references': sum(len(calls) for calls in index.calls.values()),
            'external_references': sum(external_counts.values()),
            'external_libraries': external_counts,
            'unresolved_references': len(index.diagnostics)}


def build(source, output):
    # A fresh directory prevents old module pages from masking broken links.
    output.mkdir(parents=True, exist_ok=False)
    rst_dir, html_dir, reports = (output / name for name in ('rst', 'html', 'reports'))
    for directory in (rst_dir, html_dir, reports):
        directory.mkdir()
    summary = generate_reference(source, rst_dir, reports)
    with (reports / 'render.log').open('w', encoding='utf-8') as warnings:
        render_reference(rst_dir, html_dir, warnings)
    links, cross_file = check_links(html_dir)
    summary.update(html_pages=len(list(html_dir.glob('*.html'))),
                   checked_links=links, cross_file_links=cross_file)
    (html_dir / '.nojekyll').touch()
    (reports / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(summary, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, action='append',
                        help='Source directory; repeat to combine roots (replaces the MOM6 defaults)')
    parser.add_argument('--output', type=Path, default=Path('docs/_build/flint-api'))
    args = parser.parse_args()
    build(args.source or [Path(root) for root in DEFAULT_SOURCES], args.output)


if __name__ == '__main__':
    main()
