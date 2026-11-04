#!/usr/bin/env python3
"""Build the experimental Flint API reference, then validate local HTML links."""
import argparse
import contextlib
from dataclasses import replace
from html.parser import HTMLParser
import json
import shutil
import subprocess
import tomllib
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.parse import unquote, urlsplit

import flint
from flint.docgen.model import iter_module_documentation
from flint.docgen.references import HtmlDestinations, Reference, ReferenceIndex, anchor
from flint.docgen.graphs import call_graph, write_dot
from flint.docgen.rst import write_module
from render_api_docs import GUIDE_PAGES, render_reference
from render_guide_docs import SCIENCE_SOURCES
from science_sources import science_text


DEFAULT_SOURCES = (
    'src',
    'config_src/drivers/solo_driver',
    'config_src/external',
    'config_src/infra/FMS2',
    'config_src/memory/dynamic_symmetric',
)
REPOSITORY = Path(__file__).resolve().parents[1]
DEPENDENCY_CONFIG = Path(__file__).with_name('api-dependencies.toml')
DEPENDENCY_SOURCES = {'GSW': 'pkg/GSW-Fortran', 'CVMix': 'pkg/CVMix-src'}


class HtmlLinks(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.links = []

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'id' in attrs:
            self.ids.add(attrs['id'])
        if tag == 'a':
            href = attrs.get('href', attrs.get('xlink:href'))
            if href is not None:
                self.links.append(href)
        if tag in ('img', 'script') and attrs.get('src'):
            self.links.append(attrs['src'])
        if tag == 'link' and attrs.get('href'):
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
            if (not (directory / target).is_file()
                    or (fragment and (target not in pages or fragment not in pages[target].ids))):
                errors.append(f'{filename}: {href}')
    if errors:
        raise ValueError('Broken HTML links:\n' + '\n'.join(errors))
    return count, cross_file


def generate_graphs(index, directory, depth=1):
    """Render directional graphs, including private implementation procedures."""
    dot = shutil.which('dot')
    if dot is None:
        raise RuntimeError('Call graphs require Graphviz; install the dot executable.')
    images = {}
    ET.register_namespace('', 'http://www.w3.org/2000/svg')
    ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')
    for key, unit in sorted(index.units.items()):
        if unit.kind not in ('function', 'subroutine', 'interface'):
            continue
        for direction in ('callees', 'callers'):
            graph = call_graph(index, key, depth=depth, direction=direction)
            if not graph.edges:
                continue
            path = directory / (direction + '-' + anchor(key) + '.dot')
            with path.open('w', encoding='utf-8') as stream:
                write_dot(stream, graph, index, html_links=True)
            image = path.with_suffix('.svg')
            # Graphviz uses points: 12pt matches the site's 16px reading size.
            subprocess.run([dot, '-Tsvg', '-Nfontsize=12', '-Gfontsize=12',
                            '-Efontsize=12', str(path), '-o', str(image)], check=True)
            # Graphviz repeats graph0/node1/etc. across files. Inline SVGs share
            # the HTML document, so give their internal IDs unique prefixes.
            tree = ET.parse(image)
            for node in tree.iter():
                if 'id' in node.attrib:
                    node.set('id', path.stem + '-' + node.get('id'))
                for attr in ('href', '{http://www.w3.org/1999/xlink}href'):
                    value = node.get(attr, '')
                    if value.startswith('#'):
                        node.set(attr, '#' + path.stem + '-' + value[1:])
            tree.write(image, encoding='utf-8', xml_declaration=True)
            images.setdefault(key, []).append((direction.title(), image.name))
    return images


def generate_reference(source, rst_dir, reports, *, graph_depth=1):
    """Parse one file at a time; retain only documentation snapshots for linking."""
    dependency_file = DEPENDENCY_CONFIG
    if dependency_file.exists():
        with dependency_file.open('rb') as stream:
            dependencies = tomllib.load(stream)
    else:
        dependencies = {}
    library_roots = {library: (REPOSITORY / path).resolve()
                     for library, path in dependencies.get('source_libraries', DEPENDENCY_SOURCES).items()}
    sources = [source] if isinstance(source, Path) else list(source)
    for root in sources:
        if not root.is_dir():
            raise FileNotFoundError(f'Source directory does not exist: {root}')
    paths = sorted({path for root in sources for path in root.rglob('*')})
    project_paths = []
    excluded_counts = {}
    with (reports / 'excluded-dependencies.txt').open('w', encoding='utf-8') as stream:
        for path in paths:
            target = path.resolve()
            library = next((name for name, root in library_roots.items()
                            if target.is_relative_to(root)), None)
            if library is None:
                project_paths.append(path)
            elif path.suffix in ('.F90', '.f90'):
                excluded_counts[library] = excluded_counts.get(library, 0) + 1
                stream.write(f'{library}: {path} -> {target}\n')
    paths = project_paths
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
    external_modules = {name: library for library, names in dependencies.get('libraries', {}).items()
                        for name in names}
    index = ReferenceIndex(modules, external_modules=external_modules)
    (reports / 'unresolved.txt').write_text(
        ''.join(message + '\n' for message in index.diagnostics), encoding='utf-8')
    external_counts = {}
    potential_external = (sum(bool(ref.unindexed_module) for calls in index.calls.values() for ref in calls)
                          if hasattr(Reference, 'unindexed_module') else None)
    if potential_external:
        print(f'Warning: {potential_external} call references import symbols from unindexed modules; '
              f'they may be external dependencies or omitted project sources. See {reports / "unresolved.txt"}.')
    elif potential_external is None and index.diagnostics:
        print(f'Warning: {len(index.diagnostics)} unresolved call references; some may be external '
              'dependencies. This Flint revision does not report unindexed-module provenance. '
              f'See {reports / "unresolved.txt"}.')
    with (reports / 'external.txt').open('w', encoding='utf-8') as stream:
        for scope, calls in index.calls.items():
            for ref in calls:
                if ref.external is not None:
                    symbol = ref.external
                    external_counts[symbol.library] = external_counts.get(symbol.library, 0) + 1
                    stream.write(f'{scope}: {ref.label}: {symbol.library}: {symbol.module}::{symbol.name}\n')
    destinations = HtmlDestinations(index)
    graph_images = generate_graphs(index, rst_dir, depth=graph_depth)
    for module in modules:
        with (rst_dir / (module.name + '.rst')).open('w', encoding='utf-8') as stream:
            write_module(stream, module, references=index, destinations=destinations,
                         graph_images=graph_images)
    shutil.copyfile(Path(__file__).with_name('api-home.rst'), rst_dir / 'index.rst')
    for name in GUIDE_PAGES:
        shutil.copyfile(Path(__file__).with_name(name + '.rst'), rst_dir / (name + '.rst'))
    for name, path in SCIENCE_SOURCES.items():
        text = science_text(path, Path(__file__).parent, 'images/')
        (rst_dir / (name + '.rst')).write_text(text, encoding='utf-8')
    shutil.copytree(Path(__file__).parent / 'images', rst_dir / 'images')
    with (rst_dir / 'api-index.rst').open('w', encoding='utf-8') as stream:
        stream.write('MOM6 API reference\n==================\n\n')
        stream.write('Experimental reference generated from these source directories:\n\n')
        for root in sources:
            stream.write(f'* ``{root}``\n')
        stream.write('\nExternal calls are grouped by dependency. Unknown references are retained '
                     'in build diagnostics, not published as calls.\n\n')
        if library_roots:
            stream.write('Dependency implementation sources excluded (including symlinked files): '
                         + ', '.join(library_roots) + '.\n\n')
        stream.write('Modules\n-------\n\n')
        for module in sorted(modules, key=lambda entry: entry.name.lower()):
            stream.write(f'* `{module.display_name or module.name} <{module.name}.html>`__\n')
    return {'source_files': len(files), 'module_pages': len(modules),
            'call_graphs': sum(len(images) for images in graph_images.values()),
            'excluded_dependency_files': excluded_counts,
            'call_references': sum(len(calls) for calls in index.calls.values()),
            'external_references': sum(external_counts.values()),
            'external_libraries': external_counts,
            'potential_external_references': potential_external,
            'unresolved_references': len(index.diagnostics)}


def build(source, output, *, graph_depth=1):
    if graph_depth < 0:
        raise ValueError('Graph depth must be nonnegative')
    # A fresh directory prevents old module pages from masking broken links.
    output.mkdir(parents=True, exist_ok=False)
    rst_dir, html_dir, reports = (output / name for name in ('rst', 'html', 'reports'))
    for directory in (rst_dir, html_dir, reports):
        directory.mkdir()
    summary = generate_reference(source, rst_dir, reports, graph_depth=graph_depth)
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
    parser.add_argument('--graph-depth', type=int, default=1,
                        help='Caller/callee graph depth (default: 1)')
    args = parser.parse_args()
    build(args.source or [Path(root) for root in DEFAULT_SOURCES], args.output,
          graph_depth=args.graph_depth)


if __name__ == '__main__':
    main()
