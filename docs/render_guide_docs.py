"""Render the handwritten manual with Sphinx, then reuse the MOM6 site theme."""
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

from sphinx.util.inventory import InventoryFile
from science_sources import science_sources


GUIDE_PAGES = {
    'about': 'About this documentation',
    'equations': 'Equations',
    'discrete_space': 'Spatial Discretization',
    'discrete_time': 'Time Discretization',
    'tracers': 'Tracers in MOM6',
    'grids': 'Grids',
    'parameterizations': 'Parameterizations',
    'parameterizations_vertical': 'Vertical Parameterizations',
    'parameterizations_lateral': 'Lateral Parameterizations',
    'other_physics': 'Other Physics',
    'working_with_MOM6': 'Working with MOM6',
    'code_organization': 'Organization of the code',
    'forcing': 'Forcing',
    'parallel': 'Parallel Implementation',
    'testing': 'Testing of MOM6',
    'zzbibliography': 'Bibliography',
}
MANUAL_URL = 'https://mom6.readthedocs.io/en/main/'
DOCS = Path(__file__).resolve().parent
SCIENCE_SOURCES = science_sources(DOCS)
SCIENCE_PAGES = {name: path.read_text(encoding='utf-8').splitlines()[2]
                 for name, path in SCIENCE_SOURCES.items()}


def external_chapters(source, inventory):
    """Resolve legacy chapter paths to local pages or published destinations."""
    def replace(match):
        name = match[2].removesuffix('.rst')
        if name.startswith('api/generated/pages/') and name.rsplit('/', 1)[-1] in SCIENCE_PAGES:
            return match[1] + name.rsplit('/', 1)[-1]
        entry = inventory['std:doc'].get(name)
        if entry is None:
            raise ValueError(f'No published manual destination for {name}')
        return f'{match[1]}{entry.display_name} <{entry.uri}>'
    return re.sub(r'^([ \t]+)(api/generated/\S+)[ \t]*$', replace, source, flags=re.MULTILINE)


def render_guides(rst_dir, warning_stream, html_dir=None):
    """Return HTML fragments; Sphinx handles toctrees, citations and references.

    Build in a subprocess so Sphinx's Docutils registrations cannot leak into
    the standalone Flint renderer. The checked-in inventory makes reference
    resolution independent of network availability during builds.
    """
    names = [name for name in {**GUIDE_PAGES, **SCIENCE_PAGES}
             if (rst_dir / (name + '.rst')).is_file()]
    if not names:
        return {}
    inventory = InventoryFile.load(io.BytesIO((DOCS / 'guide-inventory.inv').read_bytes()),
                                   MANUAL_URL, lambda base, path: base + path)
    with tempfile.TemporaryDirectory(prefix='mom6-guides-') as directory:
        root = Path(directory)
        source, output = root / 'source', root / 'html'
        source.mkdir()
        shutil.copyfile(DOCS / 'guide_sphinx/conf.py', source / 'conf.py')
        shutil.copyfile(DOCS / 'guide-inventory.inv', source / 'guide-inventory.inv')
        shutil.copytree(DOCS / 'images', source / 'images')
        for name in names:
            text = (rst_dir / (name + '.rst')).read_text(encoding='utf-8')
            (source / (name + '.rst')).write_text(external_chapters(text, inventory), encoding='utf-8')
        for name in ('ocean.bib', 'references.bib', 'zotero.bib'):
            shutil.copyfile(DOCS / name, source / name)
        (source / 'guide-index.rst').write_text(
            'Manual\n======\n\n.. toctree::\n\n' + ''.join(f'   {name}\n' for name in names),
            encoding='utf-8')
        result = subprocess.run(
            [sys.executable, '-m', 'sphinx', '-b', 'json', '-q', '-W', '--keep-going',
             str(source), str(output)],
            capture_output=True, text=True)
        warning_stream.write(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError('Manual rendering failed:\n' + result.stdout + result.stderr)
        if html_dir is not None and (output / '_images').exists():
            shutil.copytree(output / '_images', html_dir / '_images', dirs_exist_ok=True)
        return {name: json.loads((output / (name + '.fjson')).read_text(encoding='utf-8'))['body']
                for name in names}
