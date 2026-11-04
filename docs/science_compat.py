"""Expose maintained science reST at the legacy manual's generated paths."""
from pathlib import Path
import re
from science_sources import science_sources, science_text


def write_doxygen_targets(docs):
    """Give Fortran comments link destinations, without duplicating science prose."""
    destination = docs / '_build/science-targets.dox'
    destination.parent.mkdir(parents=True, exist_ok=True)
    blocks = []
    for path in science_sources(docs).values():
        text = path.read_text(encoding='utf-8')
        title = text.splitlines()[2]
        targets = re.findall(r'^\.\. _([^\n]+):$', text, re.MULTILINE)
        blocks.append('/*! \\page ' + path.stem + ' ' + title + '\n'
                      + ''.join('\\anchor ' + target + '\n' for target in targets
                                if target != path.stem) + '*/\n')
    destination.write_text('\n'.join(blocks), encoding='utf-8')
    return destination


def copy_science_pages(app):
    """Run after Doxygen autosummary so the maintained sources take precedence."""
    docs = Path(app.srcdir)
    destination = docs / 'api/generated/pages'
    destination.mkdir(parents=True, exist_ok=True)
    for path in science_sources(docs).values():
        text = science_text(path, docs, '../../../images/')
        # Doxygen's references to section anchors include the owning page ID.
        text = re.sub(r'^\.\. _([^\n]+):$',
                      lambda m: (m[0] if m[1] == path.stem else
                                 f'.. _{path.stem}_1{m[1]}:\n\n' + m[0]),
                      text, flags=re.MULTILINE)
        (destination / path.name).write_text(text, encoding='utf-8')
