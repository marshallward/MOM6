"""Locate maintained science chapters alongside their implementation sources."""
import os
import re


def science_sources(docs):
    pages = {}
    for path in sorted((docs.parent / 'src').rglob('*.rst')):
        if path.stem in pages:
            raise ValueError(f'Duplicate science page {path.stem}: {pages[path.stem]} and {path}')
        pages[path.stem] = path
    return pages


def science_text(path, docs, image_prefix):
    """Adjust figure paths when staging a chapter away from its source directory."""
    prefix = os.path.relpath(docs / 'images', path.parent) + '/'
    text = path.read_text(encoding='utf-8')
    return re.sub(r'(?m)^(\s*\.\. (?:figure|image):: )' + re.escape(prefix),
                  lambda match: match[1] + image_prefix, text)
