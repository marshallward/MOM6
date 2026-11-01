#!/usr/bin/env python3
"""Report real variables without bracketed units in associated documentation.

Flint supplies declared variables, comment associations, and source locations.
MOM6 accepts bracketed unit expressions such as [m] or [L T-1 ~> m s-1],
and reports competing forward/backward Doxygen documentation.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import re
from typing import NamedTuple

import flint


UNIT_RE = re.compile(r'\[[^\]\n]+\]')


class Issue(NamedTuple):
    path: Path
    line_number: int
    name: str
    original: str
    comments: tuple[str, ...]
    message: str = 'missing bracketed units'


def collect_units(unit):
    """Visit implementation scopes and explicit interface bodies."""
    yield from unit.walk()


def real_unit_issues(source) -> list[Issue]:
    """Apply MOM6's units-documentation policy to Flint variables."""
    path = Path(source.path)
    lines = path.read_text(encoding='utf-8', errors='ignore').splitlines()
    issues = []
    for root in source.units:
        for unit in collect_units(root):
            for variable in unit.variables:
                if variable.type != 'real' or not getattr(variable.stmt, 'source_visible', True):
                    continue
                location = variable.location
                if location is None:
                    continue
                documentation = variable.documentation
                if documentation is None:
                    selected = []
                    message = 'unsupported declaration documentation'
                else:
                    direct = [block for block in documentation if block.placement == 'trailing']
                    selected = direct + [block for block in documentation if block.placement == 'leading']
                    if not direct:
                        selected += [block for block in documentation if block.placement == 'following']
                    kinds = {block.kind for block in documentation}
                    competing = 'forward' in kinds and 'backward' in kinds
                    if not competing and any(UNIT_RE.search(' '.join(block.comments)) for block in selected):
                        continue
                    message = 'competing Doxygen comments' if competing else 'missing bracketed units'
                issues.append(Issue(
                    path, location.line, str(variable.name), lines[location.line - 1],
                    selected[0].comments if selected else (), message,
                ))
    return issues


def report_real_units(paths, report=False, excludes=()) -> int:
    project = flint.parse(*(str(path) for path in paths), excludes=excludes)
    found = False
    for source in project.sources:
        issues = real_unit_issues(source)
        found |= bool(issues)
        if not report:
            if issues:
                print(source.path)
            continue
        for issue in issues:
            print(f'{issue.path}:{issue.line_number}: {issue.name}: {issue.message}')
            print(f"  declaration: {issue.original.split('!', 1)[0].rstrip()}")
            if issue.comments:
                print('  comments:')
                for comment in issue.comments:
                    print(f'    {comment}')
            else:
                print('  comments: <no associated comment>')
    return int(found)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='+', type=Path, help='Files to process')
    parser.add_argument('--report', action='store_true', help='Report variables missing units with line numbers')
    parser.add_argument('--exclude', action='append', default=[], help='Directory to exclude; may be repeated')
    args = parser.parse_args()
    return report_real_units(args.paths, report=args.report, excludes=args.exclude)


if __name__ == '__main__':
    raise SystemExit(main())
