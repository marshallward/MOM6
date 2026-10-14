#!/usr/bin/env python3
"""Report Fortran operator spacing inconsistencies.

The default mode is a dry run over explicit source paths. Binary operators
require one or more spaces on each side. Unary signs require no following
space, and ``.not.`` requires a following space. Concatenation operators,
exponentiation operators, arithmetic inside callable or subscript groups, and
do-loop controls, and ``operator(...)`` generic specifiers are left unchanged.

Options:
  ``--report`` prints each changed source line with its file and line number.
  ``--diff`` prints unified diffs for proposed changes.
  ``--write`` applies changes in place.
  ``--exclude`` skips matching paths when parsing directory trees.
  ``--strict`` requires exactly one space around binary operators.
"""

from __future__ import annotations

import argparse
import difflib
from pathlib import Path
from typing import NamedTuple

import flint
from flint.formatting.spacing import format_spacing, spacing_edits


class Change(NamedTuple):
    """A line diagnostic with its source path attached for reporting."""

    path: Path
    line_number: int
    original: str


SKIPPED_OPERATOR_ROLES = {
    'callable_or_subscript_operator',
    'concatenation',
    'declaration_specifier',
    'do_control_operator',
    'exponentiation',
    'format_specifier',
    'generic_operator_specifier',
}


def collect_spacing_edits(statement, strict=False):
    """Return operator spacing proposals without changing parsed tokens."""
    edits = []
    for token in statement.tokens:
        role = getattr(token, 'operator_role', None)
        if role is None or role in SKIPPED_OPERATOR_ROLES:
            continue

        if role == 'unary_operator' and token == '.not.':
            edits.extend(spacing_edits(token, after=' ', minimum=not strict))
        elif role == 'unary_operator':
            edits.extend(spacing_edits(token, after=''))
        elif role == 'binary_operator':
            edits.extend(spacing_edits(token, ' ', ' ', minimum=not strict))
    return edits


def format_statements(
    path: Path,
    statements,
    source_lines=None,
    strict=False,
) -> tuple[str, list[Change]]:
    """Format parsed statements and return rendered source plus changes."""
    formatted, changes = format_spacing(
        statements,
        lambda statement: collect_spacing_edits(statement, strict=strict),
        source_lines,
    )
    return formatted, [
        Change(path, change.line_number, change.original) for change in changes
    ]


def report_operator_spacing(
    paths,
    write=False,
    report=False,
    diff=False,
    excludes=(),
    strict=False,
) -> int:
    """Report or fix operator spacing in parsed flint sources."""
    changed: list[Path] = []
    all_changes: list[Change] = []
    project = flint.parse(*(str(path) for path in paths), excludes=excludes)

    for source in project.sources:
        path = Path(source.path)
        old = path.read_text(encoding='utf-8', errors='ignore')
        source_lines = old.splitlines() if report else None

        new, file_changes = format_statements(
            path,
            source.statements,
            source_lines,
            strict=strict,
        )

        if new != old:
            changed.append(path)
            all_changes.extend(file_changes)

            if diff:
                print(
                    ''.join(
                        difflib.unified_diff(
                            old.splitlines(True),
                            new.splitlines(True),
                            fromfile=str(path),
                            tofile=str(path),
                        )
                    ),
                    end='',
                )

            if write:
                path.write_text(new, encoding='utf-8')

    if report:
        seen = set()
        for change in all_changes:
            key = (change.path, change.line_number)
            if key in seen:
                continue
            seen.add(key)
            print(f'{change.path}:{change.line_number}: {change.original}')
    elif not diff:
        for path in changed:
            print(path)

    return 1 if changed and not write else 0


def main() -> int:
    """Parse command-line arguments and run the operator checker."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='+', type=Path, help='Files to process')
    parser.add_argument(
        '--write', action='store_true', help='Write changes in place',
    )
    parser.add_argument(
        '--report',
        action='store_true',
        help='Report changed statements with line numbers',
    )
    parser.add_argument(
        '--diff',
        action='store_true',
        help='Print unified diffs for proposed changes',
    )
    parser.add_argument(
        '--strict',
        action='store_true',
        help='Require exactly one space around binary operators',
    )
    parser.add_argument(
        '--exclude',
        action='append',
        default=[],
        help='Directory to exclude; may be repeated',
    )

    args = parser.parse_args()
    return report_operator_spacing(
        args.paths,
        write=args.write,
        report=args.report,
        diff=args.diff,
        excludes=args.exclude,
        strict=args.strict,
    )


if __name__ == '__main__':
    raise SystemExit(main())
