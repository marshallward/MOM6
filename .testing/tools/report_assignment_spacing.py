#!/usr/bin/env python3
"""Report Fortran assignment spacing inconsistencies.

The default mode is a dry run over explicit source paths. Assignment operators
require one or more spaces on each side. Name-value pairs, such as named
arguments and declaration specifiers, require no spaces around the equals sign.

Options:
  ``--report`` prints each changed source line with its file and line number.
  ``--diff`` prints unified diffs for proposed changes.
  ``--write`` applies changes in place.
  ``--exclude`` skips matching paths when parsing directory trees.
  ``--strict`` requires exactly one space around assignment operators.
  ``--do-controls`` also checks do-loop controls such as ``do i=1,n``.
"""

from __future__ import annotations

import argparse
import difflib
from pathlib import Path
from typing import NamedTuple

import flint


class Change(NamedTuple):
    """Description of one source line with assignment spacing changes."""

    path: Path
    line_number: int
    original: str


def simple_liminals(liminals: list[str]) -> bool:
    """Return True if liminals are simple same-line whitespace."""
    return all(
        (item == '' or item.isspace()) and '\n' not in item
        for item in liminals
    )


def set_spacing(token, before: str, after: str) -> None:
    """Set spacing around a token when current liminals are simple."""
    if simple_liminals(token.head):
        token.head[:] = [before] if before else []
    if simple_liminals(token.tail):
        token.tail[:] = [after] if after else []


def set_minimum_spacing(token, before: str, after: str) -> None:
    """Ensure at least the requested spacing around a token."""
    if simple_liminals(token.head) and not ''.join(token.head):
        token.head[:] = [before]
    if simple_liminals(token.tail) and not ''.join(token.tail):
        token.tail[:] = [after]


def render(statements) -> str:
    """Render a sequence of flint statements back to source text."""
    output: list[str] = []
    first = True
    for statement in statements:
        if first:
            output.append(''.join(statement[0].head))
            first = False
        for token in statement:
            output.append(str(token))
            output.append(''.join(token.tail))
    return ''.join(output)


def apply_spacing_policy(statement, do_controls=False, strict=False) -> None:
    """Apply MOM6 assignment spacing policy to one statement."""
    for token in statement.tokens:
        if token.syntax_role == 'assignment':
            if strict:
                set_spacing(token, ' ', ' ')
            else:
                set_minimum_spacing(token, ' ', ' ')
        elif token.syntax_role == 'name_value':
            set_spacing(token, '', '')
        elif token.syntax_role == 'do_control' and do_controls:
            set_spacing(token, '', '')

    if do_controls and statement.is_do_statement():
        for token in statement:
            if str(token) == ',' and simple_liminals(token.tail):
                token.tail[:] = []


def changed_line(statement, original, formatted, source_lines=None):
    """Return the first physical source line changed by formatting."""
    original_lines = original.splitlines()
    formatted_lines = formatted.splitlines()

    for offset, original_line in enumerate(original_lines):
        if (
            offset >= len(formatted_lines)
            or original_line != formatted_lines[offset]
        ):
            line_number = statement.line_number + offset
            if source_lines and line_number <= len(source_lines):
                return line_number, source_lines[line_number - 1]
            return line_number, original_line

    return statement.line_number, original_lines[0] if original_lines else ''


def format_statements(
    path: Path,
    statements,
    source_lines=None,
    do_controls=False,
    strict=False,
) -> tuple[str, list[Change]]:
    """Format parsed statements and return rendered source plus changes."""
    changes: list[Change] = []

    for statement in statements:
        original_statement = statement.source_text()
        original_body = statement.source_line()

        apply_spacing_policy(statement, do_controls=do_controls, strict=strict)

        formatted_statement = statement.source_text()
        if formatted_statement != original_statement:
            line_number, original_line = changed_line(
                statement,
                original_body,
                statement.source_line(),
                source_lines,
            )
            changes.append(
                Change(
                    path=path,
                    line_number=line_number,
                    original=original_line,
                )
            )

    return render(statements), changes


def report_assignment_spacing(
    paths,
    write=False,
    report=False,
    diff=False,
    excludes=(),
    do_controls=False,
    strict=False,
) -> int:
    """Report or fix assignment spacing in parsed flint sources."""
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
            do_controls=do_controls,
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
    """Parse command-line arguments and run the assignment spacing checker."""
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
        '--do-controls',
        action='store_true',
        help='Also check do-loop control spacing',
    )
    parser.add_argument(
        '--strict',
        action='store_true',
        help='Require exactly one space around assignment operators',
    )
    parser.add_argument(
        '--exclude',
        action='append',
        default=[],
        help='Directory to exclude; may be repeated',
    )

    args = parser.parse_args()
    return report_assignment_spacing(
        args.paths,
        write=args.write,
        report=args.report,
        diff=args.diff,
        excludes=args.exclude,
        do_controls=args.do_controls,
        strict=args.strict,
    )


if __name__ == '__main__':
    raise SystemExit(main())
