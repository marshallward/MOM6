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


SKIPPED_OPERATOR_ROLES = {
    'callable_or_subscript_operator',
    'concatenation',
    'declaration_specifier',
    'do_control_operator',
    'exponentiation',
    'format_specifier',
    'generic_operator_specifier',
}


class Change(NamedTuple):
    """Description of one source line with operator spacing changes."""

    path: Path
    line_number: int
    original: str


def simple_liminals(liminals: list[str]) -> bool:
    """Return True if liminals are simple same-line whitespace."""
    return all(
        (item == '' or item.isspace()) and '\n' not in item
        for item in liminals
    )


def set_spacing(token, before: str | None, after: str | None) -> None:
    """Set spacing around a token when current liminals are simple."""
    if before is not None and simple_liminals(token.head):
        token.head[:] = [before] if before else []
    if after is not None and simple_liminals(token.tail):
        token.tail[:] = [after] if after else []


def set_minimum_spacing(token, before: str | None, after: str | None) -> None:
    """Ensure at least the requested spacing around a token."""
    if (
        before is not None
        and simple_liminals(token.head)
        and not ''.join(token.head)
    ):
        token.head[:] = [before]
    if (
        after is not None
        and simple_liminals(token.tail)
        and not ''.join(token.tail)
    ):
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


def apply_spacing_policy(statement, strict=False) -> None:
    """Apply operator spacing policy to one statement."""
    for token in statement.tokens:
        role = getattr(token, 'operator_role', None)
        if role is None or role in SKIPPED_OPERATOR_ROLES:
            continue

        if role == 'unary_operator' and token == '.not.':
            set_spacing(token, None, ' ') if strict else set_minimum_spacing(
                token, None, ' ',
            )
        elif role == 'unary_operator':
            set_spacing(token, None, '')
        elif role == 'binary_operator' and strict:
            set_spacing(token, ' ', ' ')
        elif role == 'binary_operator':
            set_minimum_spacing(token, ' ', ' ')


def changed_lines(statement, original, formatted, source_lines=None):
    """Return physical source lines changed by formatting."""
    original_lines = original.splitlines()
    formatted_lines = formatted.splitlines()

    for offset, original_line in enumerate(original_lines):
        if (
            offset >= len(formatted_lines)
            or original_line != formatted_lines[offset]
        ):
            line_number = statement.line_number + offset
            if source_lines and line_number <= len(source_lines):
                yield line_number, source_lines[line_number - 1]
            else:
                yield line_number, original_line


def format_statements(
    path: Path,
    statements,
    source_lines=None,
    strict=False,
) -> tuple[str, list[Change]]:
    """Format parsed statements and return rendered source plus changes."""
    changes: list[Change] = []

    for statement in statements:
        if not getattr(statement, 'source_visible', True):
            continue

        original_statement = statement.source_text()
        if not original_statement:
            continue

        original_body = statement.source_line()
        apply_spacing_policy(statement, strict=strict)

        formatted_statement = statement.source_text()
        if formatted_statement != original_statement:
            for line_number, original_line in changed_lines(
                statement,
                original_body,
                statement.source_line(),
                source_lines,
            ):
                changes.append(
                    Change(
                        path=path,
                        line_number=line_number,
                        original=original_line,
                    )
                )

    return render(statements), changes


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
