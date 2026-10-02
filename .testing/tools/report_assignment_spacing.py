#!/usr/bin/env python3
"""Report Fortran assignment spacing inconsistencies using flint statements.

Default behavior is a dry run over explicit ``*.F90`` paths. Use ``--write`` to
apply changes. The default excludes unmaintained cap directories.
"""

from __future__ import annotations

import argparse
import difflib
from pathlib import Path
from typing import NamedTuple

import flint

MOM6_SPACING = {
    'assignment': (' ', ' '),
    'name_value': ('', ''),
    'do_control': ('', ''),
}


class Change(NamedTuple):
    path: Path
    line_number: int
    original: str


def simple_liminals(liminals: list[str]) -> bool:
    return all((item == '' or item.isspace()) and '\n' not in item for item in liminals)


def set_spacing(token, before: str, after: str) -> None:
    if simple_liminals(token.head):
        token.head[:] = [before] if before else []
    if simple_liminals(token.tail):
        token.tail[:] = [after] if after else []


def render(statements) -> str:
    output: list[str] = []
    first = True
    for statement in statements:
        if not statement:
            continue
        if first:
            output.append(''.join(statement[0].head))
            first = False
        for token in statement:
            output.append(str(token))
            output.append(''.join(token.tail))
    return ''.join(output)


def apply_spacing_policy(statement) -> None:
    for token in statement.tokens:
        if token.syntax_role in MOM6_SPACING:
            set_spacing(token, *MOM6_SPACING[token.syntax_role])

    if statement.is_do_statement():
        for token in statement:
            if str(token) == ',' and simple_liminals(token.tail):
                token.tail[:] = []


def changed_line(statement, original, formatted, source_lines=None):
    original_lines = original.splitlines()
    formatted_lines = formatted.splitlines()

    for offset, original_line in enumerate(original_lines):
        if offset >= len(formatted_lines) or original_line != formatted_lines[offset]:
            line_number = statement.line_number + offset
            if source_lines and line_number <= len(source_lines):
                return line_number, source_lines[line_number - 1]
            return line_number, original_line

    return statement.line_number, original_lines[0] if original_lines else ''


def format_statements(path: Path, statements, source_lines=None) -> tuple[str, list[Change]]:
    changes: list[Change] = []

    for statement in statements:
        if not statement:
            continue
        original_statement = statement.source_text()
        if not original_statement:
            continue

        original_body = statement.source_line()

        apply_spacing_policy(statement)

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


def report_assignment_spacing(paths, write=False, report=False, diff=False, excludes=()) -> int:
    changed: list[Path] = []
    all_changes: list[Change] = []
    project = flint.parse(*(str(path) for path in paths), excludes=excludes)

    for source in project.sources:
        path = Path(source.path)
        old = path.read_text(encoding='utf-8', errors='ignore')
        source_lines = old.splitlines() if report else None

        new, file_changes = format_statements(path, source.statements, source_lines)
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='+', type=Path, help='Files to process')
    parser.add_argument('--write', action='store_true', help='Write changes in place')
    parser.add_argument('--report', action='store_true', help='Report changed statements with line numbers')
    parser.add_argument('--diff', action='store_true', help='Print unified diffs for proposed changes')
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
    )


if __name__ == '__main__':
    raise SystemExit(main())
