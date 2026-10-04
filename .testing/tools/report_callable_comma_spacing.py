#!/usr/bin/env python3
"""Report Fortran callable comma spacing inconsistencies.

The default mode is a dry run over explicit source paths. Commas separating
known function or subroutine arguments require one following space. Commas in
known array references, such as ``A(i,j,k)``, are left unchanged.

Options:
  ``--report`` prints each changed source line with its file and line number.
  ``--diff`` prints unified diffs for proposed changes.
  ``--write`` applies changes in place.
  ``--exclude`` skips matching paths when parsing directory trees.
"""

from __future__ import annotations

import argparse
import difflib
from pathlib import Path
from typing import NamedTuple

import flint
from flint.intrinsics import intrinsic_fns


NON_CALLABLE_PAREN_NAMES = {
    'backspace',
    'close',
    'endfile',
    'flush',
    'inquire',
    'open',
    'read',
    'rewind',
    'wait',
    'write',
}


class Change(NamedTuple):
    """Description of one source line with callable comma changes."""

    path: Path
    line_number: int
    original: str


class Context(NamedTuple):
    """Known data objects and callables in scope for a statement."""

    variable_names: set[str]
    callables: set[str]


def simple_liminals(liminals: list[str]) -> bool:
    """Return True if liminals are simple same-line whitespace."""
    return all(
        (item == '' or item.isspace()) and '\n' not in item
        for item in liminals
    )


def set_comma_spacing(token) -> None:
    """Require exactly one same-line space after a comma token."""
    newline_indices = [
        index for index, item in enumerate(token.tail) if '\n' in item
    ]
    if newline_indices:
        newline_index = newline_indices[0]
        prefix = token.tail[:newline_index]
        if '&' in prefix and all(
            item.isspace() or item == '&' for item in prefix
        ):
            token.tail[:] = [' ', '&'] + token.tail[newline_index:]
        return
    if simple_liminals(token.tail):
        token.tail[:] = [' ']


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


def collect_contexts(unit, contexts: dict[int, Context]) -> None:
    """Collect per-statement array and callable context from a parsed unit."""
    for subprogram in unit.subprograms:
        collect_contexts(subprogram, contexts)

    context = Context(
        variable_names={str(var.name).lower() for var in unit.variables},
        callables={str(name).lower() for name in unit.callees},
    )
    for statement in unit.statements:
        contexts.setdefault(id(statement), context)


def statement_contexts(source) -> dict[int, Context]:
    """Return the best known callable context for each source statement."""
    contexts: dict[int, Context] = {}
    for unit in source.units:
        collect_contexts(unit, contexts)
    return contexts


def is_known_callable(name: str, context: Context) -> bool:
    """Return True if name is known to be callable in this context."""
    lowered = name.lower()
    return lowered in context.callables or lowered in intrinsic_fns


def opens_callable_group(statement, index: int, context: Context) -> bool:
    """Return True if parenthesis opens a callable argument list."""
    if statement[index] != '(' or index == 0:
        return False

    if statement.is_call_statement():
        depth = 0
        for tok in statement[statement.code_index + 1:index]:
            if tok == '(':
                depth += 1
            elif tok == ')':
                depth -= 1
        if depth == 0:
            return True

    prior = statement[index - 1]
    if not getattr(prior, 'is_name', False):
        return False
    if str(prior).lower() in NON_CALLABLE_PAREN_NAMES:
        return False
    if index >= 2 and statement[index - 2] == '%':
        return False
    if str(prior).lower() in context.variable_names:
        return False
    return is_known_callable(str(prior), context)


def apply_spacing_policy(statement, context: Context) -> None:
    """Apply callable comma spacing policy to one statement."""
    paren_stack: list[bool] = []
    for index, token in enumerate(statement.tokens):
        if token == '(':
            paren_stack.append(opens_callable_group(statement, index, context))
        elif token in ('[', '{', '(/'):
            paren_stack.append(False)
        elif token in (')', ']', '}', '/)'):
            if paren_stack:
                paren_stack.pop()
        elif token == ',' and paren_stack and paren_stack[-1]:
            set_comma_spacing(token)


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
    contexts: dict[int, Context],
    source_lines=None,
) -> tuple[str, list[Change]]:
    """Format parsed statements and return rendered source plus changes."""
    changes: list[Change] = []

    for statement in statements:
        if not getattr(statement, 'source_visible', True):
            continue

        original_statement = statement.source_text()
        if not original_statement:
            continue
        if statement.kind == 'declaration':
            continue

        context = contexts.get(id(statement), Context(set(), set()))
        original_body = statement.source_line()

        apply_spacing_policy(statement, context)

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


def report_callable_comma_spacing(
    paths,
    write=False,
    report=False,
    diff=False,
    excludes=(),
) -> int:
    """Report or fix callable comma spacing in parsed flint sources."""
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
            statement_contexts(source),
            source_lines,
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
    """Parse command-line arguments and run the callable comma checker."""
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
        '--exclude',
        action='append',
        default=[],
        help='Directory to exclude; may be repeated',
    )

    args = parser.parse_args()
    return report_callable_comma_spacing(
        args.paths,
        write=args.write,
        report=args.report,
        diff=args.diff,
        excludes=args.exclude,
    )


if __name__ == '__main__':
    raise SystemExit(main())
