#!/usr/bin/env python3
"""Report Fortran callable comma spacing inconsistencies.

The default mode is a dry run over explicit source paths. Commas separating
known function or subroutine arguments require at least one following space.
Extra alignment spaces and spacing before a continuation marker are tolerated.
Commas in known array references, such as ``A(i,j,k)``, are left unchanged.

Options:
  ``--report`` prints each changed source line with its file and line number.
  ``--diff`` prints unified diffs for proposed changes.
  ``--write`` applies changes in place.
  ``--exclude`` skips matching paths when parsing directory trees.
  ``--strict`` requires exactly one space, including before continuation ``&``.
  ``--exclude-function`` skips comma checks for a named callable; repeatable.
"""

from __future__ import annotations

import argparse
import difflib
from pathlib import Path
from typing import NamedTuple

import flint
from flint.intrinsics import intrinsic_fns
from flint.formatting.spacing import format_spacing, propose_edit, spacing_edits


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
    """A line diagnostic with its source path attached for reporting."""

    path: Path
    line_number: int
    original: str


class Context(NamedTuple):
    """Known data objects and callables in scope for a statement."""

    variable_names: set[str]
    callables: set[str]


def comma_spacing_edits(token, strict=False):
    """Propose minimum spacing, or exact spacing when strict is requested.

    Non-strict checks leave continuation gaps untouched, including ``,&``.
    Strict checks normalize simple spacing before a continuation marker but
    preserve comments and the subsequent newline and indentation.
    """
    newline_indices = [
        index for index, item in enumerate(token.tail) if '\n' in item
    ]
    if newline_indices:
        if not strict:
            return []
        newline_index = newline_indices[0]
        prefix = token.tail[:newline_index]
        if '&' in prefix and all(
            item.isspace() or item == '&' for item in prefix
        ):
            return propose_edit(token.tail, [' ', '&'] + token.tail[newline_index:])
        return []
    return spacing_edits(token, after=' ', minimum=not strict)


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


def collect_spacing_edits(
    statement, context: Context, strict=False, excluded_functions=(),
):
    """Return callable comma proposals without changing parsed tokens.

    Exclusions match callable names case-insensitively and only skip that
    call's own commas. Nested calls are checked independently.
    """
    if statement.kind == 'declaration':
        return []
    edits = []
    excluded = {name.lower() for name in excluded_functions}
    paren_stack: list[bool] = []
    for index, token in enumerate(statement.tokens):
        if token == '(':
            check_commas = opens_callable_group(statement, index, context)
            if check_commas and str(statement[index - 1]).lower() in excluded:
                check_commas = False
            paren_stack.append(check_commas)
        elif token in ('[', '{', '(/'):
            paren_stack.append(False)
        elif token in (')', ']', '}', '/)'):
            if paren_stack:
                paren_stack.pop()
        elif token == ',' and paren_stack and paren_stack[-1]:
            edits.extend(comma_spacing_edits(token, strict=strict))
    return edits


def format_statements(
    path: Path,
    statements,
    contexts: dict[int, Context],
    source_lines=None,
    strict=False,
    excluded_functions=(),
) -> tuple[str, list[Change]]:
    """Format parsed statements and return rendered source plus changes."""
    formatted, changes = format_spacing(
        statements,
        lambda statement: collect_spacing_edits(
            statement, contexts.get(id(statement), Context(set(), set())),
            strict=strict, excluded_functions=excluded_functions,
        ),
        source_lines,
    )
    return formatted, [
        Change(path, change.line_number, change.original) for change in changes
    ]


def report_callable_comma_spacing(
    paths,
    write=False,
    report=False,
    diff=False,
    excludes=(),
    strict=False,
    excluded_functions=(),
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
            strict=strict,
            excluded_functions=excluded_functions,
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
        '--strict',
        action='store_true',
        help='Require exactly one space after commas, including before continuation &',
    )
    parser.add_argument(
        '--exclude-function',
        action='append',
        default=[],
        metavar='NAME',
        help='Skip comma checks for this callable name (case-insensitive); may be repeated',
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
        strict=args.strict,
        excluded_functions=args.exclude_function,
    )


if __name__ == '__main__':
    raise SystemExit(main())
