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
from flint.token import TokenKind


BINARY_OPERATOR_KINDS = {
    TokenKind.POINTER_ASSIGNMENT,
    TokenKind.ARITHMETIC_OPERATOR,
    TokenKind.RELATIONAL_OPERATOR,
    TokenKind.LOGICAL_OPERATOR,
    TokenKind.DEFINED_OPERATOR,
}

GROUP_STARTERS = {'(', '[', '{', '(/'}
GROUP_ENDERS = {')': '(', ']': '[', '}': '{', '/)': '(/'}
NON_CALLABLE_GROUP_NAMES = {
    'associate', 'do', 'else', 'elseif', 'forall', 'if', 'select', 'where',
    'while',
}

UNARY_SIGN_CONTEXT = {
    '(', '[', '{', '(/', ',', ':', '::', '=', '=>', '+', '-', '*', '/', '//',
    '<', '>', '<=', '>=', '==', '/=', '.and.', '.or.', '.eqv.', '.neqv.',
    '.not.',
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


def is_unary_sign(statement, index: int) -> bool:
    """Return True if a plus or minus token is acting as a unary sign."""
    token = statement[index]
    if token not in ('+', '-'):
        return False
    if index <= statement.code_index:
        return True
    prior = statement[index - 1]
    return prior in UNARY_SIGN_CONTEXT or prior.is_operator


def in_declaration_prefix(statement, index: int, declaration_end: int | None) -> bool:
    """Return True for tokens before the declarator section of declarations."""
    if statement.kind != 'declaration':
        return False
    if declaration_end is None:
        return True
    return index < declaration_end


def opens_callable_or_subscript_group(statement, index: int) -> bool:
    """Return True if a delimiter opens a likely argument/subscript group."""
    token = statement[index]
    if token in ('[', '{', '(/'):
        return True
    if token != '(' or index == 0:
        return False

    prior = statement[index - 1]
    return (
        getattr(prior, 'is_name', False)
        and str(prior).lower() not in NON_CALLABLE_GROUP_NAMES
    )


def in_do_control(statement, index: int) -> bool:
    """Return True for tokens in a do-loop control clause."""
    return statement.is_do_statement() and index > statement.code_index


def is_print_format_star(statement, index: int) -> bool:
    """Return True for the star in print *, output statements."""
    return (
        index > statement.code_index
        and statement[statement.code_index] == 'print'
        and statement[index] == '*'
        and getattr(statement[index - 1], 'is_name', False)
        and statement[index - 1] == 'print'
        and len(statement) > index + 1
        and statement[index + 1] == ','
    )


def apply_spacing_policy(statement, strict=False) -> None:
    """Apply operator spacing policy to one statement."""
    group_stack: list[tuple[str, bool, bool]] = []
    callable_or_subscript_depth = 0
    operator_spec_depth = 0
    try:
        declaration_end = statement.index('::')
    except ValueError:
        declaration_end = None

    for index, token in enumerate(statement.tokens):
        if token in GROUP_ENDERS and group_stack:
            _, in_callable_or_subscript, in_operator_spec = group_stack.pop()
            if in_callable_or_subscript:
                callable_or_subscript_depth -= 1
            if in_operator_spec:
                operator_spec_depth -= 1

        if token.kind not in BINARY_OPERATOR_KINDS:
            if token in GROUP_STARTERS:
                in_callable_or_subscript = opens_callable_or_subscript_group(
                    statement, index,
                )
                in_operator_spec = (
                    token == '('
                    and index > 0
                    and statement[index - 1] == 'operator'
                )
                group_stack.append((
                    str(token), in_callable_or_subscript, in_operator_spec,
                ))
                if in_callable_or_subscript:
                    callable_or_subscript_depth += 1
                if in_operator_spec:
                    operator_spec_depth += 1
            continue

        if (
            token == '**'
            or operator_spec_depth > 0
            or in_declaration_prefix(statement, index, declaration_end)
            or (
                token.kind == TokenKind.ARITHMETIC_OPERATOR
                and (
                    callable_or_subscript_depth > 0
                    or in_do_control(statement, index)
                    or is_print_format_star(statement, index)
                )
            )
        ):
            continue
        if token == '.not.':
            set_spacing(token, None, ' ') if strict else set_minimum_spacing(
                token, None, ' ',
            )
        elif is_unary_sign(statement, index):
            set_spacing(token, None, '')
        elif strict:
            set_spacing(token, ' ', ' ')
        else:
            set_minimum_spacing(token, ' ', ' ')

        if token in GROUP_STARTERS:
            in_callable_or_subscript = opens_callable_or_subscript_group(
                statement, index,
            )
            in_operator_spec = (
                token == '('
                and index > 0
                and statement[index - 1] == 'operator'
            )
            group_stack.append((
                str(token), in_callable_or_subscript, in_operator_spec,
            ))
            if in_callable_or_subscript:
                callable_or_subscript_depth += 1
            if in_operator_spec:
                operator_spec_depth += 1


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
