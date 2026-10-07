#!/usr/bin/env python3
"""Report real declarations without bracketed units in trailing comments.

The default mode is a dry run over explicit source paths. A real variable is
considered documented when the comment associated with that variable contains a
bracketed unit expression, such as ``[m]`` or ``[L T-1 ~> m s-1]``.

Options:
  ``--exclude`` skips matching paths when parsing directory trees.
  ``--report`` prints each variable missing units with its file and line number.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import NamedTuple

import flint


UNIT_RE = re.compile(r'\[[^\]\n]+\]')


class Issue(NamedTuple):
    """Description of one real variable documentation issue."""

    path: Path
    line_number: int
    name: str
    original: str
    comments: tuple[str, ...]
    message: str = 'missing bracketed units'


def is_real_declaration(statement) -> bool:
    """Return True if a statement declares real variables."""
    return (
        getattr(statement, 'kind', None) == 'declaration'
        and len(statement) > statement.code_index
        and statement[statement.code_index] == 'real'
    )


def token_line_number(statement, token) -> int:
    """Return the physical source line for a token in a statement."""
    pattern = re.compile(rf'\b{re.escape(str(token))}\b', re.IGNORECASE)
    for offset, line in enumerate(statement.source_line().splitlines()):
        code = line.split('!', 1)[0]
        if pattern.search(code):
            return statement.line_number + offset
    return statement.line_number + ''.join(token.head).count('\n')


def comments_in_liminals(liminals: list[str]) -> list[str]:
    """Return comments stored in a token head or tail."""
    return [item for item in liminals if item.startswith('!')]


def comments_starting_with(comments: list[str], prefix: str) -> list[str]:
    """Return a comment block starting with a given prefix."""
    for index, comment in enumerate(comments):
        if comment.startswith(prefix):
            return comments[index:]
    return []


def leading_forward_comment(statement) -> str:
    """Return a Doxygen block that documents the following statement."""
    return ' '.join(leading_forward_comments(statement))


def leading_forward_comments(statement) -> list[str]:
    """Return Doxygen comments that document the following statement."""
    return comments_starting_with(comments_in_liminals(statement[0].head), '!>')


def comment_after_token(statement, index: int) -> str:
    """Return the comment most directly associated with a token."""
    return ' '.join(comments_after_token(statement, index))


def comments_after_token(statement, index: int) -> list[str]:
    """Return ordinary comments most directly associated with a token."""
    depth = 0
    comments: list[str] = []
    for offset, token in enumerate(statement[index:], start=index):
        if token in ('(', '[', '{', '(/'):
            depth += 1
        elif token in (')', ']', '}', '/)') and depth > 0:
            depth -= 1

        comments.extend(comments_in_liminals(token.tail))

        if offset > index and token == ',' and depth == 0:
            break

        if comments:
            break

    return comments


def doxygen_tail_comments(statement, index: int) -> list[str]:
    """Return a backward Doxygen comment block for a declaration token."""
    comments = comments_starting_with(comments_in_liminals(statement[index].tail), '!<')
    if comments:
        return comments

    return comments_starting_with(statement_trailing_comments(statement), '!<')


def statement_trailing_comment(statement) -> str:
    """Return trailing comments on a declaration statement."""
    return ' '.join(statement_trailing_comments(statement))


def statement_trailing_comments(statement) -> list[str]:
    """Return trailing comments on a declaration statement."""
    comments: list[str] = []
    for token in statement:
        comments.extend(comments_in_liminals(token.tail))
    return comments


def associated_comments(statement, index: int) -> list[str]:
    """Return the best available comment block for a declaration token."""
    comments = doxygen_tail_comments(statement, index)
    if comments:
        return comments

    comments = comments_after_token(statement, index)
    if comments:
        return comments

    comments = statement_trailing_comments(statement)
    if comments:
        return comments

    return leading_forward_comments(statement)


def comment_blocks(lines: list[str], statement, index: int) -> list[list[str]]:
    """Return possible documentation blocks for a declaration token."""
    doxygen_tail = doxygen_tail_comments(statement, index)
    direct = comments_after_token(statement, index)
    leading = leading_forward_comments(statement)

    if doxygen_tail or direct:
        blocks = [doxygen_tail, direct, leading]
    else:
        blocks = [
            statement_trailing_comments(statement),
            following_comment_block(lines, statement),
            leading,
        ]

    result: list[list[str]] = []
    for block in blocks:
        if block and block not in result:
            result.append(block)
    return result


def has_competing_doxygen_comments(statement, index: int) -> bool:
    """Return True when both forward and backward Doxygen forms are present."""
    return bool(leading_forward_comments(statement) and doxygen_tail_comments(statement, index))


def declarator_name_indices(statement) -> list[int]:
    """Return token indices for variable names after a declaration separator."""
    try:
        index = statement.index('::') + 1
    except ValueError:
        return []

    indices: list[int] = []
    depth = 0
    expect_name = True
    while index < len(statement):
        token = statement[index]
        if token in ('(', '[', '{', '(/'):
            depth += 1
        elif token in (')', ']', '}', '/)') and depth > 0:
            depth -= 1
        elif token == ',' and depth == 0:
            expect_name = True
        elif expect_name and getattr(token, 'is_name', False):
            indices.append(index)
            expect_name = False
        index += 1

    return indices


def has_units(comment: str) -> bool:
    """Return True if a comment contains a bracketed unit expression."""
    return bool(UNIT_RE.search(comment))


def following_comment_block(lines: list[str], statement) -> list[str]:
    """Return an indented plain comment block following a declaration."""
    statement_lines = statement.source_line().splitlines()
    if not statement_lines:
        return []

    declaration_index = statement.line_number - 1
    last_code_offset = 0
    for offset, line in enumerate(statement_lines):
        if line.split('!', 1)[0].strip():
            last_code_offset = offset

    next_index = declaration_index + last_code_offset + 1
    if next_index >= len(lines):
        return []

    first = lines[next_index]
    stripped = first.lstrip(' ')
    if not stripped.startswith('!'):
        return []
    comments: list[str] = []
    while next_index < len(lines):
        line = lines[next_index]
        stripped = line.lstrip(' ')
        if not stripped.startswith('!'):
            break
        comments.append(stripped.rstrip())
        next_index += 1

    return comments


def source_line(path: Path, line_number: int) -> str:
    """Return one source line from a path."""
    lines = path.read_text(encoding='utf-8', errors='ignore').splitlines()
    if 1 <= line_number <= len(lines):
        return lines[line_number - 1]
    return ''


def declaration_text(line: str) -> str:
    """Return declaration source without an inline comment."""
    return line.split('!', 1)[0].rstrip()


def collect_units(unit):
    """Yield a unit and its contained subprograms and derived types."""
    yield unit
    for derived_type in getattr(unit, 'derived_types', []):
        yield from collect_units(derived_type)
    for subprogram in getattr(unit, 'subprograms', []):
        yield from collect_units(subprogram)


def real_unit_issues(source) -> list[Issue]:
    """Return real variables missing bracketed units in one flint source."""
    path = Path(source.path)
    lines = path.read_text(encoding='utf-8', errors='ignore').splitlines()
    issues: list[Issue] = []

    seen: set[tuple[int, str]] = set()

    for statement in source.statements:
        if not getattr(statement, 'source_visible', True):
            continue
        if not is_real_declaration(statement):
            continue

        for index in declarator_name_indices(statement):
            token = statement[index]
            seen.add((id(statement), str(token).lower()))

            blocks = comment_blocks(lines, statement, index)
            comments = blocks[0] if blocks else []

            competing_doxygen = has_competing_doxygen_comments(statement, index)
            if any(has_units(' '.join(block)) for block in blocks) and not competing_doxygen:
                continue

            line_number = token_line_number(statement, token)
            issues.append(
                Issue(
                    path=path,
                    line_number=line_number,
                    name=str(token),
                    original=source_line(path, line_number),
                    comments=tuple(comments),
                    message=(
                        'competing Doxygen comments'
                        if competing_doxygen
                        else 'missing bracketed units'
                    ),
                )
            )

    for unit in source.units:
        for subunit in collect_units(unit):
            for variable in subunit.variables:
                if variable.type != 'real':
                    continue

                statement = variable.stmt
                if (id(statement), str(variable.name).lower()) in seen:
                    continue
                if '::' in statement:
                    continue
                try:
                    index = statement.index(variable.name)
                except ValueError:
                    continue

                blocks = comment_blocks(lines, statement, index)
                comments = blocks[0] if blocks else []

                competing_doxygen = has_competing_doxygen_comments(statement, index)
                if any(has_units(' '.join(block)) for block in blocks) and not competing_doxygen:
                    continue

                line_number = token_line_number(statement, variable.name)
                issues.append(
                    Issue(
                        path=path,
                        line_number=line_number,
                        name=str(variable.name),
                        original=source_line(path, line_number),
                        comments=tuple(comments),
                        message=(
                            'competing Doxygen comments'
                            if competing_doxygen
                            else 'missing bracketed units'
                        ),
                    )
                )

    return issues


def report_real_units(paths, report=False, excludes=()) -> int:
    """Report real variables missing bracketed units."""
    all_issues: list[Issue] = []
    project = flint.parse(*(str(path) for path in paths), excludes=excludes)

    for source in project.sources:
        all_issues.extend(real_unit_issues(source))

    if report:
        for issue in all_issues:
            print(
                f'{issue.path}:{issue.line_number}: '
                f'{issue.name}: {issue.message}'
            )
            print(f'  declaration: {declaration_text(issue.original)}')
            if issue.comments:
                print('  comments:')
                for comment in issue.comments:
                    print(f'    {comment}')
            else:
                print('  comments: <no associated comment>')
    else:
        seen = set()
        for issue in all_issues:
            if issue.path in seen:
                continue
            seen.add(issue.path)
            print(issue.path)

    return 1 if all_issues else 0


def main() -> int:
    """Parse command-line arguments and run the real units checker."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='+', type=Path, help='Files to process')
    parser.add_argument(
        '--report',
        action='store_true',
        help='Report variables missing units with line numbers',
    )
    parser.add_argument(
        '--exclude',
        action='append',
        default=[],
        help='Directory to exclude; may be repeated',
    )

    args = parser.parse_args()
    return report_real_units(
        args.paths,
        report=args.report,
        excludes=args.exclude,
    )


if __name__ == '__main__':
    raise SystemExit(main())
