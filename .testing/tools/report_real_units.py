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
    """Description of one real variable missing bracketed units."""

    path: Path
    line_number: int
    name: str
    original: str
    comments: tuple[str, ...]


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


def associated_liminal_comments(liminals: list[str]) -> list[str]:
    """Return comments that look associated with preceding source code."""
    comments: list[str] = []
    newline_before_comment = False
    for item in liminals:
        if item.startswith('!'):
            if (
                not comments
                and newline_before_comment
                and not item.startswith(('!<', '!>', '!!'))
            ):
                return []
            comments.append(item)
        elif '\n' in item and not comments:
            newline_before_comment = True
    return comments


def leading_forward_comment(statement) -> str:
    """Return a Doxygen block that documents the following statement."""
    return ' '.join(leading_forward_comments(statement))


def leading_forward_comments(statement) -> list[str]:
    """Return Doxygen comments that document the following statement."""
    comments = comments_in_liminals(statement[0].head)
    for index, comment in enumerate(comments):
        if comment.startswith('!>'):
            return comments[index:]
    return []


def comment_after_token(statement, index: int) -> str:
    """Return the comment most directly associated with a token."""
    return ' '.join(comments_after_token(statement, index))


def comments_after_token(statement, index: int) -> list[str]:
    """Return comments most directly associated with a token."""
    depth = 0
    comments: list[str] = []
    for offset, token in enumerate(statement[index:], start=index):
        if token in ('(', '[', '{', '(/'):
            depth += 1
        elif token in (')', ']', '}', '/)') and depth > 0:
            depth -= 1

        comments.extend(associated_liminal_comments(token.tail))

        if offset > index and token == ',' and depth == 0:
            break

        if comments:
            break

    return comments


def statement_trailing_comment(statement) -> str:
    """Return trailing comments on a declaration statement."""
    return ' '.join(statement_trailing_comments(statement))


def statement_trailing_comments(statement) -> list[str]:
    """Return trailing comments on a declaration statement."""
    comments: list[str] = []
    for token in statement:
        comments.extend(associated_liminal_comments(token.tail))
    return comments


def associated_comments(statement, index: int) -> list[str]:
    """Return the best available comment block for a declaration token."""
    comments = comments_after_token(statement, index)
    if comments:
        return comments

    comments = statement_trailing_comments(statement)
    if comments:
        return comments

    return leading_forward_comments(statement)


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

            comments = associated_comments(statement, index)
            comment = ' '.join(comments)

            if has_units(comment):
                continue

            line_number = token_line_number(statement, token)
            issues.append(
                Issue(
                    path=path,
                    line_number=line_number,
                    name=str(token),
                    original=source_line(path, line_number),
                    comments=tuple(comments),
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

                comments = associated_comments(statement, index)
                comment = ' '.join(comments)

                if has_units(comment):
                    continue

                line_number = token_line_number(statement, variable.name)
                issues.append(
                    Issue(
                        path=path,
                        line_number=line_number,
                        name=str(variable.name),
                        original=source_line(path, line_number),
                        comments=tuple(comments),
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
                f'{issue.name}: missing bracketed units'
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
