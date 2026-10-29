#!/usr/bin/env python3
"""Report Fortran indentation inconsistencies using Flint statements.

Block bodies use two spaces. Module contents are not indented, matching MOM6
style. Continuation code is at least four spaces beyond its statement's first
line; the action of a continued single-line IF may instead use two spaces.
Declaration continuations allow at least two extra spaces by default;
``--strict`` requires at least four for declarations too.
Blank lines, standalone comments (including Doxygen documentation), preprocessor
directives, and inactive preprocessor branches are not checked. Comments between
continuation lines do not end the continued statement.

The default dry run lists affected paths. ``--report`` prints line diagnostics;
``--exclude`` skips matching paths when parsing directory trees. No files are
modified. This initial checker covers structured free-form Fortran, not legacy
label-terminated DO loops.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import NamedTuple

import flint


class Issue(NamedTuple):
    """An indentation diagnostic attached to a physical source line."""

    path: Path
    line_number: int
    original: str
    message: str


class Block(NamedTuple):
    """Expected header and body indentation for an open construct."""

    kind: str
    indent: int
    body: int


def code_line_indices(statement):
    """Map physical code lines to their first token, excluding liminal comments.

    Flint keeps whitespace, comments, and continuation markers in token tails.
    Count their newlines for locations, but only tokens introduce code lines.
    A token itself can span lines, for example a continued character literal.
    """
    line_number = statement.line_number
    lines = {}
    for index, token in enumerate(statement):
        text = str(token)
        if text:
            lines.setdefault(line_number, index)
            for _ in range(text.count('\n')):
                line_number += 1
                lines.setdefault(line_number, index)
        line_number += ''.join(token.tail).count('\n')
    return lines


def action_index(statement):
    """Return the action token of a single-line IF, or None for block IFs."""
    start = statement.code_index
    if statement[start:start + 2] != ['if', '(']:
        return None
    depth = 0
    for index in range(start + 1, len(statement)):
        token = statement[index]
        if token == '(':
            depth += 1
        elif token == ')':
            depth -= 1
            if depth == 0:
                following = index + 1
                if following < len(statement) and statement[following] != 'then':
                    return following
                return None
    return None


def block_event(statement):
    """Classify structural statements without treating keyword variables as blocks."""
    words = [str(token).lower() for token in statement.code_tokens]
    if not words:
        return None, None
    first = words[0]
    if len(words) > 1 and words[1] in ('=', '=>', '%'):
        return None, None
    if first == 'end':
        return 'close', words[1] if len(words) > 1 else None
    if first in {
        'endif', 'enddo', 'endwhere', 'endforall', 'endselect', 'endassociate',
        'endblock', 'endcritical', 'endtype', 'endinterface', 'endenum',
        'endsubroutine', 'endfunction', 'endprogram', 'endmodule', 'endsubmodule',
    }:
        return 'close', first[3:]
    if first in ('else', 'elseif', 'elsewhere'):
        return 'branch', 'where' if first == 'elsewhere' else 'if'
    if first == 'case' or words[:2] in (['type', 'is'], ['class', 'is'], ['class', 'default']):
        return 'case', 'select'
    if first == 'contains':
        return 'contains', None
    if first == 'if' and words[-1] == 'then':
        return 'open', 'if'
    if first in ('where', 'forall'):
        depth = 0
        for index, word in enumerate(words[1:], 1):
            if word == '(':
                depth += 1
            elif word == ')':
                depth -= 1
                if depth == 0:
                    return ('open', first) if index == len(words) - 1 else (None, None)
        return None, None
    if first in ('do', 'select', 'associate', 'critical', 'enum', 'interface'):
        return 'open', first
    if words[:2] == ['abstract', 'interface']:
        return 'open', 'interface'
    if first == 'block' and words[1:2] != ['data']:
        return 'open', 'block'
    if first == 'type' and words[1:2] != ['(']:
        return 'open', 'type'
    # Prefixes such as PURE, RECURSIVE, or REAL(...) may precede a procedure.
    prefixes = {
        'subroutine', 'function', 'pure', 'impure', 'elemental', 'recursive',
        'non_recursive', 'module', 'real', 'integer', 'logical', 'complex',
        'character', 'double', 'type', 'class',
    }
    if first in prefixes:
        depth = 0
        for index, word in enumerate(words):
            if word == '(':
                depth += 1
            elif word == ')':
                depth -= 1
            elif depth == 0 and word in ('::', '='):
                break
            elif depth == 0 and word in ('subroutine', 'function'):
                if index + 1 < len(words):
                    return 'open', word
    if first == 'module' and words[1:2] not in (['procedure'], ['subroutine'], ['function']):
        return 'open', 'module'
    if first in ('program', 'submodule'):
        return 'open', first
    return None, None


def indentation_issues(source, strict=False) -> list[Issue]:
    """Check block indentation and physical continuations in a Flint source."""
    path = Path(source.path)
    lines = path.read_text(encoding='utf-8', errors='ignore').splitlines()
    blocks: list[Block] = []
    issues: list[Issue] = []
    checked = set()
    line_indents = {}

    for statement in source.statements:
        if not statement or not getattr(statement, 'source_visible', True):
            continue
        event, kind = block_event(statement)
        expected = blocks[-1].body if blocks else 0
        if event in ('close', 'branch', 'case'):
            # Match by construct, so a case body is not confused with SELECT.
            match = next(
                (index for index in range(len(blocks) - 1, -1, -1)
                 if kind is None or blocks[index].kind == kind),
                None,
            )
            if match is not None:
                block = blocks[match]
                expected = block.body if event == 'case' else block.indent
                keep = match if event == 'close' else match + 1
                del blocks[keep:]
        elif event == 'contains' and blocks:
            block = blocks[-1]
            expected = block.body if block.kind in ('subroutine', 'function') else block.indent

        physical = code_line_indices(statement)
        first_line = statement.line_number
        # MOM6 often puts paired loop headers on one line. These share a
        # visual indentation level even though they open two logical blocks.
        if event == 'open' and first_line in line_indents:
            expected = line_indents[first_line]
        line_indents.setdefault(first_line, expected)
        original = lines[first_line - 1]
        base = len(original) - len(original.lstrip(' '))
        action = action_index(statement)
        declaration = (
            statement.kind == 'declaration'
            and event is None
            and statement.code_tokens[1:2] not in (['='], ['=>'], ['%'])
        )
        continuation_step = 2 if declaration and not strict else 4
        for line_number, token_index in physical.items():
            if line_number in checked:
                continue  # Multiple semicolon-separated statements share a line.
            checked.add(line_number)
            original = lines[line_number - 1]
            prefix = original[:len(original) - len(original.lstrip())]
            indent = len(prefix)
            if '\t' in prefix:
                message = 'indentation must use spaces, not tabs'
            elif line_number == first_line:
                if indent == expected:
                    continue
                message = f'expected {expected} spaces, found {indent}'
            else:
                if indent >= base + continuation_step or (token_index == action and indent == base + 2):
                    continue
                expectation = f'at least {base + continuation_step} spaces'
                if token_index == action:
                    expectation += f' (or {base + 2} for the IF action)'
                message = f'continuation: expected {expectation}, found {indent}'
            issues.append(Issue(path, line_number, original, message))

        if event == 'open':
            step = 0 if kind in ('module', 'submodule') else 2
            blocks.append(Block(kind, expected, expected + step))
        elif event == 'case':
            blocks.append(Block('case', expected, expected + 2))

    return issues


def report_indentation(paths, report=False, excludes=(), strict=False) -> int:
    """Report indentation issues without modifying source files."""
    project = flint.parse(*(str(path) for path in paths), excludes=excludes)
    found = False
    for source in project.sources:
        issues = indentation_issues(source, strict=strict)
        if not issues:
            continue
        found = True
        if report:
            for issue in issues:
                print(f'{issue.path}:{issue.line_number}: {issue.message}')
                print(f'  {issue.original}')
        else:
            print(source.path)
    return int(found)


def main() -> int:
    """Parse command-line arguments and run the indentation checker."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='+', type=Path, help='Files to process')
    parser.add_argument('--report', action='store_true', help='Report issues with line numbers')
    parser.add_argument('--exclude', action='append', default=[], help='Directory to exclude; repeatable')
    parser.add_argument(
        '--strict', action='store_true',
        help='Require at least four extra spaces for declaration continuations',
    )
    args = parser.parse_args()
    return report_indentation(
        args.paths, report=args.report, excludes=args.exclude, strict=args.strict,
    )


if __name__ == '__main__':
    raise SystemExit(main())
