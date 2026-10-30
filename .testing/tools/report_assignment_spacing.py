#!/usr/bin/env python3
"""Report Fortran assignment spacing inconsistencies.

The default mode is a dry run over explicit source paths. Assignment operators
require one or more spaces on each side. Name-value pairs, such as named
arguments and declaration specifiers, require no spaces around the equals sign.
Explicit and implied-do control equals signs are left unchanged by default.

Options:
  ``--report`` prints each changed source line with its file and line number.
  ``--diff`` prints unified diffs for proposed changes.
  ``--write`` applies changes in place.
  ``--exclude`` skips matching paths when parsing directory trees.
  ``--strict`` requires exactly one space around assignment operators.
  ``--do-controls`` also checks explicit and implied-do control equals signs.
"""

from __future__ import annotations

import argparse
from functools import partial
from pathlib import Path
from flint.formatting.runner import Change, run_spacing
from flint.formatting.spacing import format_spacing, spacing_edits


def prepare_check(source, do_controls=False, strict=False):
    """Prepare this rule for a shared Flint spacing run."""
    return partial(collect_spacing_edits, do_controls=do_controls, strict=strict)


def collect_spacing_edits(statement, do_controls=False, strict=False):
    """Return assignment spacing proposals without changing parsed tokens."""
    edits = []
    for token in statement.tokens:
        if token.syntax_role == 'assignment':
            edits.extend(spacing_edits(token, ' ', ' ', minimum=not strict))
        elif token.syntax_role == 'name_value':
            edits.extend(spacing_edits(token, '', ''))
        elif token.syntax_role == 'do_control' and do_controls:
            edits.extend(spacing_edits(token, '', ''))

    if do_controls and statement.is_do_statement():
        for token in statement:
            if str(token) == ',':
                edits.extend(spacing_edits(token, after=''))
    return edits


def format_statements(
    path: Path,
    statements,
    source_lines=None,
    do_controls=False,
    strict=False,
) -> tuple[str, list[Change]]:
    """Format parsed statements and return rendered source plus changes."""
    formatted, changes = format_spacing(
        statements,
        lambda statement: collect_spacing_edits(
            statement, do_controls=do_controls, strict=strict,
        ),
        source_lines,
    )
    return formatted, [
        Change(path, change.line_number, change.original) for change in changes
    ]


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
    return run_spacing(
        paths, checks=(partial(prepare_check, do_controls=do_controls, strict=strict),),
        excludes=excludes, write=write, report=report, diff=diff,
    )


def main() -> int:
    """Parse command-line arguments and run the assignment spacing checker."""
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        'paths',
        nargs='+',
        type=Path,
        help='Files to process',
    )
    parser.add_argument(
        '--write',
        action='store_true',
        help='Write changes in place',
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
