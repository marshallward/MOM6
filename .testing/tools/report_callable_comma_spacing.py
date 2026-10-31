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
from pathlib import Path

from flint.formatting.runner import Change
from flint.formatting.callables import (
    Context, NON_CALLABLE_PAREN_NAMES, collect_contexts, collect_spacing_edits,
    comma_spacing_edits, is_known_callable, opens_callable_group, prepare_check,
    statement_contexts, format_statements, run as report_callable_comma_spacing,
)


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
