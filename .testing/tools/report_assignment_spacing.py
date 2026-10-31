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
from pathlib import Path
from flint.formatting.runner import Change
from flint.formatting.assignment import (
    collect_spacing_edits, format_statements, prepare_check,
    run as report_assignment_spacing,
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
