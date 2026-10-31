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
from pathlib import Path
from flint.formatting.runner import Change
from flint.formatting.operators import (
    SKIPPED_OPERATOR_ROLES, collect_spacing_edits, format_statements, prepare_check,
    run as report_operator_spacing,
)


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
