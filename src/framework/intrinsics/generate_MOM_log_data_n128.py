#!/usr/bin/env python3
"""Generate lookup tables for MOM_log_data_n128.F90.

This script documents and regenerates a log_repro lookup table using the
table-reduction strategy described in MOM_log_data_n128.F90:

  * cover the conventional normalized-significand interval [1, 2),
  * divide one positive-normal binary64 exponent stride into 128 equal bins in
    the ordered floating-point bit grid,
  * search representable reciprocals near each nominal bin center's reciprocal,
  * define c = 1/invc and choose candidates whose -log(invc) nearly lies on a
    2**-43 grid, keeping the tabulated log(c) on that grid,
  * derive a high/low split of c for the non-FMA argument reduction.

The default reciprocal-first strategy keeps the reciprocal and split center
consistent while retaining grid-aligned log(c).  It screens candidate grid
crossings in extended precision and selects the best high-precision result
among the finalists; the search is not exhaustive over all reciprocal values.
Use --strategy log-first for the original method: choose log(c) on the grid,
then derive c = exp(log(c)) and round 1/c.  The --search-radius and --objective
options explore nearby grid points for that method.
Reciprocal searches require NumPy with extended-precision longdouble for
screening; final candidates are evaluated with Decimal at --precision digits.

The generated tables are printed to stdout as a complete Fortran module.
"""

from __future__ import annotations

import argparse
import struct
from decimal import Decimal, ROUND_HALF_EVEN, getcontext


# Number of lookup-table subdivisions.
LOG_NDIV = 128
# Bit-pattern step between adjacent binary64 exponents.
EXP_STRIDE = 2**52
# Spacing of the grid used to choose log(c).
LOGC_GRID = Decimal(2) ** -43


def float_to_bits(x: float) -> int:
    """Return the binary64 bit pattern of x as an unsigned integer."""

    return struct.unpack(">Q", struct.pack(">d", x))[0]


def bits_to_float(xb: int) -> float:
    """Return the binary64 value represented by xb."""

    return struct.unpack(">d", struct.pack(">Q", xb))[0]


def nearest_integer(x: Decimal) -> int:
    """Return the nearest integer to x, with ties to even."""

    return int(x.to_integral_value(rounding=ROUND_HALF_EVEN))


def bin_midpoint(j: int) -> Decimal:
    """Return the exact midpoint of table bin j."""

    start_bits = float_to_bits(1.0)
    step_bits = EXP_STRIDE // LOG_NDIV
    center_bits = start_bits + j * step_bits + step_bits // 2
    return Decimal.from_float(bits_to_float(center_bits))


def split_high_low(c: Decimal) -> tuple[float, float]:
    """Split c into a binary64 high part and residual low part."""

    c_hi = float(c)
    c_lo = float(c - Decimal.from_float(c_hi))
    return c_hi, c_lo


def choose_logc(c0: Decimal, search_radius: int, objective: str) -> Decimal:
    """Choose log(c) on the 2**-43 grid near the nominal bin center."""

    logc0 = getcontext().ln(c0)
    q0 = nearest_integer(logc0 / LOGC_GRID)

    best_q = q0
    best_score = None
    for q in range(q0 - search_radius, q0 + search_radius + 1):
        logc = Decimal(q) * LOGC_GRID
        c = getcontext().exp(logc)
        c_hi, c_lo = split_high_low(c)

        if objective == "nearest-log":
            score = abs(logc - logc0)
        elif objective == "nearest-c":
            score = abs(c - c0)
        elif objective == "nearest-c-hi":
            score = abs(Decimal.from_float(c_hi) - c0)
        elif objective == "smallest-c-lo":
            score = abs(Decimal.from_float(c_lo))
        else:
            raise ValueError(f"unknown objective {objective!r}")

        if best_score is None or score < best_score:
            best_score = score
            best_q = q

    return Decimal(best_q) * LOGC_GRID


def choose_reciprocal(c0: Decimal, radius: int):
    """Screen log-grid crossings, then refine with Decimal."""

    import numpy as np

    if np.finfo(np.longdouble).nmant < 63:
        raise ValueError(
            "reciprocal searches need extended-precision NumPy longdouble"
        )
    center_bits = float_to_bits(float(1 / c0))
    first_bits = center_bits - radius
    last_bits = center_bits + radius
    qlo = nearest_integer(
        -getcontext().ln(Decimal.from_float(bits_to_float(last_bits)))
        / LOGC_GRID
    )
    qhi = nearest_integer(
        -getcontext().ln(Decimal.from_float(bits_to_float(first_bits)))
        / LOGC_GRID
    )

    # Near a log-grid crossing, only the closest binary64 reciprocals can win.
    # Screen crossings instead of examining every reciprocal bit pattern.
    grid = (
        np.arange(qlo, qhi + 1, dtype=np.int64).astype(np.longdouble)
        * np.longdouble(2) ** -43
    )
    target = np.exp(-grid)
    rounded = target.astype(np.float64)
    scores = np.abs(target - rounded.astype(np.longdouble)) / target
    nbest = min(256, len(scores))
    shortlist = np.argpartition(scores, nbest - 1)[:nbest]

    candidates = {first_bits, center_bits, last_bits}
    for index in shortlist:
        bits = float_to_bits(float(rounded[index]))
        candidates.update(
            b
            for b in (bits - 1, bits, bits + 1)
            if first_bits <= b <= last_bits
        )

    best = None
    best_key = None
    for bits in sorted(candidates):
        v_float = bits_to_float(bits)
        v = Decimal.from_float(v_float)
        exact_logc = -getcontext().ln(v)
        logc = Decimal(nearest_integer(exact_logc / LOGC_GRID)) * LOGC_GRID
        key = (abs(exact_logc - logc), abs(bits - center_bits))
        if best_key is None or key < best_key:
            best_key = key
            best = (v_float, logc, 1 / v)
    return best


def generate_tables(
    search_radius: int,
    objective: str,
    strategy="reciprocal-first",
    reciprocal_radius=2**29,
):
    """Generate log table arrays."""

    invc = []
    logc = []
    c_hi = []
    c_lo = []

    for j in range(LOG_NDIV):
        c0 = bin_midpoint(j)
        if strategy == "reciprocal-first":
            invcj, logcj, cj = choose_reciprocal(c0, reciprocal_radius)
        else:
            logcj = choose_logc(c0, search_radius, objective)
            cj = getcontext().exp(logcj)
            invcj = float(1 / cj)
        c_hi_j, c_lo_j = split_high_low(cj)

        logc.append(float(logcj))
        invc.append(invcj)
        c_hi.append(c_hi_j)
        c_lo.append(c_lo_j)

    return invc, logc, c_hi, c_lo


def format_array(name: str, values: list[float], comment: str) -> str:
    """Return a Fortran parameter array definition."""

    lines = [comment]
    lines.append(f"real, parameter :: {name}(0:ndiv-1) = [ &")
    for i in range(0, len(values), 4):
        chunk = values[i : i + 4]
        suffix = ", &" if i + 4 < len(values) else " &"
        body = ", ".join(f"{v: .17e}" for v in chunk)
        lines.append(f"    {body}{suffix}")
    lines.append("]")
    return "\n".join(lines)


def print_module(
    invc,
    logc,
    c_hi,
    c_lo,
    objective: str,
    search_radius: int,
    strategy="reciprocal-first",
    reciprocal_radius=2**29,
):
    """Print a complete Fortran table module."""

    print("! This file is part of MOM6, the Modular Ocean Model version 6.")
    print("! See the LICENSE file for licensing information.")
    print("! SPDX-License-Identifier: Apache-2.0")
    print()
    print(
        "!> Lookup tables for the reproducible log() function, "
        "assuming real64 precision."
    )
    print("!! Generated by generate_MOM_log_data_n128.py.")
    if strategy == "reciprocal-first":
        print(
            f"!! Reciprocal-first search; radius: {reciprocal_radius} "
            "binary64 steps on each side."
        )
    else:
        print(
            f"!! log(c) search objective: {objective}; "
            f"search radius: {search_radius} grid steps."
        )
    print("module MOM_log_data_n128")
    print()
    print("implicit none ; private")
    print()
    print("public :: ndiv")
    print("public :: inv_c, log_c, c_hi, c_lo")
    print()
    print("!> Number of subdivisions in the log table interval.")
    print(f"integer, parameter :: ndiv = {LOG_NDIV}")
    print()
    print(
        format_array(
            "inv_c",
            invc,
            "!> Lookup table of 1 / c for log_repro() "
            "table interval centers [nondim]",
        )
    )
    print()
    print(
        format_array(
            "log_c",
            logc,
            "!> Lookup table of log(c) for log_repro() "
            "table interval centers [nondim]",
        )
    )
    print()
    print(
        format_array(
            "c_hi",
            c_hi,
            "!> Lookup table of high parts of c for log_repro() "
            "table interval centers [nondim]",
        )
    )
    print()
    print(
        format_array(
            "c_lo",
            c_lo,
            "!> Lookup table of low parts of c for log_repro() "
            "table interval centers [nondim]",
        )
    )
    print()
    print("end module MOM_log_data_n128")


def main() -> None:
    """Parse command-line options and print the generated Fortran module."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--precision", type=int, default=100, help="decimal working precision"
    )
    parser.add_argument(
        "--strategy",
        choices=["log-first", "reciprocal-first"],
        default="reciprocal-first",
        help="choose grid-aligned logarithms first, "
        "or search representable reciprocals",
    )
    parser.add_argument(
        "--reciprocal-radius",
        type=int,
        default=2**29,
        help="binary64 reciprocal steps to search on each side "
        "for reciprocal-first",
    )
    parser.add_argument(
        "--search-radius",
        type=int,
        default=0,
        help="number of 2**-43 log(c) grid points to search on each side",
    )
    parser.add_argument(
        "--objective",
        default="nearest-log",
        choices=["nearest-log", "nearest-c", "nearest-c-hi", "smallest-c-lo"],
        help="criterion for selecting among searched log(c) grid points",
    )
    args = parser.parse_args()

    if args.reciprocal_radius < 0 or args.search_radius < 0:
        parser.error("search radii must be nonnegative")
    if args.strategy == "reciprocal-first" and (
        args.search_radius != 0 or args.objective != "nearest-log"
    ):
        parser.error("--search-radius and --objective apply only to log-first")

    if EXP_STRIDE % LOG_NDIV != 0:
        raise SystemExit(
            "LOG_NDIV must evenly divide one binary64 exponent stride"
        )

    getcontext().prec = args.precision
    invc, logc, c_hi, c_lo = generate_tables(
        args.search_radius,
        args.objective,
        args.strategy,
        args.reciprocal_radius,
    )
    print_module(
        invc,
        logc,
        c_hi,
        c_lo,
        args.objective,
        args.search_radius,
        args.strategy,
        args.reciprocal_radius,
    )


if __name__ == "__main__":
    main()
