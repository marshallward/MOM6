#!/usr/bin/env python
import argparse
import collections
import json
import sys

record_type = collections.defaultdict(lambda: float)
for rec in ('grain', 'pemin', 'pemax',):
    record_type[rec] = int


def main():
    desc = 'Parse MOM6 model stdout and return clock data in JSON format.'

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('--format', '-f', action='store_true')
    parser.add_argument('input')
    args = parser.parse_args()

    with open(args.input) as log:
        clocks = parse_clocks(log)

    if args.format:
        print(json.dumps(clocks, indent=4))
    else:
        print(json.dumps(clocks))


def parse_clocks(log):
    clock_start_msg = 'Tabulating mpp_clock statistics across'
    clock_end_msg = 'MPP_STACK high water mark='

    fields = []
    for line in log:
        if line.startswith(clock_start_msg):
            npes = line.lstrip(clock_start_msg).split()[0]

            # Get records
            fields = []
            line = next(log)

            # Skip blank lines
            while line.isspace():
                line = next(log)

            fields = line.split()

            # Exit this loop, begin clock parsing
            break

    clocks = {}
    for line in log:
        # Treat MPP_STACK usage as end of clock reports
        if line.lstrip().startswith(clock_end_msg):
            break

        record = line.split()[-len(fields):]

        clk = line.split(record[0])[0].strip()
        clocks[clk] = {}

        for fld, rec in zip(fields, record):
            rtype = record_type[fld]
            clocks[clk][fld] = rtype(rec)

    return clocks


if __name__ == '__main__':
    main()
