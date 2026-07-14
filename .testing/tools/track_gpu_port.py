#!/usr/bin/env python3
"""Cross-reference gcov execution coverage against GPU-ported source regions.

Ported-region detection is automatic for:
  - `do concurrent` ... `enddo`/`end do` blocks
  - `!$omp target` / `!$omp target teams` block constructs (paired with their
    `end target[ teams]`)
  - `!$omp target teams loop`, `!$omp target teams distribute parallel do`,
    and bare `!$omp loop` combined constructs (attached to the following
    do-loop, which has no separate closing directive)

`!$omp target data` / `end target data` brackets a data environment, not a
compute region, so its contents are NOT auto-marked ported (only nested
compute directives inside it are). `!$omp target enter/exit data` and
`!$omp target update` are single-statement data-movement directives, counted
separately from compute regions.

For code with no directive at all (e.g. pure/elemental helpers invoked from
device code), wrap it manually:

    ! @gpu-ported-start
    ...
    ! @gpu-ported-end

Executed-but-unported lines are further split into "portable" and
"not portable" so the completion metric isn't diluted by code that will
never move to the GPU:
  - portable: lines inside any `do` loop nest (concurrent or not — a plain
    do-loop is exactly the kind of code a port targets), or whole-array
    Fortran assignment syntax (`h(:,:,:) = 0.0`)
  - not portable: `allocate`/`deallocate`, I/O statements, block-if
    control-flow keywords (`if (...) then`, `else if (...) then`, `else`,
    `end if`), and `call` statements — always, even if they appear inside a
    loop — plus everything else that executes outside any loop nest —
    scalar setup, control-flow bookkeeping, subroutine-entry index
    computation, and the like. A single-line `if (cond) stmt` is NOT
    excluded — the statement half may be real work.

The structural heuristic above is necessarily approximate, so it can be
overridden by hand in either direction:

    !@start noport         a loop (or any region) that will never be ported
    ...                    — e.g. serial halo-exchange bookkeeping, a tiny
    !@end noport           setup loop not worth offloading. Removed from
                            "portable" even if it structurally looks like a
                            loop/array-assignment.

    !@start toport          non-loop code that genuinely needs porting work
    ...                     — e.g. a scalar reduction to restructure for the
    !@end toport            device. Added to "portable" even though it
                            doesn't match the structural heuristic.

If `!@start noport`/`!@start toport` is immediately followed by a do-loop
(nothing but comments/blank lines in between), it attaches to that loop and
needs no matching `!@end` at all — the loop's own `end do` closes it, same
as how `!$omp loop` attaches to the following loop with no closing
directive of its own:

    !@start noport
    do i = 1, n
      ...
    enddo

Anything else (not immediately followed by a loop) falls back to requiring
an explicit, balanced `!@end noport`/`!@end toport`. Writing an explicit end
after a marker that already auto-attached to a loop is an error (it's
orphaned — the loop already closed it).

`noport` always wins over both the structural heuristic and `toport` if
they overlap. Markers nest like `@gpu-ported-start`/`-end` and must be
balanced within a file.

Usage:
    track_gpu_port.py --src-root <dir> --gcov-dir <dir> [--out-md FILE] [--out-json FILE]
    track_gpu_port.py --src-root <dir> --validate-only

Both the markdown and JSON reports carry exact line RANGES (not just counts)
for "portable, not yet ported" and "already ported" code, per file and per
routine — this is the authoritative "which lines" answer.

  - `--github-annotations` emits GitHub Actions `::notice
    file=...,line=...,endLine=...::` workflow commands straight from the
    same line ranges — these render natively on a PR's "Files changed" tab
    with zero dependency on hosting an HTML artifact, and are the
    recommended CI integration.
"""

import argparse
import html
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

DO_OPEN_RE = re.compile(r'^\s*(?:\w+\s*:\s*)?do\b', re.IGNORECASE)
DO_CONCURRENT_RE = re.compile(r'^\s*(?:\w+\s*:\s*)?do\s+concurrent\b', re.IGNORECASE)
END_DO_RE = re.compile(r'^\s*end\s*do\b', re.IGNORECASE)

OMP_RE = re.compile(r'^\s*!\$omp\s+(.*)$', re.IGNORECASE)
TAG_START_RE = re.compile(r'^\s*!\s*@gpu-ported-start\s*$', re.IGNORECASE)
TAG_END_RE = re.compile(r'^\s*!\s*@gpu-ported-end\s*$', re.IGNORECASE)
NOPORT_START_RE = re.compile(r'^\s*!\s*@start\s+noport\s*$', re.IGNORECASE)
NOPORT_END_RE = re.compile(r'^\s*!\s*@end\s+noport\s*$', re.IGNORECASE)
TOPORT_START_RE = re.compile(r'^\s*!\s*@start\s+toport\s*$', re.IGNORECASE)
TOPORT_END_RE = re.compile(r'^\s*!\s*@end\s+toport\s*$', re.IGNORECASE)

SUBPROG_START_RE = re.compile(
    r'^\s*(?:recursive\s+)?(?:pure\s+)?(?:elemental\s+)?(?:module\s+)?'
    r'(subroutine|function)\s+(\w+)', re.IGNORECASE)
SUBPROG_END_RE = re.compile(r'^\s*end\s*(subroutine|function)\s*(\w+)?', re.IGNORECASE)

# "Not portable" overrides: memory management and I/O never belong on the
# device, even when they happen to sit inside a loop nest.
ALLOC_RE = re.compile(r'^\s*(allocate|deallocate)\s*\(', re.IGNORECASE)
IO_RE = re.compile(
    r'^\s*(write|read|open|close|inquire|rewind|backspace|endfile)\s*[\(\*]', re.IGNORECASE)
PRINT_RE = re.compile(r'^\s*print\b', re.IGNORECASE)

# Block-if control-flow lines are pure branching, not work — never a porting target on
# their own, even inside a loop. Single-line "if (cond) stmt" is NOT matched here since
# the statement half may be real work; only the block-if keywords themselves are excluded.
IF_THEN_RE = re.compile(r'^\s*(?:\w+\s*:\s*)?if\s*\(.*\)\s*then\s*$', re.IGNORECASE)
ELSEIF_RE = re.compile(r'^\s*else\s*if\s*\(.*\)\s*then\s*$', re.IGNORECASE)
ELSE_RE = re.compile(r'^\s*else\s*$', re.IGNORECASE)
ENDIF_RE = re.compile(r'^\s*end\s*if(?:\s+\w+)?\s*$', re.IGNORECASE)

# A "call" statement invokes a subroutine outside the current unit — not itself GPU
# work, even inside a loop (the callee, if it needs porting, is tracked separately
# where it's defined).
CALL_RE = re.compile(r'^\s*call\s+\w', re.IGNORECASE)

# Whole-array Fortran assignment, e.g. "h(:,:,:) = 0.0" — a portable
# construct even outside an explicit do-loop.
ARRAY_ASSIGN_RE = re.compile(r'^\s*[A-Za-z_]\w*\s*\([^=()]*:[^=()]*\)\s*=[^=]')


@dataclass
class FileRegions:
    ported: set = field(default_factory=set)
    directive: set = field(default_factory=set)
    manual: set = field(default_factory=set)
    loop: set = field(default_factory=set)          # any do-loop body, ported or not
    array_syntax: set = field(default_factory=set)  # whole-array assignment lines
    nonportable: set = field(default_factory=set)   # allocate/deallocate/IO override
    manual_noport: set = field(default_factory=set)  # "!@start noport" — force-excluded
    manual_toport: set = field(default_factory=set)  # "!@start toport" — force-included
    warnings: list = field(default_factory=list)

    @property
    def portable(self):
        """Structurally portable lines (loop bodies or array syntax), overrides excluded.

        Manual markers win over the structural heuristic: `toport` force-includes lines
        the heuristic would otherwise skip (e.g. scalar code); `noport` force-excludes
        lines the heuristic would otherwise include, and wins if both are present.
        """
        auto = (self.loop | self.array_syntax) - self.nonportable
        return (auto | self.manual_toport) - self.manual_noport


class StructuralError(Exception):
    pass


def strip_comment(text):
    """Remove an unquoted trailing '!' comment from a plain Fortran statement."""
    in_quote = None
    for i, ch in enumerate(text):
        if in_quote:
            if ch == in_quote:
                in_quote = None
        elif ch in ("'", '"'):
            in_quote = ch
        elif ch == '!':
            return text[:i]
    return text


def split_statements(text):
    """Split a physical line on unquoted ';' into individual statements."""
    parts = []
    buf = []
    in_quote = None
    for ch in text:
        if in_quote:
            buf.append(ch)
            if ch == in_quote:
                in_quote = None
        elif ch in ("'", '"'):
            in_quote = ch
            buf.append(ch)
        elif ch == ';':
            parts.append(''.join(buf))
            buf = []
        else:
            buf.append(ch)
    parts.append(''.join(buf))
    return parts


def iter_logical_lines(lines):
    """Join '&'-continued physical lines. Yields (start, end, joined_text), 1-indexed inclusive."""
    i = 0
    n = len(lines)
    while i < n:
        start = i
        raw = lines[i].rstrip('\n')
        text = raw
        stripped = raw.strip()
        is_omp = bool(OMP_RE.match(raw)) or stripped.lower().startswith('!$omp')
        while text.rstrip().endswith('&') and i + 1 < n:
            i += 1
            nxt = lines[i].rstrip('\n')
            nxt_body = nxt
            if is_omp:
                m = re.match(r'^\s*!\$omp\s?', nxt, re.IGNORECASE)
                if m:
                    nxt_body = nxt[m.end():]
            text = text.rstrip()[:-1] + ' ' + nxt_body.strip()
        yield start + 1, i + 1, text
        i += 1


def classify_omp(body):
    """Classify the body of an '!$omp ...' directive (text after the sentinel)."""
    norm = re.sub(r'\s+', ' ', body).strip().lower()
    if norm.startswith('end target teams'):
        return 'end_teams'
    if norm.startswith('end target data'):
        return 'end_data'
    if norm.startswith('end target'):
        return 'end_target'
    if norm.startswith('target data'):
        return 'open_data'
    if norm.startswith('target enter data') or norm.startswith('target exit data') \
            or norm.startswith('target update'):
        return 'directive_line'
    if norm.startswith('target teams') and re.search(r'\b(loop|distribute)\b', norm):
        return 'combined'
    if norm.startswith('target teams'):
        return 'open_teams'
    if norm.startswith('target parallel') and re.search(r'\b(loop|do)\b', norm):
        return 'combined'
    if norm.startswith('target') and re.search(r'\bloop\b', norm):
        return 'combined'
    if norm.startswith('target'):
        return 'open_target'
    if norm.startswith('loop'):
        return 'loop_only'
    return 'other'


def parse_ported_regions(path):
    """Return a FileRegions describing one source file's ported/portable/nonportable lines."""
    text = path.read_text(errors='replace')
    lines = text.split('\n')
    if lines and lines[-1] == '':
        lines = lines[:-1]

    fr = FileRegions()
    manual_regions = []
    noport_regions = []
    toport_regions = []

    do_stack = []       # entries: {'concurrent', 'start', 'parent_start', 'combined_start',
                         #           'noport_start', 'toport_start'}
    region_stack = []   # entries: (kind, start_line) for target/teams/data blocks
    tag_stack = []       # entries: start_line for manual tags
    noport_stack = []    # entries: start_line for "!@start noport" not yet closed/attached
    toport_stack = []    # entries: start_line for "!@start toport" not yet closed/attached
    pending_combined_start = None
    # If "!@start noport"/"!@start toport" is immediately followed by a do-loop, it attaches
    # to that loop (closed implicitly by "end do") and needs no matching "!@end" — mirrors how
    # "!$omp loop" attaches to the following loop. Otherwise it falls back to requiring an
    # explicit "!@end noport"/"!@end toport" (tracked via noport_stack/toport_stack above).
    pending_noport_start = None
    pending_toport_start = None
    # Every closed do-loop's (start, end, parent_start), used to exclude nested loops from a
    # noport/toport marking — the marker should cover the loop(s) it directly names, not loops
    # nested inside them (those keep their own independent classification).
    all_loop_records = []

    def mark(lo, hi, dest):
        for ln in range(lo, hi + 1):
            dest.add(ln)

    def shallow_span(lo, hi):
        """[lo, hi] minus lines belonging to any loop nested inside another loop that is
        itself within [lo, hi] — i.e. keep the loop(s) directly named by the marker, drop
        loops nested inside those."""
        excluded = set()
        for rec in all_loop_records:
            s, e, ps = rec['start'], rec['end'], rec['parent_start']
            if s < lo or e > hi:
                continue
            if ps is not None and ps >= lo:
                excluded.update(range(s, e + 1))
        return set(range(lo, hi + 1)) - excluded

    for start, end, text_joined in iter_logical_lines(lines):
        stripped = text_joined.strip()
        if not stripped or stripped.startswith('!') and not OMP_RE.match(stripped):
            if stripped.startswith('!') and not OMP_RE.match(stripped):
                # comment-only line; still check for manual tags
                if TAG_START_RE.match(stripped):
                    tag_stack.append(start)
                elif TAG_END_RE.match(stripped):
                    if not tag_stack:
                        raise StructuralError(
                            f'{path}:{start}: @gpu-ported-end with no matching -start')
                    tstart = tag_stack.pop()
                    mark(tstart, end, fr.manual)
                    manual_regions.append((tstart, end))
                elif NOPORT_START_RE.match(stripped):
                    noport_stack.append(start)
                    pending_noport_start = start
                elif NOPORT_END_RE.match(stripped):
                    if not noport_stack:
                        raise StructuralError(
                            f'{path}:{start}: "!@end noport" with no matching "!@start noport"')
                    nstart = noport_stack.pop()
                    span = shallow_span(nstart, end)
                    fr.manual_noport.update(span)
                    noport_regions.append((nstart, end, span))
                    if pending_noport_start == nstart:
                        pending_noport_start = None
                elif TOPORT_START_RE.match(stripped):
                    toport_stack.append(start)
                    pending_toport_start = start
                elif TOPORT_END_RE.match(stripped):
                    if not toport_stack:
                        raise StructuralError(
                            f'{path}:{start}: "!@end toport" with no matching "!@start toport"')
                    tpstart = toport_stack.pop()
                    span = shallow_span(tpstart, end)
                    fr.manual_toport.update(span)
                    toport_regions.append((tpstart, end, span))
                    if pending_toport_start == tpstart:
                        pending_toport_start = None
            continue

        omp_m = OMP_RE.match(text_joined)
        if omp_m:
            kind = classify_omp(omp_m.group(1))
            if kind == 'end_teams':
                if not region_stack or region_stack[-1][0] != 'teams':
                    raise StructuralError(f'{path}:{start}: "end target teams" without matching open')
                _, rstart = region_stack.pop()
                mark(rstart, end, fr.ported)
            elif kind == 'end_data':
                if not region_stack or region_stack[-1][0] != 'data':
                    raise StructuralError(f'{path}:{start}: "end target data" without matching open')
                region_stack.pop()  # data region contents are not auto-ported
            elif kind == 'end_target':
                if not region_stack or region_stack[-1][0] != 'target':
                    raise StructuralError(f'{path}:{start}: "end target" without matching open')
                _, rstart = region_stack.pop()
                mark(rstart, end, fr.ported)
            elif kind == 'open_data':
                region_stack.append(('data', start))
            elif kind == 'open_teams':
                region_stack.append(('teams', start))
            elif kind == 'open_target':
                region_stack.append(('target', start))
            elif kind == 'directive_line':
                mark(start, end, fr.directive)
            elif kind in ('combined', 'loop_only'):
                if kind == 'loop_only' and region_stack and region_stack[-1][0] in ('teams', 'target'):
                    # already inside an open ported block; nothing further to attach
                    pass
                else:
                    if pending_combined_start is not None:
                        fr.warnings.append(
                            f'{path}:{pending_combined_start}: combined directive not '
                            f'immediately followed by a do-loop')
                    pending_combined_start = start
            # 'other' (plain CPU omp directives): ignored
            continue

        # plain Fortran: may contain several ';'-separated statements
        for stmt in split_statements(strip_comment(text_joined)):
            s = stmt.strip()
            if not s:
                continue
            if END_DO_RE.match(s):
                if not do_stack:
                    fr.warnings.append(
                        f'{path}:{end}: "end do" with no matching "do" (unbalanced; skipping)')
                    continue
                entry = do_stack.pop()
                mark(entry['start'], end, fr.loop)  # any do-loop body is a portable candidate
                if entry['concurrent']:
                    mark(entry['start'], end, fr.ported)
                if entry['combined_start'] is not None:
                    mark(entry['combined_start'], end, fr.ported)
                all_loop_records.append({
                    'start': entry['start'], 'end': end, 'parent_start': entry['parent_start']})
                if entry['noport_start'] is not None:
                    span = shallow_span(entry['noport_start'], end)
                    fr.manual_noport.update(span)
                    noport_regions.append((entry['noport_start'], end, span))
                if entry['toport_start'] is not None:
                    span = shallow_span(entry['toport_start'], end)
                    fr.manual_toport.update(span)
                    toport_regions.append((entry['toport_start'], end, span))
            elif DO_OPEN_RE.match(s):
                do_stack.append({
                    'concurrent': bool(DO_CONCURRENT_RE.match(s)),
                    'start': start,
                    'parent_start': do_stack[-1]['start'] if do_stack else None,
                    'combined_start': pending_combined_start,
                    'noport_start': pending_noport_start,
                    'toport_start': pending_toport_start,
                })
                pending_combined_start = None
                if pending_noport_start is not None:
                    noport_stack.pop()  # attached to this loop; no "!@end noport" required
                    pending_noport_start = None
                if pending_toport_start is not None:
                    toport_stack.pop()  # attached to this loop; no "!@end toport" required
                    pending_toport_start = None
            else:
                # marker(s) weren't immediately followed by a loop; fall back to requiring
                # an explicit "!@end noport"/"!@end toport" (still open on the stack)
                pending_noport_start = None
                pending_toport_start = None
                if (ALLOC_RE.match(s) or IO_RE.match(s) or PRINT_RE.match(s)
                        or IF_THEN_RE.match(s) or ELSEIF_RE.match(s)
                        or ELSE_RE.match(s) or ENDIF_RE.match(s) or CALL_RE.match(s)):
                    mark(start, end, fr.nonportable)
                elif ARRAY_ASSIGN_RE.match(s):
                    mark(start, end, fr.array_syntax)

    if region_stack:
        kinds = ', '.join(f'{k}@{ln}' for k, ln in region_stack)
        raise StructuralError(f'{path}: unclosed omp target region(s): {kinds}')
    if tag_stack:
        raise StructuralError(f'{path}: unclosed @gpu-ported-start at line(s) {tag_stack}')
    if noport_stack:
        raise StructuralError(f'{path}: unclosed "!@start noport" at line(s) {noport_stack}')
    if toport_stack:
        raise StructuralError(f'{path}: unclosed "!@start toport" at line(s) {toport_stack}')
    if do_stack:
        fr.warnings.append(
            f'{path}: {len(do_stack)} unbalanced "do" opener(s) at EOF (parser limitation)')

    # Flag manual tags that are already fully covered by automatic detection.
    for tstart, tend in manual_regions:
        span = set(range(tstart + 1, tend))
        if span and span.issubset(fr.ported):
            fr.warnings.append(
                f'{path}:{tstart}-{tend}: @gpu-ported tag is redundant '
                f'(region already detected automatically)')

    # Flag "noport" regions that contradict already-ported code, and "toport"
    # regions that are redundant with the structural heuristic. Both checks use the
    # *shallow* span (nested loops excluded) — that's what actually got marked.
    auto_portable = (fr.loop | fr.array_syntax) - fr.nonportable
    for nstart, nend, span in noport_regions:
        check = span - {nstart, nend}  # exclude the marker comment line(s) themselves
        if check & fr.ported:
            fr.warnings.append(
                f'{path}:{nstart}-{nend}: "noport" region overlaps automatically-detected '
                f'ported code (contradiction)')
    for tpstart, tpend, span in toport_regions:
        check = span - {tpstart, tpend}
        if check and check.issubset(auto_portable):
            fr.warnings.append(
                f'{path}:{tpstart}-{tpend}: "toport" tag is redundant '
                f'(region already structurally portable)')

    return fr


def parse_routines(path):
    """Best-effort (name, start, end) for subroutines/functions, for the per-routine report."""
    text = path.read_text(errors='replace')
    lines = text.split('\n')
    stack = []
    routines = []
    for lineno, raw in enumerate(lines, start=1):
        s = raw.strip()
        if not s or s.startswith('!'):
            continue
        m = SUBPROG_END_RE.match(s)
        if m and stack:
            name, start = stack.pop()
            routines.append((name, start, lineno))
            continue
        m = SUBPROG_START_RE.match(s)
        if m:
            stack.append((m.group(2), lineno))
    return routines


def parse_gcov(gcov_path):
    """Return (source_path_str, {line_no: exec_count_or_None})."""
    source = None
    counts = {}
    line_re = re.compile(r'^\s*([0-9]+|#####|-)\*?:\s*([0-9]+):')
    with open(gcov_path, errors='replace') as f:
        for raw in f:
            m = line_re.match(raw)
            if not m:
                continue
            count_s, lineno_s = m.group(1), m.group(2)
            lineno = int(lineno_s)
            if lineno == 0:
                sm = re.search(r'Source:(.*)$', raw)
                if sm:
                    source = sm.group(1).strip()
                continue
            if count_s == '-':
                continue  # non-executable line
            count = 0 if count_s == '#####' else int(count_s)
            counts[lineno] = max(count, counts.get(lineno, 0))
    return source, counts


def resolve_source(gcov_source, src_root, index):
    """Map a gcov Source: path onto a real file under src_root."""
    p = Path(gcov_source)
    if p.is_absolute() and p.exists():
        try:
            return p.resolve()
        except OSError:
            return p
    candidates = index.get(p.name, [])
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        # prefer the one whose path shares the longest suffix with gcov_source
        parts = p.parts
        best, best_score = None, -1
        for c in candidates:
            cparts = c.parts
            score = 0
            for a, b in zip(reversed(parts), reversed(cparts)):
                if a == b:
                    score += 1
                else:
                    break
            if score > best_score:
                best, best_score = c, score
        return best
    return None


def ranges_from_lines(line_set):
    """Collapse a set of line numbers into sorted, contiguous [start, end] ranges."""
    if not line_set:
        return []
    lines = sorted(line_set)
    ranges = []
    start = prev = lines[0]
    for ln in lines[1:]:
        if ln == prev + 1:
            prev = ln
        else:
            ranges.append([start, prev])
            start = prev = ln
    ranges.append([start, prev])
    return ranges


def format_ranges(ranges):
    return ', '.join(str(lo) if lo == hi else f'{lo}-{hi}' for lo, hi in ranges)


HTML_STYLE = """
:root { color-scheme: light dark; }
body { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
       margin: 0; padding: 1rem 1.5rem; }
h1 { font-size: 1rem; font-weight: 600; }
a { color: #2b6cb0; }
p.back { margin-top: 0; }
p.legend span { display: inline-block; padding: 1px 8px; margin-right: 4px; border-radius: 3px; }
table.summary { border-collapse: collapse; width: 100%; }
table.summary th, table.summary td { padding: 4px 10px; border-bottom: 1px solid #8883; text-align: left; }
table.summary td.num, table.summary th.num { text-align: right; }
table.src { border-collapse: collapse; font-size: 12.5px; width: 100%; }
table.src td { padding: 0 6px; white-space: pre; vertical-align: top; }
table.src td.ln { color: #888; text-align: right; user-select: none; }
table.src td.cnt { color: #888; text-align: right; min-width: 2.5em; user-select: none; }
.ported { background: #d9f2d9; }
.portable { background: #ffe6a8; }
.executed { background: #eaeaea; }
.nothit { background: #fad2d2; }
tr.nodata { background: transparent; }
@media (prefers-color-scheme: dark) {
  body { background: #1e1e1e; color: #ddd; }
  .ported { background: #234d23; }
  .portable { background: #5c4a13; }
  .executed { background: #333333; }
  .nothit { background: #5c2323; }
}
"""

LEGEND_HTML = (
    '<p class="legend">'
    '<span class="ported">ported</span>'
    '<span class="portable">portable, not yet ported</span>'
    '<span class="executed">executed, not portable</span>'
    '<span class="nothit">executable, not hit by this run</span>'
    '</p>'
)


def html_slug(rel_path):
    """Turn a src-root-relative path into a flat, collision-safe .html filename."""
    return re.sub(r'[^A-Za-z0-9_.-]', '_', str(rel_path).replace('/', '__')) + '.html'


def render_file_html(rel_path, resolved, counts, fr):
    """Render one source file as an lcov-style line-by-line coverage/port-status page."""
    try:
        lines = resolved.read_text(errors='replace').splitlines()
    except OSError:
        lines = []
    ported_all = fr.ported | fr.manual
    portable = fr.portable
    rows = []
    for lineno, text in enumerate(lines, start=1):
        count = counts.get(lineno)
        if count is None:
            cls, marker = 'nodata', ''
        elif count == 0:
            cls, marker = 'nothit', '0'
        elif lineno in ported_all:
            cls, marker = 'ported', str(count)
        elif lineno in portable:
            cls, marker = 'portable', str(count)
        else:
            cls, marker = 'executed', str(count)
        rows.append(f'<tr class="{cls}"><td class="ln">{lineno}</td><td class="cnt">{marker}</td>'
                     f'<td class="src">{html.escape(text)}</td></tr>')
    title = html.escape(str(rel_path))
    return (f'<!doctype html><html><head><meta charset="utf-8"><title>{title}</title>'
            f'<style>{HTML_STYLE}</style></head><body>'
            f'<p class="back"><a href="index.html">&larr; back to index</a></p>'
            f'<h1>{title}</h1>{LEGEND_HTML}'
            f'<table class="src">{"".join(rows)}</table></body></html>')


def render_index_html(per_file, total_ported, total_portable, total_executed,
                       total_not_portable, overall_pct_portable, overall_pct_executed):
    rows = []
    for r in per_file:
        pct = f"{r['pct_of_portable']:.1f}%" if r['pct_of_portable'] is not None else '—'
        link = r.get('html_file')
        file_cell = f'<a href="{link}">{html.escape(r["file"])}</a>' if link else html.escape(r['file'])
        rows.append(f'<tr><td>{file_cell}</td><td class="num">{r["executed_lines"]}</td>'
                     f'<td class="num">{r["ported_lines"]}</td>'
                     f'<td class="num">{r["portable_remaining_lines"]}</td>'
                     f'<td class="num">{r["not_portable_lines"]}</td>'
                     f'<td class="num">{pct}</td></tr>')
    return (f'<!doctype html><html><head><meta charset="utf-8"><title>GPU Port Coverage Report</title>'
            f'<style>{HTML_STYLE}</style></head><body>'
            f'<h1>GPU Port Coverage Report</h1>'
            f'<p><strong>Of portable code: {total_ported} / {total_portable} executed lines ported '
            f'({overall_pct_portable:.1f}%)</strong><br>'
            f'({total_ported} / {total_executed} of all executed lines, {overall_pct_executed:.1f}% '
            f'&mdash; {total_not_portable} executed lines are allocate/I/O/scalar setup code with no '
            f'porting target and are excluded from the portable metric above)</p>'
            f'{LEGEND_HTML}'
            f'<table class="summary"><tr><th>File</th><th class="num">Executed</th>'
            f'<th class="num">Ported</th><th class="num">Portable remaining</th>'
            f'<th class="num">Not portable</th><th class="num">% of portable</th></tr>'
            f'{"".join(rows)}</table></body></html>')


def build_index(src_root):
    index = {}
    for p in Path(src_root).rglob('*.F90'):
        index.setdefault(p.name, []).append(p)
    for p in Path(src_root).rglob('*.f90'):
        index.setdefault(p.name, []).append(p)
    return index


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--src-root', required=True, help='Root of the MOM6 source tree to scan for ported regions')
    ap.add_argument('--gcov-dir', help='Directory containing .gcov files (recursively searched)')
    ap.add_argument('--out-md', help='Write markdown report here (default: stdout)')
    ap.add_argument('--out-json', help='Write JSON snapshot here')
    ap.add_argument('--validate-only', action='store_true',
                     help='Only run the source-tree parser and report structural errors/warnings')
    ap.add_argument('--github-annotations', action='store_true',
                     help='Emit GitHub Actions ::notice file=...,line=...,endLine=...:: workflow '
                          'commands for each executed, portable, not-yet-ported line range — shows '
                          'up natively on the PR "Files changed" diff, no HTML artifact needed')
    ap.add_argument('--baseline-json', help='A prior --out-json snapshot (e.g. from the PR base '
                     'branch) to diff the current totals against, reported as "since baseline" '
                     'in --out-comment')
    ap.add_argument('--changed-files', help='File with one src-root-relative source path per '
                     'line (e.g. from `git diff --name-only`) — scopes an additional "this PR\'s '
                     'files" subset of the totals in --out-comment')
    ap.add_argument('--out-comment', help='Write a short PR-comment-sized markdown summary here: '
                     'overall %% ported, delta vs --baseline-json, and the --changed-files subset '
                     'if given. Intended to be posted/updated on the PR itself; --out-md remains '
                     'the full per-file/per-routine report for the job summary.')
    ap.add_argument('--verbose', action='store_true',
                     help='List every routine with executed, portable, not-yet-ported code (and '
                          'its exact line ranges) in --out-md, instead of just the top 25 by '
                          'portable-remaining lines')
    ap.add_argument('--out-html', help='Write a self-contained, lcov-style line-by-line HTML '
                     'report (index.html plus one page per executed file, colored by ported / '
                     'portable-remaining / executed-not-portable / not-hit) into this directory '
                     '— e.g. for publishing to GitHub Pages')
    args = ap.parse_args()

    src_root = Path(args.src_root).resolve()
    all_warnings = []
    ported_by_file = {}   # resolved Path -> FileRegions
    errors = []

    src_files = sorted(set(src_root.rglob('*.F90')) | set(src_root.rglob('*.f90')))
    for f in src_files:
        if not f.is_file():
            # Broken symlink — typically an uninitialized git submodule
            # (e.g. src/equation_of_state/TEOS10 -> pkg/GSW-Fortran). Nothing
            # to scan; not a parse error.
            all_warnings.append(f'{f}: skipped (target does not exist, likely '
                                 f'an uninitialized submodule)')
            continue
        try:
            fr = parse_ported_regions(f)
        except StructuralError as e:
            errors.append(str(e))
            continue
        ported_by_file[f] = fr
        all_warnings.extend(fr.warnings)

    if errors:
        print('STRUCTURAL ERRORS (must fix before this report is trustworthy):', file=sys.stderr)
        for e in errors:
            print(f'  {e}', file=sys.stderr)

    if args.validate_only:
        for w in all_warnings:
            print(f'WARNING: {w}', file=sys.stderr)
        total_ported = sum(len(fr.ported) for fr in ported_by_file.values())
        total_directive = sum(len(fr.directive) for fr in ported_by_file.values())
        total_manual = sum(len(fr.manual) for fr in ported_by_file.values())
        total_portable = sum(len(fr.portable) for fr in ported_by_file.values())
        total_nonportable = sum(len(fr.nonportable) for fr in ported_by_file.values())
        total_noport = sum(len(fr.manual_noport) for fr in ported_by_file.values())
        total_toport = sum(len(fr.manual_toport) for fr in ported_by_file.values())
        print(f'Parsed {len(ported_by_file)} files: '
              f'{total_ported} auto-ported compute lines, '
              f'{total_directive} data-movement directive lines, '
              f'{total_manual} manually tagged lines, '
              f'{total_portable} portable-structure lines (loop bodies / array syntax), '
              f'{total_nonportable} allocate/IO override lines, '
              f'{total_noport} manually excluded (noport) lines, '
              f'{total_toport} manually included (toport) lines, '
              f'{len(all_warnings)} warnings, {len(errors)} structural errors.')
        sys.exit(1 if errors else 0)

    if not args.gcov_dir:
        ap.error('--gcov-dir is required unless --validate-only is given')

    index = build_index(src_root)
    gcov_files = sorted(Path(args.gcov_dir).rglob('*.gcov'))
    if not gcov_files:
        print(f'No .gcov files found under {args.gcov_dir}', file=sys.stderr)
        sys.exit(1)

    exec_by_file = {}  # resolved Path -> {line: count}
    unresolved = []
    for gf in gcov_files:
        source, counts = parse_gcov(gf)
        if source is None or not counts:
            continue
        resolved = resolve_source(source, src_root, index)
        if resolved is None:
            unresolved.append(source)
            continue
        merged = exec_by_file.setdefault(resolved, {})
        for ln, cnt in counts.items():
            merged[ln] = max(cnt, merged.get(ln, 0))

    def split_executed(executed, fr):
        """Partition a file's executed lines into (ported, portable-remaining, not-portable)."""
        ported_all = fr.ported | fr.manual
        executed_ported = executed & ported_all
        executed_portable_remaining = (executed & fr.portable) - executed_ported
        executed_not_portable = executed - executed_ported - executed_portable_remaining
        return executed_ported, executed_portable_remaining, executed_not_portable

    html_dir = None
    if args.out_html:
        html_dir = Path(args.out_html)
        html_dir.mkdir(parents=True, exist_ok=True)

    per_file = []
    total_executed = 0
    total_ported = 0
    total_portable_remaining = 0
    total_not_portable = 0
    routine_rows = []
    empty_fr = FileRegions()

    for resolved, counts in exec_by_file.items():
        executed = {ln for ln, c in counts.items() if c > 0}
        if not executed:
            continue
        fr = ported_by_file.get(resolved, empty_fr)
        executed_ported, executed_portable, executed_notportable = split_executed(executed, fr)
        try:
            rel = resolved.relative_to(src_root)
        except ValueError:
            rel = resolved
        portable_total = len(executed_ported) + len(executed_portable)
        pct_portable = 100.0 * len(executed_ported) / portable_total if portable_total else None
        per_file.append({
            'file': str(rel),
            'executed_lines': len(executed),
            'ported_lines': len(executed_ported),
            'portable_remaining_lines': len(executed_portable),
            'not_portable_lines': len(executed_notportable),
            'pct_of_portable': round(pct_portable, 1) if pct_portable is not None else None,
            'portable_remaining_ranges': ranges_from_lines(executed_portable),
            'ported_ranges': ranges_from_lines(executed_ported),
        })
        if html_dir is not None:
            html_name = html_slug(rel)
            (html_dir / html_name).write_text(render_file_html(rel, resolved, counts, fr))
            per_file[-1]['html_file'] = html_name
        total_executed += len(executed)
        total_ported += len(executed_ported)
        total_portable_remaining += len(executed_portable)
        total_not_portable += len(executed_notportable)

        for name, rstart, rend in parse_routines(resolved):
            rexec = {ln for ln in executed if rstart <= ln <= rend}
            if not rexec:
                continue
            rported = rexec & (fr.ported | fr.manual)
            rportable = (rexec & fr.portable) - rported
            rtotal = len(rported) + len(rportable)
            routine_rows.append({
                'file': str(rel),
                'routine': name,
                'executed_lines': len(rexec),
                'ported_lines': len(rported),
                'portable_remaining_lines': len(rportable),
                'pct_of_portable': round(100.0 * len(rported) / rtotal, 1) if rtotal else None,
                'portable_remaining_ranges': ranges_from_lines(rportable),
                'ported_ranges': ranges_from_lines(rported),
            })

    # Rank by the biggest *portable* gap, not raw unported lines — allocate/IO/scalar
    # setup code that will never be ported shouldn't crowd out real porting targets.
    per_file.sort(key=lambda r: r['portable_remaining_lines'], reverse=True)
    routine_rows.sort(key=lambda r: r['portable_remaining_lines'], reverse=True)
    total_portable = total_ported + total_portable_remaining
    overall_pct_portable = 100.0 * total_ported / total_portable if total_portable else 0.0
    overall_pct_executed = 100.0 * total_ported / total_executed if total_executed else 0.0

    baseline = None
    if args.baseline_json:
        try:
            baseline = json.loads(Path(args.baseline_json).read_text())
        except (OSError, json.JSONDecodeError) as e:
            all_warnings.append(
                f'--baseline-json {args.baseline_json}: could not load ({e}); skipping delta')

    delta = None
    if baseline is not None:
        b = baseline.get('overall', {})
        delta = {
            'ported_lines': total_ported - b.get('ported_lines', 0),
            'pct_of_portable': round(overall_pct_portable - b.get('pct_of_portable', 0.0), 1),
        }

    changed_files = None
    pr_rows = []
    pr_ported = pr_portable_total = 0
    pr_pct = None
    if args.changed_files:
        changed_files = {line.strip() for line in Path(args.changed_files).read_text().splitlines()
                          if line.strip()}
        pr_rows = [r for r in per_file if r['file'] in changed_files]
        pr_ported = sum(r['ported_lines'] for r in pr_rows)
        pr_portable_total = pr_ported + sum(r['portable_remaining_lines'] for r in pr_rows)
        pr_pct = 100.0 * pr_ported / pr_portable_total if pr_portable_total else None

    lines_out = []
    lines_out.append('# GPU Port Coverage Report\n')
    lines_out.append(f'**Of portable code: {total_ported} / {total_portable} executed lines '
                      f'ported ({overall_pct_portable:.1f}%)**  ')
    lines_out.append(f'({total_ported} / {total_executed} of *all* executed lines, '
                      f'{overall_pct_executed:.1f}% — {total_not_portable} executed lines are '
                      f'`allocate`/I/O/scalar setup code with no porting target and are excluded '
                      f'from the portable metric above)\n')
    if errors:
        lines_out.append(f'\n:warning: **{len(errors)} structural error(s) in ported-region '
                          f'annotations** — see job log. Report below excludes affected files.\n')
    lines_out.append('\n## By file (biggest portable gap first)\n')
    lines_out.append('| File | Executed | Ported | Portable remaining | Not portable | % of portable |')
    lines_out.append('|---|---:|---:|---:|---:|---:|')
    for r in per_file:
        pct = f"{r['pct_of_portable']:.1f}%" if r['pct_of_portable'] is not None else '—'
        lines_out.append(f"| {r['file']} | {r['executed_lines']} | {r['ported_lines']} | "
                          f"{r['portable_remaining_lines']} | {r['not_portable_lines']} | {pct} |")

    shown_routines = routine_rows if args.verbose else routine_rows[:25]
    heading = 'All routines with portable-remaining lines' if args.verbose \
        else 'Top porting opportunities (by executed portable-remaining lines)'
    lines_out.append(f'\n## {heading}\n')
    lines_out.append('| Routine | File | Executed | Ported | Portable remaining | % of portable |')
    lines_out.append('|---|---|---:|---:|---:|---:|')
    for r in shown_routines:
        pct = f"{r['pct_of_portable']:.1f}%" if r['pct_of_portable'] is not None else '—'
        lines_out.append(f"| {r['routine']} | {r['file']} | {r['executed_lines']} | "
                          f"{r['ported_lines']} | {r['portable_remaining_lines']} | {pct} |")

    range_heading = 'Line ranges for all porting opportunities' if args.verbose \
        else 'Line ranges for top porting opportunities'
    lines_out.append(f'\n## {range_heading}\n')
    lines_out.append('Exact executed, portable, not-yet-ported line ranges — this is the '
                      'authoritative "which lines" answer; any HTML/CI visualization is '
                      'derived from these same ranges, not the other way around.\n')
    for r in shown_routines:
        if not r['portable_remaining_ranges']:
            continue
        lines_out.append(f"- **{r['routine']}** (`{r['file']}`): lines "
                          f"{format_ranges(r['portable_remaining_ranges'])}")

    if unresolved:
        lines_out.append(f'\n_{len(unresolved)} .gcov file(s) could not be matched to a source '
                          f'file under {src_root} and were skipped._\n')
    if all_warnings:
        lines_out.append(f'\n_{len(all_warnings)} parser warning(s); see job log._\n')

    report = '\n'.join(lines_out)
    if args.out_md:
        Path(args.out_md).write_text(report + '\n')
    else:
        print(report)

    if args.out_json:
        Path(args.out_json).write_text(json.dumps({
            'overall': {
                'executed_lines': total_executed,
                'ported_lines': total_ported,
                'portable_remaining_lines': total_portable_remaining,
                'not_portable_lines': total_not_portable,
                'pct_of_portable': round(overall_pct_portable, 1),
                'pct_of_executed': round(overall_pct_executed, 1),
            },
            'by_file': per_file,
            'by_routine': routine_rows,
            'warnings': all_warnings,
            'errors': errors,
            'unresolved_gcov_sources': unresolved,
        }, indent=2))

    if html_dir is not None:
        (html_dir / 'index.html').write_text(render_index_html(
            per_file, total_ported, total_portable, total_executed,
            total_not_portable, overall_pct_portable, overall_pct_executed))

    if args.out_comment:
        c = []
        c.append('### GPU Port Coverage\n')
        c.append(f'**Overall: {total_ported} / {total_portable} portable executed lines ported '
                  f'({overall_pct_portable:.1f}%)**')
        if delta is not None:
            sign = '+' if delta['ported_lines'] >= 0 else ''
            psign = '+' if delta['pct_of_portable'] >= 0 else ''
            c.append(f"Since base branch: {sign}{delta['ported_lines']} ported lines "
                      f"({psign}{delta['pct_of_portable']:.1f} pp)")
        if changed_files is not None:
            c.append('')
            if pr_rows:
                pct_s = f'{pr_pct:.1f}%' if pr_pct is not None else '—'
                c.append(f'**Files touched by this PR: {pr_ported} / {pr_portable_total} '
                          f'portable executed lines ported ({pct_s})**')
            else:
                c.append('_No executed, GPU-portable lines in the files touched by this PR '
                          '(nothing in the diff was exercised by the coverage run, or none of '
                          "it is GPU-portable)._")
        c.append('')
        c.append('<sub>Full per-file / per-routine breakdown: see the "gpu-port-report" job '
                  'summary and artifact.</sub>')
        Path(args.out_comment).write_text('\n'.join(c) + '\n')

    if args.github_annotations:
        for r in per_file:
            for lo, hi in r['portable_remaining_ranges']:
                if lo == hi:
                    print(f"::notice file={r['file']},line={lo}::GPU-portable, not yet ported")
                else:
                    print(f"::notice file={r['file']},line={lo},endLine={hi}::"
                          f"GPU-portable, not yet ported ({hi - lo + 1} lines)")

    for w in all_warnings:
        print(f'WARNING: {w}', file=sys.stderr)

    if errors:
        sys.exit(1)


if __name__ == '__main__':
    main()
