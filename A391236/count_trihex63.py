#!/usr/bin/env python3
"""Count order-6 trihexagonal figure logs; inspect triples and complement pairs.

Usage:
  python count_trihex63.py /path/to/logs
  python count_trihex63.py /path/to/logs --triple 63 32 1
  python count_trihex63.py /path/to/logs --strict
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

N = 63
MAGIC = 96
FILE_TEMPLATE = 'trihex63_a47_{x}_all_a48_{lo}_{hi}.txt'
HEADER = re.compile(r'^a47\s*=\s*(\d+)\s*,\s*a48\s*=\s*(\d+)\s*,\s*a37\s*=\s*(\d+)\s*$')
SOLUTION = re.compile(r'^\s*(\d+)(?:st|nd|rd|th):\s*solution:\s*(.*)$')
TOTAL = re.compile(r'^total:\s*(\d+)\s*$')


def bounds(x: int) -> tuple[int, int] | None:
    """All y with x>y>z>=1, x+y+z=96 and x,y,z <=63."""
    lo = max(1, (MAGIC - x) // 2 + 1, MAGIC - x - N)
    hi = min(x - 1, MAGIC - x - 1, N)
    return (lo, hi) if lo <= hi else None


def triples():
    for x in range(1, N + 1):
        b = bounds(x)
        if b:
            lo, hi = b
            for y in range(lo, hi + 1):
                yield (x, y, MAGIC - x - y)


def complement(t: tuple[int, int, int]) -> tuple[int, int, int]:
    """Complement values 64-v, then reorder descending."""
    x, y, z = t
    return (N + 1 - z, N + 1 - y, N + 1 - x)


def read_log(path: Path, x: int, expected: set[tuple[int, int, int]]):
    counts = Counter()
    errors = []
    seen_headers = set()
    current = None
    total_declared = None
    solutions_seen = 0
    for line_no, raw in enumerate(path.read_text(encoding='utf-8-sig').splitlines(), 1):
        line = raw.strip()
        h = HEADER.fullmatch(line)
        if h:
            t = tuple(map(int, h.groups()))
            if t not in expected:
                errors.append(f'{path.name}:{line_no}: unexpected triple {t}')
            if t in seen_headers:
                errors.append(f'{path.name}:{line_no}: repeated triple header {t}')
            seen_headers.add(t)
            current = t
            continue
        s = SOLUTION.fullmatch(line)
        if s:
            solutions_seen += 1
            if current is None:
                errors.append(f'{path.name}:{line_no}: solution before triple header')
                continue
            numbers = [int(v) for v in s.group(2).split()]
            if len(numbers) != N or set(numbers) != set(range(1, N + 1)):
                errors.append(f'{path.name}:{line_no}: invalid 1..63 permutation')
            elif (numbers[25], numbers[26], numbers[32]) != current:
                errors.append(f'{path.name}:{line_no}: solution triple does not match header')
            counts[current] += 1
            continue
        m = TOTAL.fullmatch(line)
        if m:
            total_declared = int(m.group(1))
    if total_declared is None:
        errors.append(f'{path.name}: missing final total')
    elif total_declared != solutions_seen:
        errors.append(f'{path.name}: total={total_declared}, actual solutions={solutions_seen}')
    missing_headers = expected - seen_headers
    for t in sorted(missing_headers):
        errors.append(f'{path.name}: missing triple header {t}')
    return counts, errors, len(seen_headers), total_declared


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def analyze(folder: Path, output: Path, query=None, strict=False):
    output.mkdir(parents=True, exist_ok=True)
    all_triples = list(triples())
    valid = set(all_triples)
    counts = Counter()
    by_x = Counter()
    missing_files = []
    errors = []
    files = []
    for x in range(1, N + 1):
        b = bounds(x)
        if b is None:
            continue
        lo, hi = b
        filename = FILE_TEMPLATE.format(x=x, lo=lo, hi=hi)
        path = folder / filename
        expected = {(x, y, MAGIC-x-y) for y in range(lo, hi+1)}
        if not path.is_file():
            missing_files.append(filename)
            continue
        local, local_errors, nheaders, declared = read_log(path, x, expected)
        counts.update(local)
        by_x[x] = sum(local.values())
        errors.extend(local_errors)
        files.append(dict(a47=x, a48_min=lo, a48_max=hi, combinations=len(expected),
                          headers=nheaders, solutions=by_x[x], declared_total=declared, filename=filename))

    # Do not silently treat unsearched triples as proven zero.
    covered_x = {r['a47'] for r in files}
    covered = {t for t in all_triples if t[0] in covered_x}
    triple_rows = [dict(a47=x, a48=y, a37=z,
                        count=counts[(x,y,z)] if (x,y,z) in covered else '',
                        status='read' if (x,y,z) in covered else 'missing_file')
                   for x,y,z in all_triples]
    x_rows = [dict(a47=x, a48_min=bounds(x)[0], a48_max=bounds(x)[1],
                   combinations=bounds(x)[1]-bounds(x)[0]+1,
                   count=by_x[x] if x in covered_x else '',
                   status='read' if x in covered_x else 'missing_file')
              for x in range(1,N+1) if bounds(x)]
    pair_rows = []
    unequal = []
    unknown = []
    for t in all_triples:
        c = complement(t)
        if t > c:  # each unordered pair once
            continue
        available = t in covered and c in covered
        a = counts[t] if t in covered else None
        b = counts[c] if c in covered else None
        match = 'MATCH' if available and a == b else ('MISMATCH' if available else 'UNKNOWN')
        row = dict(a47=t[0], a48=t[1], a37=t[2], count=a if a is not None else '',
                   complement_a47=c[0], complement_a48=c[1], complement_a37=c[2],
                   complement_count=b if b is not None else '', status=match)
        pair_rows.append(row)
        if match == 'MISMATCH': unequal.append(row)
        if match == 'UNKNOWN': unknown.append(row)

    write_csv(output/'counts_by_a47.csv', list(x_rows[0]), x_rows)
    write_csv(output/'counts_by_triple.csv', list(triple_rows[0]), triple_rows)
    write_csv(output/'complement_pairs.csv', list(pair_rows[0]), pair_rows)
    write_csv(output/'files_summary.csv', list(files[0]) if files else
              ['a47','a48_min','a48_max','combinations','headers','solutions','declared_total','filename'], files)
    raw_total = sum(counts.values())
    summary = dict(expected_files=len(x_rows), read_files=len(files),
                   missing_files=missing_files, expected_triples=len(all_triples),
                   covered_triples=len(covered), observed_solutions=raw_total,
                   multiplied_by_8=raw_total*8, complement_pairs=len(pair_rows),
                   matched_pairs=sum(r['status']=='MATCH' for r in pair_rows),
                   mismatched_pairs=unequal, unknown_pairs=len(unknown), errors=errors,
                   complete=not missing_files and not errors)
    (output/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'Files: {len(files)}/{len(x_rows)}, triples: {len(covered)}/{len(all_triples)}')
    print(f'Solutions in read files: {raw_total:,}; multiplied by 8: {raw_total*8:,}')
    print(f'Complement pairs: {len(pair_rows)}; match={summary["matched_pairs"]}, mismatch={len(unequal)}, unknown={len(unknown)}')
    print(f'Missing files: {len(missing_files)}; log errors: {len(errors)}')
    for s in missing_files: print('  MISSING', s)
    for s in errors[:30]: print('  ERROR', s)
    for r in unequal[:30]: print('  MISMATCH', r)
    if query is not None:
        t = tuple(query)
        if t not in valid:
            print('Query: invalid triple (require x>y>z>=1, sum=96, values<=63)')
        else:
            c = complement(t)
            def describe(k):
                return str(counts[k]) if k in covered else 'UNKNOWN (missing file)'
            print(f'Query {t}: {describe(t)} solutions')
            print(f'Complement {c}: {describe(c)} solutions')
    print('Saved CSV/JSON in:', output.resolve())
    if strict and (missing_files or errors or unequal):
        raise SystemExit(1)
    return summary


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('folder', nargs='?', type=Path, default=Path('.'), help='log directory')
    p.add_argument('--output', type=Path, help='output directory (default: folder/order6analysis)')
    p.add_argument('--triple', nargs=3, type=int, metavar=('A47','A48','A37'))
    p.add_argument('--strict', action='store_true', help='exit nonzero if incomplete or mismatch')
    a = p.parse_args()
    analyze(a.folder, a.output or Path('order6analysis'), a.triple, a.strict)
