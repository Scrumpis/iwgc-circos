#!/usr/bin/env python3
"""Validate external interval tables and convert them to BED coordinates.

No scientific estimation lives here: bedtools still calculates coverage/means.
"""
import argparse
import math
from pathlib import Path
import re


def fail(message):
    raise ValueError(message)


def index(path):
    lengths = {}
    for number, line in enumerate(Path(path).read_text().splitlines(), 1):
        f = line.split()
        if len(f) < 2 or not f[1].isdigit() or int(f[1]) <= 0:
            fail(f'{path}:{number}: invalid sequence length')
        if not re.fullmatch(r'[A-Za-z0-9_.:-]+', f[0]) or f[0] in lengths:
            fail(f'{path}:{number}: duplicate or unsupported sequence ID: {f[0]}')
        lengths[f[0]] = int(f[1])
    if not lengths:
        fail(f'{path}: no sequences')
    return lengths


def interval(chrom, start, end, lengths, where, closed=False):
    if not re.fullmatch(r'[0-9]+', start) or not re.fullmatch(r'[0-9]+', end):
        fail(f'{where}: coordinates must be integers')
    a, b = int(start), int(end)
    if a > b:
        a, b = b, a
    if closed:
        a -= 1
    if chrom not in lengths or a < 0 or b > lengths[chrom] or a >= b:
        fail(f'{where}: invalid/out-of-bounds interval {chrom}:{start}..{end}')
    return a, b


def convert(args):
    lengths = index(args.fai)
    retained = index(args.retained)
    if args.mode == 'index':
        return
    handles = {}
    schema = None
    with open(args.input) as handle:
        for number, line in enumerate(handle, 1):
            where = f'{args.input}:{number}'
            if line.rstrip() == '##FASTA' and args.mode in ('gene', 'repeat', 'intact'):
                break
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            if args.mode in ('gene', 'repeat', 'intact'):
                f = line.rstrip('\r\n').split('\t')
                if len(f) != 9:
                    fail(f'{where}: expected nine tab-separated GFF3 columns')
                if int(f[3]) > int(f[4]):
                    fail(f'{where}: GFF3 start must not exceed end')
                a, b = interval(f[0], f[3], f[4], lengths, where, closed=True)
                if f[0] not in retained or (args.mode == 'gene' and f[2] != 'gene'):
                    continue
                if args.mode == 'intact':
                    attrs = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
                    classification = attrs.get('Classification', '')
                    parts = classification.split('/')
                    if len(parts) < 2 or not parts[1]:
                        fail(f'{where}: intact feature requires Classification=type/class')
                    family = parts[1]
                    if not re.fullmatch(r'[A-Za-z0-9_-]+', family):
                        fail(f'{where}: unsafe intact class name')
                    if family not in handles:
                        handles[family] = open(Path(args.output) / (family + '.bed'), 'w')
                    print(f[0], a, b, sep='\t', file=handles[family])
                else:
                    print(f[0], a, b, sep='\t')
            elif args.mode == 'ltr':
                f = line.rstrip('\r\n').split('\t')
                m = re.fullmatch(r'(.+):([0-9]+)\.\.([0-9]+)', f[0])
                if len(f) < 12 or not m:
                    fail(f'{where}: expected chr:start..end and age in tab column 12')
                age = float(f[11])
                if not math.isfinite(age) or age < 0:
                    fail(f'{where}: LTR age must be finite and nonnegative')
                a, b = interval(m[1], m[2], m[3], lengths, where, closed=True)
                if m[1] in retained:
                    print(m[1], a, b, age, sep='\t')
            elif args.mode == 'links':
                f = line.split()
                current = len(f)
                if current not in (6, 11, 13) or (schema is not None and schema != current):
                    fail(f'{where}: expected consistent 6, 11 or 13-column links')
                schema = current
                if current == 6:
                    c1, s1, e1, c2, s2, e2 = f
                else:
                    if not all(math.isfinite(float(x.replace(',', ''))) for x in f[:-2]):
                        fail(f'{where}: nonnumeric MUMmer field')
                    s1, e1, s2, e2 = (x.replace(',', '') for x in f[:4])
                    c1, c2 = f[-2:]
                a, b = interval(c1, s1, e1, lengths, where, closed=current != 6)
                c, d = interval(c2, s2, e2, lengths, where, closed=current != 6)
                if c1 in retained and c2 in retained:
                    print(c1, a, b, c2, c, d, sep='\t')
    for handle in handles.values():
        handle.close()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode', choices=['index', 'gene', 'repeat', 'intact', 'ltr', 'links'])
    p.add_argument('--fai', required=True)
    p.add_argument('--retained', required=True)
    p.add_argument('--input')
    p.add_argument('--output')
    args = p.parse_args()
    try:
        convert(args)
    except (ValueError, OSError) as exc:
        p.exit(1, f'[ERROR] {exc}\n')


if __name__ == '__main__':
    main()
