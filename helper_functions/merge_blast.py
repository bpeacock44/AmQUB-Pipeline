#!/usr/bin/env python3
"""Combine format-7 BLAST outputs, replacing first-pass query blocks with expanded searches."""
import argparse
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--initial', nargs='+', required=True)
p.add_argument('--expanded', required=True)
p.add_argument('--output', required=True)
a = p.parse_args()

def blocks(path):
    current = None
    block = []
    for line in Path(path).read_text().splitlines(keepends=True):
        if line.startswith('# BLAST processed'):
            continue
        if line.startswith('# Query:'):
            if current is not None:
                yield current, block
            current = line.partition(':')[2].strip().split()[0]
            block = [line]
        elif current is not None:
            block.append(line)
        elif line.strip() and not line.startswith('#'):
            raise ValueError(f'Hit outside # Query block in {path}')
    if current is not None:
        yield current, block

expanded = dict(blocks(a.expanded)) if Path(a.expanded).exists() else {}
seen = set()
with open(a.output, 'w') as out:
    for path in a.initial:
        for ident, block in blocks(path):
            if ident in seen:
                raise ValueError(f'Duplicate initial query {ident}')
            seen.add(ident)
            out.writelines(expanded.pop(ident, block))
    for ident, block in expanded.items():
        if ident not in seen:
            raise ValueError(f'Expanded query {ident} absent from initial output')
    out.write('# BLAST processed\n')
