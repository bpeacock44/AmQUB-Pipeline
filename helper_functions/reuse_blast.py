#!/usr/bin/env python3
"""One-off reuse of BLAST format-7 query blocks, matched by exact nucleotide sequence."""
import argparse
from collections import OrderedDict
from pathlib import Path


def fasta(path):
    result = OrderedDict()
    name, seq = None, []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith('>'):
                if name is not None:
                    result[name] = ''.join(seq).upper()
                name, seq = line[1:].split()[0], []
                if name in result:
                    raise ValueError(f'Duplicate FASTA ID: {name}')
            elif line:
                if name is None:
                    raise ValueError(f'Sequence before FASTA header: {path}')
                seq.append(line)
    if name is not None:
        result[name] = ''.join(seq).upper()
    return result


def blocks(path):
    """Select the largest complete query block when legacy output repeats queries."""
    result = {}
    current, lines = None, []
    def store():
        if current is None:
            return
        hits = sum(bool(x.strip()) and not x.startswith('#') for x in lines)
        if current not in result or hits > result[current][0]:
            result[current] = (hits, lines[:])
    with open(path) as fh:
        for line in fh:
            if line.startswith('# Query:'):
                store()
                current = line.partition(':')[2].strip().split()[0]
                lines = [line]
            elif current is not None:
                if line.startswith('# BLAST processed') or line.startswith('# BLASTN'):
                    continue
                if not line.startswith('#') and line.strip():
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) != 10 or fields[0] != current:
                        raise ValueError(f'Invalid BLAST hit in query {current}: {line[:100]}')
                lines.append(line)
        store()
    return {k: v[1] for k, v in result.items()}


def main():
    p = argparse.ArgumentParser()
    for arg in ('old-fasta', 'new-fasta', 'old-blast', 'new-only', 'reused-blast'):
        p.add_argument('--' + arg, required=True)
    a = p.parse_args()
    old = fasta(a.old_fasta)
    new = fasta(a.new_fasta)
    hits = blocks(a.old_blast)
    by_seq = {}
    for ident, sequence in old.items():
        if ident in hits:
            by_seq.setdefault(sequence, ident)
    reused = 0
    with open(a.new_only, 'w') as fresh, open(a.reused_blast, 'w') as output:
        for ident, sequence in new.items():
            old_ident = by_seq.get(sequence)
            if old_ident is None:
                fresh.write(f'>{ident}\n{sequence}\n')
                continue
            reused += 1
            for line in hits[old_ident]:
                if line.startswith('# Query:'):
                    output.write(f'# Query: {ident}\n')
                elif line.startswith('#'):
                    output.write(line)
                else:
                    fields = line.rstrip('\n').split('\t')
                    fields[0] = ident
                    output.write('\t'.join(fields) + '\n')
    print(f'Reused {reused}/{len(new)} sequences; {len(new)-reused} require initial BLAST.')


if __name__ == '__main__':
    main()
