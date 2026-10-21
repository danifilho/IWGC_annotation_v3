#!/usr/bin/env python3
"""Write the subject sequences of the MMseqs best hits from a (large) protein
FASTA, for agat_sp_manage_functional_annotation.pl --db.

Also repairs a known defect of the NCBI GenomesDatabase FASTA: some headers are
glued to the end of the previous sequence line (no newline before '>'), and the
first header may lack its '>'. Only the needed records are written, so the
3+ GB database never has to be rewritten.

Usage: extract_hits.py best_hit.out database.fasta > hits.fasta
"""
import sys

def main(best_hits, fasta):
    want = set()
    with open(best_hits) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) > 1:
                want.add(f[1])
    out = sys.stdout
    keep = False
    found = 0

    def header(h):
        nonlocal keep, found
        keep = h.split()[0] in want
        if keep:
            found += 1
            out.write('>' + h + '\n')

    with open(fasta) as fh:
        first = True
        for line in fh:
            line = line.rstrip('\r\n')
            if first:
                first = False
                if line and not line.startswith('>'):
                    line = '>' + line          # first header missing '>'
            if line.startswith('>'):
                header(line[1:])
                continue
            k = line.find('>')
            if k == -1:
                if keep:
                    out.write(line + '\n')
            else:                               # sequence + glued header
                if keep and k:
                    out.write(line[:k] + '\n')
                header(line[k + 1:])
    sys.stderr.write(f"best-hit subjects: {len(want)}  written: {found}\n")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
