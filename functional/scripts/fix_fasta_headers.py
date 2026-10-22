#!/usr/bin/env python3
"""Repair a protein FASTA before it is used as an MMseqs2 target database.

Fixes the two defects found in the IWGC GenomesDatabase.v1.prot:
  1. headers glued to the end of the previous sequence line (no newline before
     '>'), e.g. "...MKVLLA>AA051080-RA". MMseqs only starts a record at a '>'
     at the START of a line, so every glued record is silently merged into the
     previous sequence. In GenomesDatabase.v1.prot this turned 49,542 proteins
     into ONE 18.9 M aa "protein" (AA049522-RA), which became the best hit of
     most grass genes.
  2. a first header without its '>'.
Sequences are rewrapped at 60 characters.

Usage: fix_fasta_headers.py in.fasta out.fasta
Prints the number of records and of repaired headers to stderr.
"""
import re, sys

SEQ = re.compile(r'^[A-Za-z*.\-]*$')     # a protein sequence line

def main(inp, outp):
    out = open(outp, 'w')
    seq, n_rec, n_glued, n_first = [], 0, 0, 0

    def flush():
        s = ''.join(seq)
        for i in range(0, len(s), 60):
            out.write(s[i:i + 60] + '\n')
        seq.clear()

    def header(h):
        nonlocal n_rec
        flush()
        out.write('>' + h + '\n')
        n_rec += 1

    with open(inp) as fh:
        for i, line in enumerate(fh):
            line = line.rstrip('\r\n')
            if i == 0 and line and not line.startswith('>') and not SEQ.match(line):
                header(line)                     # first header missing its '>'
                n_first += 1
                continue
            if line.startswith('>'):
                header(line[1:])
                continue
            k = line.find('>')
            if k == -1:
                seq.append(line)
            else:                                # sequence + glued header(s)
                parts = line.split('>')
                seq.append(parts[0])
                for p in parts[1:]:
                    header(p)
                    n_glued += 1
    flush()
    out.close()
    sys.stderr.write(f"records: {n_rec}  glued headers repaired: {n_glued}  "
                     f"missing first '>' added: {n_first}\n")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
