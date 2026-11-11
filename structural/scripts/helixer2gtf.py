#!/usr/bin/env python3
"""Helixer GFF3 -> AUGUSTUS-style GTF, the format TSEBRA and finalize_gff3.py read.

Helixer CDS include the stop codon (GFF3 convention); BRAKER GTFs leave it out
and give it as a stop_codon line. So the last 3 CDS bases become stop_codon
lines (split over two exons if the codon is), and the CDS is written without
them. Also writes start_codon and intron lines, and keeps the GFF3 phase as the
GTF frame. UTRs are not used (BRAKER models have none either).

Usage: helixer2gtf.py helixer.gff3 genome.fa > helixer.gtf
"""
import sys
from collections import defaultdict

STOPS = ('TAA', 'TAG', 'TGA')
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')

def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith('>'):
            if name:
                seqs[name] = ''.join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = ''.join(buf)
    return seqs

def ranges(pos):
    """sorted genomic positions -> [(start, end)] of contiguous runs"""
    pos = sorted(pos)
    out, s, prev = [], pos[0], pos[0]
    for p in pos[1:]:
        if p != prev + 1:
            out.append((s, prev))
            s = p
        prev = p
    out.append((s, prev))
    return out

def main(gff, fasta):
    seqs = read_fasta(fasta)
    parent, cds = {}, defaultdict(list)
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) != 9:
            continue
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        if f[2] == 'mRNA':
            parent[a['ID']] = a['Parent']
        elif f[2] == 'CDS':
            cds[a['Parent']].append((f[0], f[6], int(f[3]), int(f[4]), f[7]))
    out = sys.stdout
    n = n_stop = 0
    for t in sorted(cds, key=lambda t: (cds[t][0][0], min(c[2] for c in cds[t]))):
        segs = sorted(cds[t], key=lambda c: c[2])
        seq, strand, g = segs[0][0], segs[0][1], parent[t]
        # CDS bases in transcript order (5' -> 3')
        pos = [p for c in segs for p in range(c[2], c[3] + 1)]
        if strand == '-':
            pos = pos[::-1]
        codon = ''.join(seqs[seq][p - 1] for p in pos[-3:]).upper()
        if strand == '-':
            codon = codon.translate(COMP)
        has_stop = codon in STOPS
        core = pos[:-3] if has_stop else pos
        five = segs[0] if strand == '+' else segs[-1]
        frame = five[4] if five[4] in '012' else '0'
        att = f'transcript_id "{t}"; gene_id "{g}";'
        lo, hi = min(pos), max(pos)
        out.write(f"{seq}\tHelixer\tgene\t{lo}\t{hi}\t.\t{strand}\t.\t{g}\n")
        out.write(f"{seq}\tHelixer\ttranscript\t{lo}\t{hi}\t.\t{strand}\t.\t{t}\n")
        ex = ranges(core)
        # frame of each CDS piece from the cumulative length, 5' -> 3'
        order = ex if strand == '+' else ex[::-1]
        done, frames = (3 - int(frame)) % 3, {}
        for a, b in order:
            frames[(a, b)] = (3 - done % 3) % 3
            done += b - a + 1
        for a, b in ex:
            out.write(f"{seq}\tHelixer\tCDS\t{a}\t{b}\t.\t{strand}\t{frames[(a, b)]}\t{att}\n")
            out.write(f"{seq}\tHelixer\texon\t{a}\t{b}\t.\t{strand}\t.\t{att}\n")
        for (a0, b0), (a1, b1) in zip(ex, ex[1:]):
            out.write(f"{seq}\tHelixer\tintron\t{b0 + 1}\t{a1 - 1}\t.\t{strand}\t.\t{att}\n")
        for a, b in ranges(pos[:3]):
            out.write(f"{seq}\tHelixer\tstart_codon\t{a}\t{b}\t.\t{strand}\t0\t{att}\n")
        if has_stop:
            n_stop += 1
            for a, b in ranges(pos[-3:]):
                out.write(f"{seq}\tHelixer\tstop_codon\t{a}\t{b}\t.\t{strand}\t0\t{att}\n")
        n += 1
    sys.stderr.write(f"helixer transcripts={n} with_stop_codon={n_stop}\n")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
