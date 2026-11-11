#!/usr/bin/env python3
"""Add Helixer genes at loci BRAKER missed.

TSEBRA was run with --keep_gtf braker.gtf --gtf helixer.gtf: it keeps every
BRAKER transcript and adds the Helixer transcripts that its hints support.
Of those additions, only Helixer transcripts whose CDS span overlaps no BRAKER
CDS span (either strand) are kept: Helixer isoforms at BRAKER loci are dropped,
so the BRAKER models are never changed. Their Helixer GTF lines (original IDs,
source column "Helixer") are appended to braker.gtf.

Usage: helixer_merge.py braker.gtf tsebra.gtf helixer.gtf merged.gtf added.tsv
"""
import re, sys, bisect
from collections import defaultdict

def chains(path):
    """transcript -> (seq, strand, sorted CDS tuple), and transcript -> gene"""
    cds, gene = defaultdict(list), {}
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) != 9 or f[2] != 'CDS':
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        gene[t] = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        cds[t].append((f[0], f[6], int(f[3]), int(f[4])))
    return {t: (v[0][0], v[0][1], tuple(sorted((c[2], c[3]) for c in v))) for t, v in cds.items()}, gene

def main(braker, tsebra, helixer, merged, table):
    bk, _ = chains(braker)
    ts, _ = chains(tsebra)
    hx, hgene = chains(helixer)
    bset = set(bk.values())
    by_chain = {ch: t for t, ch in hx.items()}
    spans = defaultdict(list)
    for s, _, ex in bset:
        spans[s].append((ex[0][0], ex[-1][1]))
    for s in spans:
        spans[s].sort()
    longest = max((e - b for v in spans.values() for b, e in v), default=0)

    def overlaps(s, b, e):
        iv = spans.get(s, [])
        for b2, e2 in iv[bisect.bisect_left(iv, (b - longest - 1, 0)):]:
            if b2 > e:
                return False
            if e2 >= b:
                return True
        return False

    added, at_locus = [], 0
    for t, ch in ts.items():
        if ch in bset or ch not in by_chain:
            continue  # a BRAKER model, or not from Helixer
        if overlaps(ch[0], ch[2][0][0], ch[2][-1][1]):
            at_locus += 1
            continue
        added.append(by_chain[ch])
    keep = set(added)
    with open(merged, 'w') as out:
        out.write(open(braker).read())
        for line in open(helixer):
            m = re.search(r'transcript_id "([^"]+)"', line)
            if m and m.group(1) in keep:
                out.write(line)
    with open(table, 'w') as out:
        out.write('helixer_transcript\thelixer_gene\tseqid\tstrand\tstart\tend\tcds_exons\n')
        for t in sorted(added, key=lambda t: (hx[t][0], hx[t][2][0][0])):
            s, st, ex = hx[t]
            out.write(f"{t}\t{hgene[t]}\t{s}\t{st}\t{ex[0][0]}\t{ex[-1][1]}\t{len(ex)}\n")
    sys.stderr.write(f"BRAKER transcripts={len(bset)} helixer added at new loci={len(added)} "
                     f"helixer alternatives at BRAKER loci dropped={at_locus}\n")

if __name__ == '__main__':
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
