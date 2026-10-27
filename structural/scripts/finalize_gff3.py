#!/usr/bin/env python3
"""braker.gtf -> clean, sorted GFF3 with lab-standard PhaMi IDs (same scheme as v2):

  gene  PhaMiChr01Cg000010             (counter +10 across the whole genome,
  mRNA  PhaMiChr01Cg000010-mRNA-1       never reset per chromosome; genes in
  exon  PhaMiChr01Cg000010-mRNA-1:1     genome FASTA order, then by start)
  CDS   PhaMiChr01Cg000010-mRNA-1:cds

braker.gtf mixes AUGUSTUS, GeneMark.hmm3 and gmst records (TSEBRA output).
Only CDS + stop_codon lines are used: BRAKER predicts no UTRs, so exons = CDS.
The stop codon is folded into the CDS when a source left it out (GeneMark
excludes it, AUGUSTUS includes it), so every CDS ends with its stop codon.
Isoforms are numbered by start coordinate, then by original transcript ID.
CDS phases are recomputed from the cumulative length, starting from the GTF
frame of the 5'-most CDS: 0 for complete models, 1 or 2 for the few GeneMark
models that are partial at the 5' end (start mid-codon).

Optional EXCLUDE file: BRAKER transcript IDs (one per line) to leave out, e.g.
models whose braker.aa protein has an internal stop codon (GeneMark.hmm3 emits
a few). They are removed BEFORE numbering, so gene IDs stay contiguous; genes
left without transcripts disappear.

Optional SOURCE: column-2 label (default BRAKER3; BRAKER4 for the BRAKER4 predictor).

Usage: finalize_gff3.py braker.gtf genome.fa.fai ORG [EXCLUDE] [SOURCE] > final.gff3
"""
import re, sys
from collections import defaultdict

def main(gtf, fai, org, exclude=None, source='BRAKER3'):
    skip = set()
    if exclude:
        skip = {l.split()[0] for l in open(exclude) if l.strip()}
    order = [l.split('\t')[0] for l in open(fai) if l.strip()]
    rank = {c: i for i, c in enumerate(order)}
    tx = defaultdict(lambda: {'cds': [], 'stop': [], 'frame': {}, 'seq': None, 'strand': None, 'gene': None})
    for line in open(gtf):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) != 9 or f[2] not in ('CDS', 'stop_codon'):
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        if t in skip:
            continue
        g = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        r = tx[t]
        r['seq'], r['strand'], r['gene'] = f[0], f[6], g
        iv = (int(f[3]), int(f[4]))
        (r['cds'] if f[2] == 'CDS' else r['stop']).append(iv)
        if f[2] == 'CDS':
            r['frame'][iv] = int(f[7]) if f[7] in '012' else 0

    genes = defaultdict(list)
    n_fold = 0
    for t, r in tx.items():
        if not r['cds']:
            continue
        cds = sorted(set(r['cds']))
        # 5' CDS frame, before the stop codon is folded into the 3' segment
        r['phase5'] = r['frame'][cds[0] if r['strand'] == '+' else cds[-1]]
        for s, e in r['stop']:
            if any(a <= s and e <= b for a, b in cds):
                continue
            # extend the terminal CDS segment that abuts the stop codon
            if r['strand'] == '+' and cds[-1][1] + 1 == s:
                cds[-1] = (cds[-1][0], e); n_fold += 1
            elif r['strand'] == '-' and cds[0][0] - 1 == e:
                cds[0] = (s, cds[0][1]); n_fold += 1
            else:
                sys.stderr.write(f"WARNING: stop codon of {t} not adjacent to CDS, left out\n")
        r['cds'] = cds
        genes[(r['seq'], r['gene'])].append((t, r))

    blocks = []
    for (seq, g), txs in genes.items():
        start = min(r['cds'][0][0] for _, r in txs)
        end = max(r['cds'][-1][1] for _, r in txs)
        blocks.append((rank.get(seq, len(rank)), start, end, seq, txs))
    blocks.sort(key=lambda b: (b[0], b[1], b[2]))

    out = sys.stdout
    out.write('##gff-version 3\n')
    counter = 0
    n_tx = 0
    for _, start, end, seq, txs in blocks:
        counter += 10
        gid = f"{org}{seq}g{counter:06d}"
        strand = txs[0][1]['strand']
        out.write(f"{seq}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid};Name={gid}\n")
        txs.sort(key=lambda x: (x[1]['cds'][0][0], x[1]['cds'][-1][1], x[0]))
        for i, (old, r) in enumerate(txs, 1):
            n_tx += 1
            mid = f"{gid}-mRNA-{i}"
            cds = r['cds']
            out.write(f"{seq}\t{source}\tmRNA\t{cds[0][0]}\t{cds[-1][1]}\t.\t{strand}\t.\t"
                      f"ID={mid};Name={mid};Parent={gid};braker_id={old}\n")
            # exons/CDS listed 5'->3' along the transcript; phase from cumulative length
            segs = cds if strand == '+' else cds[::-1]
            for k, (a, b) in enumerate(segs, 1):
                out.write(f"{seq}\t{source}\texon\t{a}\t{b}\t.\t{strand}\t.\tID={mid}:{k};Parent={mid}\n")
            # start so that the first segment gets the 5' frame of the GTF
            done = (3 - r['phase5']) % 3
            for a, b in segs:
                phase = (3 - done % 3) % 3
                out.write(f"{seq}\t{source}\tCDS\t{a}\t{b}\t.\t{strand}\t{phase}\tID={mid}:cds;Parent={mid}\n")
                done += b - a + 1
    sys.stderr.write(f"genes={len(blocks)} transcripts={n_tx} stop_codons_folded_into_CDS={n_fold} "
                     f"excluded_transcripts={len(skip)}\n")

if __name__ == '__main__':
    if len(sys.argv) not in (4, 5, 6):
        sys.exit(__doc__)
    main(*sys.argv[1:])
