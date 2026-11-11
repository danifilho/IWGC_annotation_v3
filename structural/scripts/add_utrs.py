#!/usr/bin/env python3
"""Add UTRs to the lab-format GFF3 written by finalize_gff3.py (exons = CDS).

UTR sources, per transcript and per side (5', 3'):
  1. BRAKER: braker_utr.gtf (UTRs from the RNA-seq / Iso-Seq StringTie
     transcripts, BRAKER4 stringtie2utr). Column 2 = the gene's label (BRAKER4).
  2. Helixer (optional): the UTR of the Helixer transcript whose CDS is
     identical, used only for a side BRAKER left without a UTR. Column 2 = Helixer.
Genes added from Helixer (column 2 Helixer) take their own Helixer UTRs.

braker_utr.gtf artefacts (Figuring out A. palmeri, test 9) are removed: UTR
features identical to the stop codon are dropped, the stop codon is clipped
from 3' UTRs that start on it (the CDS of the GFF3 includes it), and a side
whose UTR still overlaps the CDS or lies on the wrong side of it is dropped.
A source is used only if its CDS (+ stop codon) is identical to the GFF3 CDS.

Exons become the union of CDS and UTR segments (adjacent pieces joined),
renumbered 5' -> 3'; mRNA and gene coordinates grow to include the UTRs. UTR
IDs follow the IWGC v1 files: <mRNA>:five_prime_utr, <mRNA>:three_prime_utr.

Usage: add_utrs.py final.gff3 braker_utr.gtf[.gz]|- helixer.gff3|- summary.tsv > final_utr.gff3
"""
import gzip, re, sys
from collections import defaultdict, Counter

def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)

def merge(ivs):
    """sorted union of intervals, joining overlapping or adjacent ones"""
    out = []
    for a, b in sorted(ivs):
        if out and a <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out

def overlaps(iv, ivs):
    return any(iv[0] <= b and iv[1] >= a for a, b in ivs)

def clean_side(segs, cds, strand, side):
    """UTR segments of one side -> list, or [] if they are not a valid UTR for this CDS"""
    if not segs:
        return []
    if any(overlaps(s, cds) for s in segs):
        return []
    lo, hi = cds[0][0], cds[-1][1]
    upstream = (side == '5') == (strand == '+')
    if not all((e < lo) if upstream else (s > hi) for s, e in segs):
        return []
    return merge(segs)

def read_braker(path):
    """transcript -> (CDS incl. stop, 5' segments, 3' segments)"""
    raw = defaultdict(lambda: {'cds': [], 'stop': [], '5': [], '3': []})
    kind = {'CDS': 'cds', 'stop_codon': 'stop', 'five_prime_UTR': '5', 'three_prime_UTR': '3'}
    for line in opener(path):
        f = line.rstrip('\n').split('\t')
        if len(f) != 9 or f[2] not in kind:
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if m:
            raw[m.group(1)][kind[f[2]]].append((int(f[3]), int(f[4])))
    out = {}
    for t, r in raw.items():
        stop = set(r['stop'])
        cds = merge(r['cds'] + [s for s in stop if not overlaps(s, r['cds'])])
        utr = {}
        for side in '53':
            segs = []
            for s, e in set(r[side]):
                if (s, e) in stop:
                    continue  # the stop codon labelled as a UTR
                for a, b in stop:  # 3' UTR starting on the stop codon: clip it
                    if a <= s <= b:
                        s = b + 1
                    if a <= e <= b:
                        e = a - 1
                if s <= e:
                    segs.append((s, e))
            utr[side] = segs
        out[t] = (cds, utr['5'], utr['3'])
    return out

def read_helixer(path):
    """transcript -> (CDS, 5' segments, 3' segments); Helixer CDS include the stop codon"""
    raw = defaultdict(lambda: {'cds': [], '5': [], '3': []})
    kind = {'CDS': 'cds', 'five_prime_UTR': '5', 'three_prime_UTR': '3'}
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) != 9 or f[2] not in kind:
            continue
        t = re.search(r'Parent=([^;]+)', f[8]).group(1)
        raw[t][kind[f[2]]].append((int(f[3]), int(f[4])))
    return {t: (merge(r['cds']), r['5'], r['3']) for t, r in raw.items()}

def main(gff, braker_utr, helixer, summary):
    bk = read_braker(braker_utr) if braker_utr != '-' else {}
    hx = read_helixer(helixer) if helixer != '-' else {}
    hx_by_cds = {}
    for t, (cds, _, _) in hx.items():
        hx_by_cds.setdefault(tuple(cds), t)

    genes, order, mrnas = {}, [], defaultdict(list)
    lines = defaultdict(list)  # mRNA ID -> its child lines
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) != 9:
            continue
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        if f[2] == 'gene':
            genes[a['ID']] = f
            order.append(a['ID'])
        elif f[2] == 'mRNA':
            mrnas[a['Parent']].append(f)
        else:
            lines[a['Parent']].append(f)

    n = Counter()
    out = sys.stdout
    out.write('##gff-version 3\n')
    for gid in order:
        g = genes[gid]
        block, glo, ghi = [], int(g[3]), int(g[4])
        for m in mrnas[gid]:
            a = dict(x.split('=', 1) for x in m[8].split(';') if '=' in x)
            mid, src_id, label, strand = a['ID'], a.get('braker_id'), m[1], m[6]
            cds_lines = [f for f in lines[mid] if f[2] == 'CDS']
            cds = merge((int(f[3]), int(f[4])) for f in cds_lines)
            utr = {'5': ([], None), '3': ([], None)}
            n['mRNAs'] += 1
            if label == 'Helixer':
                h = hx.get(src_id)
                if h and h[0] == cds:
                    for i, side in ((1, '5'), (2, '3')):
                        utr[side] = (clean_side(h[i], cds, strand, side), 'Helixer')
                elif h:
                    n['helixer_cds_mismatch'] += 1
            else:
                b = bk.get(src_id)
                if b and b[0] == cds:
                    for i, side in ((1, '5'), (2, '3')):
                        utr[side] = (clean_side(b[i], cds, strand, side), label)
                elif b:
                    n['braker_cds_mismatch'] += 1
                h = hx.get(hx_by_cds.get(tuple(cds)))
                for i, side in ((1, '5'), (2, '3')):
                    if h and not utr[side][0]:
                        segs = clean_side(h[i], cds, strand, side)
                        if segs:
                            utr[side] = (segs, 'Helixer')
            for side in '53':
                if utr[side][0]:
                    n[f"{side}' UTR from {utr[side][1]}"] += 1
            if utr['5'][0] and utr['3'][0]:
                n['mRNAs with both UTRs'] += 1
            exons = merge(cds + utr['5'][0] + utr['3'][0])
            lo, hi = exons[0][0], exons[-1][1]
            glo, ghi = min(glo, lo), max(ghi, hi)
            seq = m[0]
            m = m[:3] + [str(lo), str(hi)] + m[5:]
            block.append(m)
            segs = exons if strand == '+' else exons[::-1]
            for k, (s, e) in enumerate(segs, 1):
                block.append([seq, label, 'exon', str(s), str(e), '.', strand, '.', f"ID={mid}:{k};Parent={mid}"])
            block += cds_lines
            for side, ftype, tag in (('5', 'five_prime_UTR', 'five_prime_utr'), ('3', 'three_prime_UTR', 'three_prime_utr')):
                ss, src = utr[side]
                for s, e in (ss if strand == '+' else ss[::-1]):
                    block.append([seq, src, ftype, str(s), str(e), '.', strand, '.', f"ID={mid}:{tag};Parent={mid}"])
        g = g[:3] + [str(glo), str(ghi)] + g[5:]
        out.write('\t'.join(g) + '\n')
        for f in block:
            out.write('\t'.join(f) + '\n')
    with open(summary, 'w') as fh:
        for k in sorted(n):
            fh.write(f"{k}\t{n[k]}\n")
    sys.stderr.write(' '.join(f"{k.replace(' ', '_')}={v}" for k, v in sorted(n.items())) + '\n')

if __name__ == '__main__':
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
