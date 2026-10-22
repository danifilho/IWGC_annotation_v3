#!/usr/bin/env python3
"""Write the IWGC functional GFF3 (same layout as the IWGC v1 releases, e.g.
SalTr_v01.0.gff) from the structural GFF3 + the functional pipeline results.

Functional information goes on the GENE line:
  ID=SalTrChr01Ag000010;Name=SalTrChr01Ag000010;
  Dbxref=XP_021720546.1:110688150:786.0            NCBI best hit : GeneID : bitscore
  Ontology_term=GO:0003700%2CGO:0005667            InterProScan GO terms
  Note=transcription factor E2FB-like [Chenopodium quinoa]:(InterPro|IPR015633|E2F Family)
       NCBI description of the best hit, then the best InterProScan match
       (an InterPro entry if any; otherwise a member-database signature)
  Genes without an NCBI hit: Note=Unknown Protein (Dbxref omitted).
Child features keep only ID and Parent; every feature gets SOURCE in column 2.
Feature order is the order of the input GFF3 (assembly order).

GeneID: the NCBI protein FASTA used by the pipeline has no GeneIDs. Give
--geneid-map (2 columns: protein accession <TAB> GeneID, e.g. made from NCBI
gene2refseq) for real GeneIDs; without it the LOC number in the description is
used when there is one (it is the GeneID), otherwise "NA".

Usage:
  iwgc_functional_gff.py --gff structural.gff3 --tsv functional_annotation.tsv \\
      --iprscan SP.iprscan.tsv --species PhaMi --version v3 --source IWGC_v3 \\
      [--geneid-map protein2geneid.tsv] > SP.functional.gff3
"""
import argparse, datetime, re, sys
from urllib.parse import quote

def enc(s):
    """GFF3 column-9 escaping of reserved characters (; = & , and controls)."""
    return quote(s, safe=" :|()[]/'-+.*_~!@#$^{}<>?\"")

def parse_attrs(s):
    d = {}
    for part in s.strip().rstrip(';').split(';'):
        if '=' in part:
            k, v = part.split('=', 1)
            d[k.strip()] = v
    return d

def read_tsv(path):
    """functional_annotation.tsv -> {mRNA id: (subject, bitscore, description, [GO])}"""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        col = {name: i for i, name in enumerate(header)}
        i_desc = col.get('NCBI DESCRIPTION', 1)
        i_ncbi = next((i for n, i in col.items() if n.startswith('NCBI Subject')), 2)
        i_go = col.get('IPRSCAN GO', 8)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) <= max(i_desc, i_ncbi, i_go):
                continue
            hit = f[i_ncbi].split(';')
            subject = hit[0] if hit and hit[0] not in ('N/A', '') else None
            bits = hit[2] if len(hit) > 2 and hit[2] not in ('N/A', '') else None
            desc = f[i_desc].strip()
            if desc in ('N/A', '', 'Description not found'):
                desc = None
            gos = []
            if f[i_go] not in ('N/A', ''):
                for g in f[i_go].split(';'):
                    m = re.match(r'(GO:\d+)', g.strip())
                    if m and m.group(1) not in gos:
                        gos.append(m.group(1))
            out[f[0]] = (subject, bits, desc, gos)
    return out

def best_interpro(path):
    """InterProScan TSV -> {protein: (db, accession, description)}.
    Prefers an InterPro entry (cols 12/13); otherwise the member-database
    signature (cols 4/5/6). Among candidates the lowest e-value wins."""
    def ev(x):
        try:
            return float(x)
        except ValueError:
            return float('inf')
    ipr, sig = {}, {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 11:
                continue
            prot, e = f[0], ev(f[8])
            if len(f) > 12 and f[11] not in ('', '-'):
                cand = (e, 'InterPro', f[11], f[12] if f[12] not in ('', '-') else '')
                if prot not in ipr or cand[0] < ipr[prot][0]:
                    ipr[prot] = cand
            else:
                d = f[5] if f[5] not in ('', '-') else ''
                cand = (e, f[3], f[4], d)
                if prot not in sig or cand[0] < sig[prot][0]:
                    sig[prot] = cand
    best = {}
    for prot in set(ipr) | set(sig):
        c = ipr.get(prot) or sig[prot]
        best[prot] = c[1:]
    return best

def read_geneid_map(path):
    m = {}
    if path:
        with open(path) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 2:
                    m[f[0]] = f[1]
                    m.setdefault(f[0].split('.')[0], f[1])   # also without version
    return m

def fmt_bits(b):
    try:
        return f"{float(b):.1f}"
    except (TypeError, ValueError):
        return "NA"

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gff', required=True)
    ap.add_argument('--tsv', required=True)
    ap.add_argument('--iprscan', required=True)
    ap.add_argument('--species', required=True)
    ap.add_argument('--version', default='v3')
    ap.add_argument('--source', default='IWGC_v3')
    ap.add_argument('--tools', default='InterProScan, MultiLoc2, MMseqs2')
    ap.add_argument('--geneid-map')
    a = ap.parse_args()

    func = read_tsv(a.tsv)
    ipr = best_interpro(a.iprscan)
    gid_map = read_geneid_map(a.geneid_map)

    # pass 1: mRNA -> gene, so each gene takes the annotation of its annotated mRNA
    gene_ann = {}
    with open(a.gff) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) == 9 and f[2] in ('mRNA', 'transcript'):
                at = parse_attrs(f[8])
                if at.get('ID') in func and at.get('Parent'):
                    gene_ann.setdefault(at['Parent'], at['ID'])

    out = sys.stdout
    out.write('##gff-version 3\n')
    out.write(f"##Generated using iwgc_functional_gff.py, {a.species}_{a.version}, {datetime.datetime.now().isoformat(sep=' ', timespec='seconds')}\n")
    out.write(f"##Project Name : {a.species}\n##Job Name  : IWGC Annotation Pipeline\n##Tool      : {a.tools}\n")

    n = {'gene': 0, 'hit': 0, 'unknown': 0, 'go': 0, 'ipr': 0, 'geneid': 0}
    with open(a.gff) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 9:
                continue
            f[1] = a.source
            at = parse_attrs(f[8])
            if f[2] == 'gene':
                n['gene'] += 1
                gid = at.get('ID', '')
                attrs = [f"ID={gid}", f"Name={at.get('Name', gid)}"]
                mid = gene_ann.get(gid)
                subject, bits, desc, gos = func.get(mid, (None, None, None, []))
                if subject:
                    n['hit'] += 1
                    geneid = gid_map.get(subject) or gid_map.get(subject.split('.')[0])
                    if not geneid and desc:
                        m = re.search(r'\bLOC(\d+)\b', desc)
                        geneid = m.group(1) if m else None
                    if geneid:
                        n['geneid'] += 1
                    attrs.append(f"Dbxref={enc(subject)}:{geneid or 'NA'}:{fmt_bits(bits)}")
                if gos:
                    n['go'] += 1
                    attrs.append("Ontology_term=" + '%2C'.join(gos))
                if subject:
                    note = enc(desc) if desc else 'Unknown Protein'
                    if mid in ipr:
                        db, acc, idesc = ipr[mid]
                        n['ipr'] += 1
                        note += f":({enc(db)}|{enc(acc)}|{enc(idesc)})"
                else:
                    n['unknown'] += 1
                    note = 'Unknown Protein'
                attrs.append(f"Note={note}")
                f[8] = ';'.join(attrs)
            else:
                keep = [f"{k}={at[k]}" for k in ('ID', 'Parent') if k in at]
                f[8] = ';'.join(keep)
            out.write('\t'.join(f) + '\n')

    g = max(n['gene'], 1)
    sys.stderr.write(
        f"genes: {n['gene']}  with NCBI hit: {n['hit']} ({100*n['hit']/g:.1f}%)  "
        f"Unknown Protein: {n['unknown']}  with GO: {n['go']}  "
        f"Note with InterProScan match: {n['ipr']}  Dbxref with GeneID: {n['geneid']}\n")

if __name__ == '__main__':
    main()
