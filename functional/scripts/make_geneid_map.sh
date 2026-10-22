#!/bin/bash
# Build the optional NCBI protein -> GeneID table used by iwgc_functional_gff.py
# (config key ncbi_geneid_map) for the GeneID field of Dbxref.
# Downloads NCBI gene2refseq (~2.4 GB), keeps only the accessions that are best
# hits, and deletes the download. Run on a node with internet access.
#
# Usage: make_geneid_map.sh <out.tsv> <best_hit.out> [<best_hit.out> ...]
#   e.g. make_geneid_map.sh inputs/ncbi_protein2geneid.tsv mmseqs_output/PhaMi_ncbi_best_hit.out
set -euo pipefail
OUT="$1"; shift
TMP=$(mktemp -d "${TMPDIR:-.}/gene2refseq.XXXX")
trap 'rm -rf "$TMP"' EXIT

cut -f2 "$@" | sort -u > "$TMP/hits.txt"
echo "distinct best-hit accessions: $(wc -l < "$TMP/hits.txt")"
wget -q -O "$TMP/gene2refseq.gz" https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz
# gene2refseq columns: 2 = GeneID, 6 = protein_accession.version
zcat "$TMP/gene2refseq.gz" \
  | awk -F'\t' 'NR==FNR { want[$1]; next } ($6 in want) && !seen[$6]++ { print $6 "\t" $2 }' "$TMP/hits.txt" - > "$OUT"
echo "accessions with a GeneID: $(wc -l < "$OUT")  ->  $OUT"
