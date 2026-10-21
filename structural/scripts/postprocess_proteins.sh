#!/bin/bash
# gffread -S protein FASTA -> (1) all isoforms without the terminal '*',
# (2) longest isoform per gene (IDs <gene>-mRNA-N), both wrapped at 60 aa.
# Usage: postprocess_proteins.sh raw.faa proteins.fa longest.proteins.fa
set -euo pipefail
RAW="$1"; ALL="$2"; LONGEST="$3"

awk '/^>/{ if (h) emit(); h=$1; s=""; next } { s = s $0 }
     END { emit() }
     function emit() { sub(/\*$/, "", s); print h; while (length(s)) { print substr(s,1,60); s = substr(s,61) } }' \
  "$RAW" > "$ALL"

awk '/^>/{ if (h) print h "\t" s; h = substr($1,2); s=""; next } { s = s $0 } END { print h "\t" s }' "$ALL" \
  | awk -F'\t' '{ g=$1; sub(/-mRNA-[0-9]+$/, "", g); if (length($2) > L[g]) { L[g]=length($2); S[g]=$2; I[g]=$1 } }
                END { for (g in S) print I[g] "\t" S[g] }' \
  | sort -k1,1 \
  | awk -F'\t' '{ print ">" $1; s=$2; while (length(s)) { print substr(s,1,60); s=substr(s,61) } }' > "$LONGEST"
