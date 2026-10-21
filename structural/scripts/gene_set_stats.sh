#!/bin/bash
# Summary statistics and sanity checks for the final gene set.
# Usage: gene_set_stats.sh final.gff3 proteins.fa cds.fa braker.log
GFF="$1"; PROT="$2"; CDS="$3"; LOG="$4"

echo "genes: $(awk -F'\t' '$3=="gene"' "$GFF" | wc -l)"
echo "mRNAs: $(awk -F'\t' '$3=="mRNA"' "$GFF" | wc -l)"
awk -F'\t' '$3=="CDS" { match($9, /Parent=[^;]+/); t = substr($9, RSTART+7, RLENGTH-7); n[t]++ }
  END { for (t in n) { c++; e += n[t]; if (n[t] == 1) s++ }
        printf "single-exon transcripts: %d / %d (%.1f%%)\nmean CDS exons per transcript: %.2f\n", s, c, 100*s/c, e/c }' "$GFF"
awk '/^>/{ if (l) print l; l=0; next } { l += length($0) } END { print l }' "$PROT" | sort -n \
  | awk '{ a[NR]=$1 } END { print "median protein length, all isoforms (aa): " a[int(NR/2)+1] "\nlongest protein (aa): " a[NR] }'
echo "proteins with internal stop: $(awk '/^>/{next} /\*/{n++} END{print n+0}' "$PROT")"
echo "proteins not starting with M: $(awk '/^>/{getline; if (substr($0,1,1) != "M") n++} END{print n+0}' "$PROT")"
echo "CDS length not multiple of 3: $(awk '/^>/{ if (h && l%3) n++; h=1; l=0; next } { l += length($0) } END { if (l%3) n++; print n+0 }' "$CDS")"
echo "genes per sequence:"
awk -F'\t' '$3=="gene" { c[$1]++; if (!($1 in o)) o[$1]=++k } END { for (s in c) print o[s] "\t" s "\t" c[s] }' "$GFF" \
  | sort -k1,1n | cut -f2- | sed 's/^/  /'
if [ -s "$LOG" ]; then
  echo "BRAKER:"
  grep -E "will be executed in|The accuracy after" "$LOG" | sed -E 's/^# .*[0-9]{4}: ?/  /'
fi
