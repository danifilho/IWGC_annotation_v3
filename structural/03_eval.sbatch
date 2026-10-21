#!/bin/bash --login
#SBATCH --job-name=v3_eval
#SBATCH --partition=general-long
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/mnt/scratch/dasilvaf/iwgc_pipeline_v3/logs/%x_%j.out
# Gene-set stats + BUSCO (poales_odb10, proteins) for one attempt, on the
# longest isoform per gene. Usage: sbatch 03_eval.sbatch A_isoseq|B_protein
# Proteins > 10,000 aa are dropped first: John's merged set had 216 proteins
# > 20k aa, and hmmsearch aborts on anything > 100k aa (that was his 0% BUSCO).
set -euo pipefail
source "${SLURM_SUBMIT_DIR}/config.sh"
A="${1:?usage: sbatch 03_eval.sbatch A_isoseq|B_protein}"
RUN="$WORK/braker_$A"
OUT="$FINAL/$A"
mkdir -p "$OUT"
[ -s "$RUN/braker.aa" ] || { echo "ERROR: $RUN/braker.aa missing"; exit 1; }
cp "$RUN/braker.gtf" "$RUN/braker.aa" "$RUN/braker.codingseq" "$OUT/"
[ -s "$RUN/braker.gff3" ] && cp "$RUN/braker.gff3" "$OUT/"

# longest isoform per gene (headers are gN.tM), drop > 10k aa
awk '/^>/{if(h)print h"\t"s; h=substr($1,2); s=""; next}{s=s $0}END{print h"\t"s}' "$OUT/braker.aa" \
  | awk -F'\t' '{g=$1; sub(/\.t[0-9]+$/,"",g); l=length($2); if(l>best[g]){best[g]=l; seq[g]=$2; id[g]=$1}}
      END{for(g in seq){ if(best[g]<=10000) print ">"id[g]"\n"seq[g]; else n++ } print n > "/dev/stderr"}' \
  > "$OUT/braker.longest.aa" 2> "$OUT/n_over_10k.txt"

{
  echo "attempt: $A"
  echo "genes: $(awk -F'\t' '$3=="gene"' "$OUT/braker.gtf" | wc -l)"
  echo "transcripts: $(grep -c '^>' "$OUT/braker.aa")"
  echo "proteins > 10k aa (dropped for BUSCO): $(cat "$OUT/n_over_10k.txt")"
  awk -F'\t' '$3=="CDS"{match($9,/transcript_id "[^"]+"/); t=substr($9,RSTART,RLENGTH); n[t]++}
    END{for(t in n){c++; e+=n[t]; if(n[t]==1)s++}
        printf "single-exon transcripts: %d / %d (%.1f%%)\nmean CDS exons per transcript: %.2f\n", s, c, 100*s/c, e/c}' "$OUT/braker.gtf"
  awk '/^>/{if(l)print l; l=0; next}{l+=length($0)}END{print l}' "$OUT/braker.longest.aa" | sort -n \
    | awk '{a[NR]=$1}END{print "longest-isoform protein length median: " a[int(NR/2)+1]}'
  grep -iE "accuracy" "$RUN/braker.log" || true
} | tee "$OUT/stats.txt"

module purge
module load "$BUSCO_MODULE"
cd "$OUT"
busco -i braker.longest.aa -m proteins -l "$BUSCO_LINEAGE" --offline -c "$SLURM_CPUS_PER_TASK" -o busco_poales -f
cp busco_poales/short_summary*.txt "$OUT/" 2>/dev/null || true
# BUSCO leaves ~10k small files per run; keep only the summaries + full table
tar -czf busco_poales.tar.gz busco_poales && rm -rf busco_poales
grep "C:" "$OUT"/short_summary*.txt
check_files
