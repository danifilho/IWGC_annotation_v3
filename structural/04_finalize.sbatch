#!/bin/bash --login
#SBATCH --job-name=v3_finalize
#SBATCH --partition=general-long
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/mnt/scratch/dasilvaf/iwgc_pipeline_v3/logs/%x_%j.out
# Final v3 release from attempt A (BRAKER3 + Iso-Seq):
#   PhaMi.v3.gff3 (lab-standard IDs, same scheme as v2), PhaMi.v3.proteins.fa,
#   PhaMi.v3.cds.fa, PhaMi.v3.longest.proteins.fa
# then BUSCO poales_odb10 (proteins) on v3 and on Daniel's v2 proteins, so
# both are measured the same way.
set -euo pipefail
source "${SLURM_SUBMIT_DIR}/config.sh"
REL="$FINAL/PhaMi_v3"
GTF="$WORK/braker_A_isoseq/braker.gtf"
GFFREAD_SIF=$V3/images/gffread.sif
V2_PROT=$V3/inputs/PhaMi.v2.proteins.fa
mkdir -p "$REL"
[ -s "$GTF" ] || { echo "ERROR: $GTF missing"; exit 1; }
[ -s "$V2_PROT" ] || { echo "ERROR: $V2_PROT missing (scp it from the laptop)"; exit 1; }

module purge
module load "$SAMTOOLS_MODULE"
[ -s "$GENOME.fai" ] || samtools faidx "$GENOME"

python3 "${SLURM_SUBMIT_DIR}/finalize_v3.py" "$GTF" "$GENOME.fai" PhaMi > "$REL/PhaMi.v3.gff3"

# -S: '*' for stop codons; -x CDS, -y protein (same order/IDs as the GFF3)
sing "$GFFREAD_SIF" gffread "$REL/PhaMi.v3.gff3" -g "$GENOME" -S \
  -x "$REL/PhaMi.v3.cds.fa" -y "$REL/PhaMi.v3.proteins.raw.fa"
# drop the terminal '*' and report internal stops (should be 0)
awk '/^>/{if(h){sub(/\*$/,"",s); print h; while(length(s)){print substr(s,1,60); s=substr(s,61)}} h=$1; s=""; next}
     {s=s $0}
     END{sub(/\*$/,"",s); print h; while(length(s)){print substr(s,1,60); s=substr(s,61)}}' \
  "$REL/PhaMi.v3.proteins.raw.fa" > "$REL/PhaMi.v3.proteins.fa"
rm "$REL/PhaMi.v3.proteins.raw.fa"
sing "$GFFREAD_SIF" gffread "$REL/PhaMi.v3.gff3" -g "$GENOME" -E -o /dev/null 2> "$REL/gffread_validation.log" || true

# longest isoform per gene (for BUSCO / one-protein-per-gene uses)
awk '/^>/{if(h)print h"\t"s; h=substr($1,2); s=""; next}{s=s $0}END{print h"\t"s}' "$REL/PhaMi.v3.proteins.fa" \
  | awk -F'\t' '{g=$1; sub(/-mRNA-[0-9]+$/,"",g); if(length($2)>L[g]){L[g]=length($2); S[g]=$2; I[g]=$1}}
                END{for(g in S) print I[g]"\t"S[g]}' | sort -k1,1 \
  | awk -F'\t' '{print ">"$1; s=$2; while(length(s)){print substr(s,1,60); s=substr(s,61)}}' \
  > "$REL/PhaMi.v3.longest.proteins.fa"

{
  echo "PhaMi v3 (BRAKER3 ETP: Iso-Seq FLNC + Lolium/OrthoDB proteins, one genome-wide model)"
  echo "genes: $(awk -F'\t' '$3=="gene"' "$REL/PhaMi.v3.gff3" | wc -l)"
  echo "mRNAs: $(awk -F'\t' '$3=="mRNA"' "$REL/PhaMi.v3.gff3" | wc -l)"
  echo "proteins: $(grep -c '^>' "$REL/PhaMi.v3.proteins.fa")  CDS: $(grep -c '^>' "$REL/PhaMi.v3.cds.fa")"
  echo "longest-isoform proteins: $(grep -c '^>' "$REL/PhaMi.v3.longest.proteins.fa")"
  echo "proteins with internal stop: $(awk '/^>/{next}/\*/{n++}END{print n+0}' "$REL/PhaMi.v3.proteins.fa")"
  echo "proteins not starting with M: $(awk '/^>/{getline; if(substr($0,1,1)!="M")n++}END{print n+0}' "$REL/PhaMi.v3.proteins.fa")"
  echo "CDS length not multiple of 3: $(awk '/^>/{if(h&&l%3)n++; h=1; l=0; next}{l+=length($0)}END{if(l%3)n++; print n+0}' "$REL/PhaMi.v3.cds.fa")"
  echo "first gene: $(awk -F'\t' '$3=="gene"{print $9; exit}' "$REL/PhaMi.v3.gff3")"
} | tee "$REL/release_stats.txt"

module purge
module load "$BUSCO_MODULE"
cd "$REL"
for pair in "v3:$REL/PhaMi.v3.longest.proteins.fa" "v2:$V2_PROT"; do
  name=${pair%%:*}; fa=${pair#*:}
  busco -i "$fa" -m proteins -l "$BUSCO_LINEAGE" --offline -c "$SLURM_CPUS_PER_TASK" -o "busco_poales_$name" -f
  cp busco_poales_$name/short_summary*.txt "$REL/busco_poales_$name.summary.txt"
  tar -czf "busco_poales_$name.tar.gz" "busco_poales_$name" && rm -rf "busco_poales_$name"
done
rm -rf busco_downloads
grep -H "C:" "$REL"/busco_poales_*.summary.txt
check_files
