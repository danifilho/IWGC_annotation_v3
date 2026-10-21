#!/bin/bash --login
#SBATCH --job-name=v3_braker
#SBATCH --partition=general-long
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=400G
#SBATCH --output=/mnt/scratch/dasilvaf/iwgc_pipeline_v3/logs/%x_%j.out
# ONE genome-wide BRAKER run (one GeneMark + one AUGUSTUS model for the whole
# genome). Usage: sbatch 02_braker.sbatch A_isoseq | B_protein
#   A_isoseq  -> braker3:isoseq container, proteins + Iso-Seq FLNC BAM (ETP mode)
#   B_protein -> braker3 container, proteins only (EP mode), control/fallback
set -euo pipefail
source "${SLURM_SUBMIT_DIR}/config.sh"
A="${1:?usage: sbatch 02_braker.sbatch A_isoseq|B_protein}"
RUN="$WORK/braker_$A"
ACFG="$WORK/augustus_config_$A"

case "$A" in
  A_isoseq)  SIF=$BRAKER_ISO_SIF; EXTRA=(--bam="$ISO_BAM")
             [ -s "$ISO_BAM" ] || { echo "ERROR: $ISO_BAM missing (run 01)"; exit 1; } ;;
  B_protein) SIF=$BRAKER_SIF; EXTRA=() ;;
  *) echo "unknown attempt $A"; exit 1 ;;
esac
for f in "$SIF" "$GENOME" "$PROT"; do [ -s "$f" ] || { echo "ERROR: $f missing (run 00)"; exit 1; }; done
[ -d "$ACFG" ] || { echo "ERROR: $ACFG missing (run 00)"; exit 1; }
if [ -s "$RUN/braker.gtf" ]; then echo "$A already finished"; exit 0; fi
# BRAKER can't resume; a half-finished run would also leave a trained species
# behind, so start clean.
rm -rf "$RUN" "$ACFG/species/${SP}_v3_$A"
mkdir -p "$RUN"
check_files

echo "$(date) BRAKER $A start"
sing "$SIF" braker.pl \
  --genome="$GENOME" \
  --prot_seq="$PROT" \
  "${EXTRA[@]}" \
  --species="${SP}_v3_$A" \
  --workingdir="$RUN" \
  --AUGUSTUS_CONFIG_PATH="$ACFG" \
  --threads=48 \
  --gff3
echo "$(date) BRAKER $A finished"
grep -iE "mode|accuracy" "$RUN/braker.log" || true
echo "files in $RUN: $(find "$RUN" | wc -l)"
check_files
