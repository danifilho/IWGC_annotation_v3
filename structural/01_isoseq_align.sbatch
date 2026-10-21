#!/bin/bash --login
#SBATCH --job-name=v3_isoalign
#SBATCH --partition=general-long
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --output=/mnt/scratch/dasilvaf/iwgc_pipeline_v3/logs/%x_%j.out
# Aligns the 4.1M Iso-Seq FLNC READS (not collapsed transcripts) as the BRAKER
# README asks for the braker3:isoseq variant: StringTie needs each intron
# covered by many alignments, which collapsed transcripts can't provide.
# minimap2 -ax splice:hq -uf, exactly as in the README; -I 16G so the 4.9 Gb
# genome is indexed in one batch.
set -euo pipefail
source "${SLURM_SUBMIT_DIR}/config.sh"
module purge
module load "$MINIMAP2_MODULE"
module load "$SAMTOOLS_MODULE" 2>/dev/null || module load SAMtools
mkdir -p "$WORK/isoseq"
if [ -s "$ISO_BAM" ] && [ -s "$ISO_BAM.bai" ]; then echo "already done"; exit 0; fi

FA="$WORK/isoseq/PhaMi.flnc.fasta"
[ -s "$FA" ] || samtools fasta -@ 8 "$FLNC_BAM" > "$FA"
echo "FLNC reads: $(grep -c '^>' "$FA")"

minimap2 -t "$SLURM_CPUS_PER_TASK" -I 16G -ax splice:hq -uf "$GENOME" "$FA" \
  | samtools sort -@ 8 -m 2G -T "$TMPDIR/isosort" -o "$ISO_BAM.partial" -
mv "$ISO_BAM.partial" "$ISO_BAM"
samtools index "$ISO_BAM"   # .bai is fine: largest chromosome (515 Mb) < 2^29
samtools flagstat -@ 8 "$ISO_BAM" | tee "$WORK/isoseq/flagstat.txt"
rm -f "$FA"
echo "$(date) done"
