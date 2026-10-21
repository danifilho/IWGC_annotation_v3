#!/bin/bash --login
#SBATCH --job-name=v3_setup
#SBATCH --partition=general-short
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/mnt/scratch/dasilvaf/iwgc_pipeline_v3/logs/%x_%j.out
# Pulls the BRAKER3 Iso-Seq container, builds the full protein DB and one
# writable AUGUSTUS config per attempt.
set -euo pipefail
source "${SLURM_SUBMIT_DIR}/config.sh"
mkdir -p "$WORK/isoseq" "$FINAL" "$TMPDIR" "$APPTAINER_CACHEDIR"
check_files

if [ ! -s "$BRAKER_ISO_SIF" ]; then
  singularity pull "$BRAKER_ISO_SIF.partial" docker://teambraker/braker3:isoseq
  mv "$BRAKER_ISO_SIF.partial" "$BRAKER_ISO_SIF"
fi
sing "$BRAKER_ISO_SIF" braker.pl --version

# Lolium proteins + OrthoDB Viridiplantae. OrthoDB headers are already simple
# (e.g. >101020_0:000000); Lolium headers get cut at the first space.
if [ ! -s "$PROT" ]; then
  { cut -d' ' -f1 "$LOLRI_PROT"; cat "$ORTHODB_PROT"; } > "$PROT.tmp"
  mv "$PROT.tmp" "$PROT"
fi
echo "proteins in DB: $(grep -c '^>' "$PROT")"

for a in A_isoseq B_protein; do
  sif=$BRAKER_SIF; [ "$a" = A_isoseq ] && sif=$BRAKER_ISO_SIF
  if [ ! -d "$WORK/augustus_config_$a" ]; then
    sing "$sif" cp -r /opt/Augustus/config "$WORK/augustus_config_$a"
  fi
done
# the docker layer cache is only needed for the pull
rm -rf "$APPTAINER_CACHEDIR"/*
check_files
echo "setup done"
