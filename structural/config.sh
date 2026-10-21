#!/bin/bash
# Config for the Phalaris minor v3 annotation (BRAKER fixes over John's v2 run).
# All job scripts source this file via ${SLURM_SUBMIT_DIR}: submit from scripts/.
#
# What v3 changes vs. John's run (/mnt/research/Weed_lab/projects/Phalaris_minor):
#  1. ONE genome-wide model instead of 15 per-chromosome models (his AUGUSTUS
#     accuracy was 0.21-0.42 per chromosome).
#  2. Full protein DB: L. rigidum + OrthoDB Viridiplantae (his braker_prot.faa
#     only had the 54,112 Lolium proteins; ORTHODB_PROT never got included).
#  3. Attempt A uses the Iso-Seq reads INSIDE BRAKER3 (teambraker/braker3:isoseq,
#     GeneMark-ETP with long-read StringTie); his run used them only in TSEBRA.
#  Attempt B = same as A but protein-only (EP mode), as control/fallback.

V3=/mnt/scratch/dasilvaf/iwgc_pipeline_v3
SP=PhaMi
J=/mnt/research/Weed_lab/projects/Phalaris_minor   # John's project (read-only)

GENOME=$V3/inputs/PhaMi.softmasked.fa               # copy of John's EDTA/RepeatMasker soft-masked genome (real file: symlinks into /mnt/research do not resolve inside the container)
LOLRI_PROT=$J/inputs/Lolium_rigidum_APGP_CSIRO_protein.fa
ORTHODB_PROT=$J/ortho_db/Viridiplantae.fa
PROT=$V3/inputs/proteins_v3.faa                     # Lolium + OrthoDB Viridiplantae (made in 00)
FLNC_BAM=$V3/inputs/PhaMi.flnc.bam                  # isoseq refine output: FLNC reads (symlink)
ISO_BAM=$V3/work/isoseq/PhaMi.flnc.minimap2.bam     # made in 01
BUSCO_LINEAGE=$V3/inputs/poales_odb10               # symlink to John's local lineage
BUSCO_MODULE="BUSCO/5.8.2-foss-2023a"
MINIMAP2_MODULE="minimap2/2.26-GCCcore-12.3.0"
SAMTOOLS_MODULE="SAMtools/1.18-GCC-12.3.0"

BRAKER_SIF=$V3/images/braker3.sif                   # teambraker/braker3 (BRAKER 3.0.8), symlink to John's
BRAKER_ISO_SIF=$V3/images/braker3_isoseq.sif        # teambraker/braker3:isoseq (pulled in 00)

WORK=$V3/work
FINAL=$V3/final

# Keep everything off $HOME (home quota is full) and bind every path the
# symlinks resolve to.
export TMPDIR=$V3/tmp
export APPTAINER_CACHEDIR=$V3/cache SINGULARITY_CACHEDIR=$V3/cache
export APPTAINER_TMPDIR=$V3/tmp SINGULARITY_TMPDIR=$V3/tmp
BINDS=/mnt/scratch,/mnt/gs21,/mnt/research,/mnt/ufs18
sing() { local sif="$1"; shift; singularity exec --cleanenv -B "$BINDS" --env TMPDIR="$TMPDIR" "$sif" "$@"; }

# File-count guard: scratch limit is 1,048,576 files per user.
FILE_LIMIT=1048576
check_files() {
  local used
  used=$(/usr/lpp/mmfs/bin/mmlsquota -u "$USER" gs21 2>/dev/null | awk '$2=="scratch"{print $10}')
  echo "$(date) scratch files in use: ${used:-unknown} / $FILE_LIMIT"
}
