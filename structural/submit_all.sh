#!/bin/bash
# Run from the scripts/ directory: bash submit_all.sh
set -euo pipefail
cd "$(dirname "$0")"
s0=$(sbatch --parsable 00_setup.sbatch)
s1=$(sbatch --parsable --dependency=afterok:$s0 01_isoseq_align.sbatch)
sA=$(sbatch --parsable --dependency=afterok:$s1 --job-name=v3_braker_A 02_braker.sbatch A_isoseq)
sB=$(sbatch --parsable --dependency=afterok:$s0 --job-name=v3_braker_B 02_braker.sbatch B_protein)
eA=$(sbatch --parsable --dependency=afterok:$sA --job-name=v3_eval_A 03_eval.sbatch A_isoseq)
eB=$(sbatch --parsable --dependency=afterok:$sB --job-name=v3_eval_B 03_eval.sbatch B_protein)
echo "setup=$s0 isoalign=$s1 braker_A=$sA braker_B=$sB eval_A=$eA eval_B=$eB" | tee -a ../logs/submitted_jobs.txt
