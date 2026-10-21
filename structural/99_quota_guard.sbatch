#!/bin/bash --login
#SBATCH --job-name=v3_quotaguard
#SBATCH --partition=general-long
#SBATCH --time=7-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=/mnt/scratch/dasilvaf/iwgc_pipeline_v3/logs/%x_%j.out
# Scratch file-quota guard (1,048,576 files). GeneMark splits the soft-masked
# 4.9 Gb genome into ~178k pieces and writes output/data + output/gmhmm files
# per piece (~357k files per GeneMark run). Once a GeneMark step has written
# its genemark.gtf, those per-piece files are no longer used: delete them.
# Last resort: if usage still passes HARD_LIMIT, cancel attempt B (the
# protein-only control) so attempt A can finish.
source "${SLURM_SUBMIT_DIR}/config.sh"
HARD_LIMIT=980000

files_used() { /usr/lpp/mmfs/bin/mmlsquota -u "$USER" gs21 2>/dev/null | awk '$2=="scratch"{print $10}'; }

while squeue -u "$USER" -h -o "%j" | grep -q '^v3_braker_'; do
  # Only top-level GeneMark-ES / GeneMark-EP, and only once genemark.gtf has
  # been untouched for 15 min (gmes writes it last, from output/gmhmm).
  # GeneMark-ETP may reuse its internal dirs, so it's only counted, never cleaned.
  find "$WORK"/braker_*/GeneMark-ES "$WORK"/braker_*/GeneMark-EP -maxdepth 1 -type f \
       -name genemark.gtf -size +0 -mmin +15 2>/dev/null | while read -r gtf; do
    d=$(dirname "$gtf")
    for sub in output/data output/gmhmm; do
      if [ -d "$d/$sub" ]; then
        n=$(find "$d/$sub" | wc -l)
        rm -rf "${d:?}/$sub" && echo "$(date) removed $n files: $d/$sub"
      fi
    done
  done
  used=$(files_used)
  etp=$(find "$WORK"/braker_A_isoseq/GeneMark-ETP 2>/dev/null | wc -l)
  echo "$(date) scratch files: $used / $FILE_LIMIT (GeneMark-ETP dir: $etp)"
  if [ -n "$used" ] && [ "$used" -gt "$HARD_LIMIT" ]; then
    jid=$(squeue -u "$USER" -h -n v3_braker_B -o %i)
    if [ -n "$jid" ]; then
      echo "$(date) OVER $HARD_LIMIT files: cancelling braker B ($jid) to protect A"
      scancel "$jid"
    fi
  fi
  sleep 600
done
echo "$(date) no BRAKER jobs left, guard exiting"
