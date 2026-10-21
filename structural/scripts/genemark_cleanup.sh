#!/bin/bash
# Background cleaner started by the Snakefile's braker rule.
# GeneMark splits a soft-masked genome at masked stretches; on a 4.9 Gb grass
# genome that is ~178k pieces, with one output/data and one output/gmhmm file
# each (~357k files per GeneMark run). Once a GeneMark-ES / GeneMark-EP run has
# written its genemark.gtf and left it untouched for 15 min, those per-piece
# files are no longer used and are removed. GeneMark-ETP is never touched.
# Usage: genemark_cleanup.sh <braker_workingdir>   (killed by the braker rule)
RUN="$1"
while true; do
  for d in "$RUN"/GeneMark-ES "$RUN"/GeneMark-EP; do
    gtf="$d/genemark.gtf"
    if [ -s "$gtf" ] && [ -n "$(find "$gtf" -mmin +15 2>/dev/null)" ]; then
      for sub in output/data output/gmhmm; do
        if [ -d "$d/$sub" ]; then
          n=$(find "$d/$sub" | wc -l)
          rm -rf "${d:?}/$sub" && echo "$(date) removed $n files: $d/$sub"
        fi
      done
    fi
  done
  echo "$(date) files under $RUN: $(find "$RUN" 2>/dev/null | wc -l)"
  sleep 600
done
