# Troubleshooting

Every entry below is an error that actually happened while building and running this pipeline on the MSU ICER HPCC (Phalaris minor, 4.9 Gb). Each one gives the exact message, the cause and the fix. Entries marked **(handled)** are already prevented by the current Snakefiles; they are kept here so the message can be recognized if it appears in an older setup.

**Quick checks first**

- Snakemake runs with `--keep-going`, so one failed rule does **not** stop the run. After every run, search the log for failed rules:
  ```bash
  grep -n "Error in rule" *.err
  ```
- Check the scratch quota for both space **and number of files**:
  ```bash
  /usr/lpp/mmfs/bin/mmlsquota -u $USER gs21    # scratch
  quota                                         # home
  ```

---

## Contents

**Containers, paths and quotas**
- [#1 `Genome file ... does not exist` inside a container](#1-genome-file--does-not-exist-inside-a-container)
- [#2 Scratch file quota explodes during BRAKER (GeneMark)](#2-scratch-file-quota-explodes-during-braker-genemark)
- [#3 Home quota full / `disk quota exceeded` while pulling images](#3-home-quota-full--disk-quota-exceeded-while-pulling-images)
- [#4 `singularity pull -q ...` prints the usage text](#4-singularity-pull--q--prints-the-usage-text)
- [#5 BUSCO prints a Python traceback from `BuscoLogger.py`](#5-busco-prints-a-python-traceback-from-buscologgerpy)

**Structural pipeline (BRAKER3)**
- [#6 `Only Protein input detected, BRAKER will be executed in EP mode (BRAKER2)`](#6-only-protein-input-detected-braker-will-be-executed-in-ep-mode-braker2)
- [#7 Very low AUGUSTUS accuracy (< 0.5) and bad BUSCO](#7-very-low-augustus-accuracy--05-and-bad-busco)
- [#8 The protein DB is much smaller than expected](#8-the-protein-db-is-much-smaller-than-expected)
- [#9 `The number of threads was set to 64, which is greater than 48`](#9-the-number-of-threads-was-set-to-64-which-is-greater-than-48)
- [#10 BRAKER rerun fails because the species already exists](#10-braker-rerun-fails-because-the-species-already-exists)
- [#11 `samtools index` fails on large chromosomes](#11-samtools-index-fails-on-large-chromosomes)
- [#12 minimap2 is slow, uses too much RAM or splits the index](#12-minimap2-is-slow-uses-too-much-ram-or-splits-the-index)
- [#13 Iso-Seq: collapsed transcripts instead of reads](#13-iso-seq-collapsed-transcripts-instead-of-reads)
- [#14 TSEBRA crashes with `invalid literal for int() ... '.'`](#14-tsebra-crashes-with-invalid-literal-for-int--)
- [#15 `.gm_key` / GeneMark license](#15-gm_key--genemark-license)
- [#25 Proteins with internal stop codons](#25-proteins-with-internal-stop-codons)

**QC (BUSCO)**
- [#16 BUSCO reports 0% (C:0%, F:0%, M:100%)](#16-busco-reports-0-c0-f0-m100)
- [#17 Very high BUSCO "Duplicated"](#17-very-high-busco-duplicated)

**Functional pipeline**
- [#18 `Could not create symlink of /tmp/...` in rule `mmseqs2`](#18-could-not-create-symlink-of-tmp-in-rule-mmseqs2)
- [#19 NCBI database FASTA with glued headers](#19-ncbi-database-fasta-with-glued-headers)
- [#28 Most NCBI best hits are the same ID (`AA049522-RA`)](#28-most-ncbi-best-hits-are-the-same-id-aa049522-ra)
- [#20 AGAT functional GFF rejects the InterProScan TSV](#20-agat-functional-gff-rejects-the-interproscan-tsv)
- [#27 AGAT functional GFF: every gene is "hypothetical protein"](#27-agat-functional-gff-every-gene-is-hypothetical-protein)
- [#21 InterProScan runs out of memory](#21-interproscan-runs-out-of-memory)
- [#22 `Failed to solve scheduling problem with ILP solver`](#22-failed-to-solve-scheduling-problem-with-ilp-solver)
- [#23 `Directory cannot be locked`](#23-directory-cannot-be-locked)
- [#26 Snakemake reruns BRAKER after a small edit](#26-snakemake-reruns-braker-after-a-small-edit)
- [#24 Fewer proteins in the functional TSV than genes in the GFF](#24-fewer-proteins-in-the-functional-tsv-than-genes-in-the-gff)

---

## Containers, paths and quotas

### 1. `Genome file ... does not exist` inside a container

```
ERROR: in file /opt/BRAKER/scripts/braker.pl at line 3662
Genome file /mnt/scratch/.../inputs/PhaMi.softmasked.fa does not exist.
```

**Cause.** The file is a symlink that points outside the bind-mounted volume. In this case it pointed to `/mnt/research/...`, an automounted path whose real location is `/mnt/ufs18/...`. The host resolves the symlink, but inside the container the target is not mounted. Binding `/mnt/research` does not help, because automounts are not visible through the bind.

**Fix.** Every input must be a real file inside `volume_name`:
- On the same filesystem, use a hard link: `ln /path/genome.fa inputs/` (no extra space or files).
- Across filesystems, use a copy: `cp -L`.

The Snakefiles bind only `volume_name`.

### 2. Scratch file quota explodes during BRAKER (GeneMark)

Symptom: within a few hours of starting BRAKER, `mmlsquota` jumps from about 140 k to more than 850 k files (the limit is 1,048,576). Past the limit, jobs fail with `Disk quota exceeded`.

**Cause.** GeneMark splits the soft-masked genome at masked stretches. A 4.9 Gb grass genome gives about **178,000 pieces**, and every GeneMark run writes one `output/data` and one `output/gmhmm` file per piece, which is about 357 k files per run. BRAKER2 runs both GeneMark-ES and GeneMark-EP.

**Fix (handled).** The `braker` rule starts `scripts/genemark_cleanup.sh` in the background. Every 10 minutes it deletes `output/data` and `output/gmhmm` from `GeneMark-ES` and `GeneMark-EP` once their `genemark.gtf` has been unchanged for 15 minutes. `GeneMark-ETP` is never touched. Progress is logged in `logs/genemark_cleanup.log`.

If file usage is still tight:
- avoid running two BRAKER jobs at the same time;
- ask ICER for a temporary increase of the file quota.

### 3. Home quota full / `disk quota exceeded` while pulling images

**Cause.** `singularity pull` and `apptainer` keep their cache in `~/.singularity` and `~/.apptainer`. These caches had reached about 16 GB, which filled the 110 GB home quota.

**Fix.** Point the cache and temporary directories to scratch before pulling or running. The `run_*.sb` scripts already do this:
```bash
export APPTAINER_CACHEDIR=$PWD/tmp SINGULARITY_CACHEDIR=$PWD/tmp TMPDIR=$PWD/tmp
```
To free the space already used:
```bash
singularity cache clean
```

### 4. `singularity pull -q ...` prints the usage text

**Cause.** `-q` is a global flag (`singularity -q pull ...`), not an option of `pull`.

**Fix.**
```bash
singularity pull busco.sif docker://ezlabgva/busco:v5.8.2_cv1
```

### 5. BUSCO prints a Python traceback from `BuscoLogger.py`

```
Traceback (most recent call last):
  File ".../busco/BuscoLogger.py", line 186, in __init__
```

**Cause.** BUSCO writes a log file to the current directory. When that directory is read-only (for example, running `busco --version` from inside `images/`), the logger fails, even though the version still prints.

**Fix.** Run BUSCO from a writable directory. The `busco` rule changes into `busco_outputs/` first.

---

## Structural pipeline (BRAKER3)

### 6. `Only Protein input detected, BRAKER will be executed in EP mode (BRAKER2)`

**Cause.** BRAKER received only `--prot_seq`, so it fell back to protein-only BRAKER2. Transcript evidence passed only to TSEBRA afterwards **does not help**: TSEBRA can only choose among existing gene models. This is exactly what produced C:67.7 % / M:30.1 % BUSCO in the earlier *P. minor* run.

**Fix.**
- Set `isoseq_reads` in `config.yaml`. The pipeline then uses the `braker3:isoseq` container with `--bam`, and the log must say:
  ```
  Both protein and RNA-Seq data in input detected. BRAKER will be executed in ETP mode (BRAKER3).
  ```
- Check this line in `logs/braker.log` or `braker_outputs/braker.log` early in the run.

### 7. Very low AUGUSTUS accuracy (< 0.5) and bad BUSCO

`braker.log` lines to read:
```
The accuracy after initial training (no optimize_augustus.pl, no CRF) is 0.24
The accuracy after training (after optimize_augustus.pl, no CRF) is 0.27
```

**Cause.** A small training set. In the earlier *P. minor* run, BRAKER was run **separately on each chromosome**, giving 15 independent models trained on a few hundred genes each (accuracy 0.21–0.42). The resulting gene set had 62 % single-exon transcripts, a median protein of 221 aa, and some "genes" spanning about 500 kb.

**Fix.** Run **one genome-wide BRAKER**, as this pipeline does (accuracy 0.75). If runtime is a problem, train once and predict per chromosome with `--skipAllTraining --species=<trained model>`; do not train per chromosome.

### 8. The protein DB is much smaller than expected

Symptom: `proteins in DB:` (log of rule `protein_db`) shows only the close relative, for example 54,112 instead of about 5.4 M.

**Cause.** The old SLURM script built the DB with `cat $LOLIUM ${ORTHODB:+"$ORTHODB"}`. This form silently skips OrthoDB whenever that variable is empty or unset, and the prep log of the earlier run (`54112 proteins for BRAKER`) shows that this is what happened.

**Fix.** List **every** file explicitly under `protein_files` in `config.yaml`, and check the count in the log. BRAKER's ProtHint needs a broad set such as OrthoDB Viridiplantae (about 5.3 M proteins), not only one species.

### 9. `The number of threads was set to 64, which is greater than 48`

```
WARNING: ... GeneMark has in the past been reported to die if you set such a high number of threads.
```

**Fix.** Keep `braker_threads: 48` or lower; the default config already does.

### 10. BRAKER rerun fails because the species already exists

**Cause.** BRAKER cannot resume an interrupted run, and a half-trained AUGUSTUS species is left behind in `augustus_config/species/`.

**Fix (handled).** The `braker` rule deletes `braker_outputs/` and `augustus_config/species/<braker_species>` before starting. To keep a trained model, copy `braker_inputs/augustus_config/species/<braker_species>` somewhere safe after a successful run.

### 11. `samtools index` fails on large chromosomes

```
[E::hts_idx_push] Region ... cannot be stored in a bai index. Try using a csi index
```

**Cause.** BAI indexes only support sequences shorter than 2^29 bp (about 537 Mb). Large grass chromosomes can be bigger (Chr01C of *P. minor* is 515 Mb, just under the limit).

**Fix (handled).** The `minimap2` rule falls back to `samtools index -c` (CSI).

### 12. minimap2 is slow, uses too much RAM or splits the index

**Cause.** By default minimap2 splits large genomes into several index batches.

**Fix.** Set `minimap2_index_batch` larger than the genome (`16G` for 4.9 Gb). On *P. minor*, 4.1 M FLNC reads aligned in **10 minutes** (32 threads, 31 GB peak RAM).

### 13. Iso-Seq: collapsed transcripts instead of reads

**Cause.** BRAKER's long-read mode assembles transcripts with StringTie inside GeneMark-ETP. Training "requires that each intron is covered by many alignments", which **collapsed or clustered transcripts cannot provide**.

**Fix.** Give `isoseq_reads` the **FLNC reads** (for example `isoseq refine` → `*.flnc.bam`) or HiFi reads, not `isoseq collapse` output. Read depth matters: the BRAKER authors report lower accuracy with shallow or high-error long-read data.

### 14. TSEBRA crashes with `invalid literal for int() ... '.'`

The crash is preceded by thousands of `no CDS nor exons` warnings.

**Cause.** Collapsed Iso-Seq GTFs passed to TSEBRA with `-k` (forced inclusion). Many near-duplicate single-exon isoforms break TSEBRA's overlap graph.

**Fix.** Don't post-filter with TSEBRA + Iso-Seq at all. Give the reads to BRAKER3 (see 6 and 13); BRAKER runs TSEBRA internally with the right inputs.

### 15. `.gm_key` / GeneMark license

The GeneMark bundled in the braker3 containers only **tests** for `~/.gm_key` and runs without it. No license key is needed for this pipeline.

### 25. Proteins with internal stop codons

Symptom: `proteins with internal stop:` is above 0 in `final_outputs/<SP>.v3.stats.txt`, or `*` characters appear inside protein sequences.

**Cause.** A few **GeneMark.hmm3** models in `braker.gtf` are translated with in-frame stops. The same stops are in BRAKER's own `braker.aa`, so the conversion to GFF3 is not the cause. On a single-chromosome test (Chr07B), 13 of 5,357 transcripts were affected; on the full *P. minor* genome, 1 of 49,216.

**Fix (handled).** With `drop_internal_stop_models: true` (the default), `finalize_gff3` drops these transcripts **before** numbering the IDs, so the IDs stay contiguous. The dropped BRAKER IDs are listed in `final_outputs/<SP>.v3.excluded_internal_stop.txt`. The functional pipeline would drop them anyway (`gffread -J`).

> Older versions of `gene_set_stats.sh` counted sequence **lines** containing `*` instead of proteins, which overstated this number (310 instead of 13).

---

## QC (BUSCO)

### 16. BUSCO reports 0% (C:0%, F:0%, M:100%)

`logs/hmmsearch_err.log`:
```
Fatal exception (source file p7_pipeline.c, line 714):
Target sequence length > 100K, over comparison pipeline limit.
```

**Cause.** Not a biological result: `hmmsearch` aborts on any protein longer than 100,000 aa. Poorly trained models (see 7) produced proteins of up to 142 k aa.

**Fix (handled).** The `busco` rule drops proteins longer than `busco_max_protein_length` (10,000 aa) and runs on the longest isoform per gene. If such giant proteins appear at all, the gene models are bad; check the AUGUSTUS accuracy.

### 17. Very high BUSCO "Duplicated"

For allopolyploids (*P. minor*: B and C subgenomes), D around 75–80 % is **expected**, because both homeologs are annotated. Compare **C** and **M** between annotations, not S.

---

## Functional pipeline

### 18. `Could not create symlink of /tmp/...` in rule `mmseqs2`

```
Could not create symlink of /tmp/11310392332275771123!
Error in rule mmseqs2: ... databases/current/GenomesDatabase.v1.prot
```

**Cause.** The three `mmseqs2` jobs (ncbi, herbres, uniref) run **in parallel**, and the old rule gave all of them `/tmp` as the MMseqs tmp directory. MMseqs creates a `latest` symlink inside its tmp directory; the runs raced, and one search, here the NCBI one, died after about 1 minute. Because of `--keep-going`, everything else kept running for about a day before the missing NCBI result was noticed.

**Fix (handled).** Each database now gets its own `tmp/mmseqs_<database>` inside the work volume.

If a workdir was started with the old Snakefile, rerun only that search, then rerun the pipeline (Snakemake picks up the finished `.ofmt6`):
```bash
singularity exec -B $PWD:/data images/mmseqs2.sif bash -c "cd /data && mmseqs easy-search \
  gff_read_output/SP.aa databases/current/GenomesDatabase.v1.prot \
  mmseqs_output/SP_ncbi.ofmt6 /data/tmp/mmseqs_ncbi --threads 20 --format-mode 0 \
  --start-sens 2 -s 7 --sens-steps 3"
```

### 19. NCBI database FASTA with glued headers

Symptoms: tools cannot find some subjects, or a "protein" contains a header in its sequence. The worst consequence is described in #28.

**Cause.** In `GenomesDatabase.v1.prot`, some headers are glued to the end of the previous sequence line (no newline before `>`), and the first header has no `>`.

**Fix (handled).** Repair the FASTA once with `scripts/fix_fasta_headers.py`, then point the `ncbi` entry in `config.yaml` (both `path` and `ncbidata`) at the repaired file:
```bash
python3 scripts/fix_fasta_headers.py databases/current/GenomesDatabase.v1.prot databases/current/GenomesDatabase.v1.fixed.prot
# records: 5989485  glued headers repaired: 49541  missing first '>' added: 1
```
The `mmseqs2` rule now refuses to search a FASTA target that has glued headers or does not start with `>`.

### 20. AGAT functional GFF rejects the InterProScan TSV

**Cause.** `--pathways` adds a 15th column that `agat_sp_manage_functional_annotation.pl` does not accept.

**Fix.** Pass `cut -f1-14` of the InterProScan TSV. This only matters if you run AGAT's functional script yourself; the pipeline no longer uses it (see #27).

### 21. InterProScan runs out of memory

Symptom: crash in the LevelDB / match-lookup stage on large proteomes.

**Fix (handled).** The `interproscan` rule runs in chunks of 2,000 proteins and skips chunks that are already done, so it can resume. Each chunk takes about 1.5 h (30 CPU), so 45 k proteins take about 34 h. The repeated `As the data in these versions is not the same, you cannot use this match lookup service` messages are only a warning: the local analyses still run.

### 22. `Failed to solve scheduling problem with ILP solver`

```
Failed to solve scheduling problem with ILP solver, falling back to greedy scheduler.
```

This is harmless: Snakemake falls back to the greedy scheduler. To silence it, add `--scheduler greedy`.

### 23. `Directory cannot be locked`

**Cause.** A previous Snakemake run was killed and left its lock behind.

**Fix.** The `run_*.sb` scripts call `snakemake --unlock` first. Never run two Snakemake processes on the same workdir.

### 26. Snakemake reruns BRAKER after a small edit

Symptom: after editing a script, the config or the Snakefile, a rerun starts again from `isoseq_fasta`/`minimap2`/`braker`, and BRAKER (hours) is redone. Since the `braker` rule clears `braker_outputs/` when it starts, killing the job at that point also loses the previous BRAKER result.

**Cause.** By default Snakemake also reruns jobs whose **code, parameters or software environment** changed, not only those whose inputs are newer. `isoseq_fasta` is a `temp()` output, so it is gone after the first run; any change upstream of it brings back the whole chain.

**Fix.**
- Always run with `--rerun-triggers mtime`, as `run_structural.sb` and `run_functional.sb` do.
- To redo only the post-BRAKER steps, force them explicitly:
  ```bash
  snakemake --cores 16 --use-singularity --singularity-args "-B $PWD" --rerun-triggers mtime -R finalize_gff3
  ```
- Always check what will run with `snakemake -n --rerun-triggers mtime` first.
- Keep a copy of `braker_outputs/` and of the trained species before experimenting.

### 27. AGAT functional GFF: every gene is "hypothetical protein"

Symptom: in a GFF made with `agat_sp_manage_functional_annotation.pl -b <ncbi best hits> --db <NCBI FASTA>`, every mRNA has `product=hypothetical protein`, and the AGAT log repeats:
```
Use of uninitialized value $desc in pattern match (m//) at .../agat_sp_manage_functional_annotation.pl line 259.
```

**Cause.** AGAT takes product names from **UniProt-style** headers (`>sp|P12345|NAME_SPECIES Description OS=... GN=...`). The NCBI headers (`>XP_020175095.1 description [species]`) do not match, so AGAT falls back to "hypothetical protein" for all genes. This silently affected the MAKER v2 functional GFF too (57,163 of 57,163 genes).

**Fix (handled).** The `functional_gff` rule (`iwgc_functional_gff.py`) writes the functional GFF3 in the IWGC release layout straight from the merged TSV and InterProScan, without AGAT. Check the result with:
```bash
grep -P "\tgene\t" <SP>.functional.gff3 | grep -c "Note=Unknown Protein"   # should be a minority
```

### 28. Most NCBI best hits are the same ID (`AA049522-RA`)

Symptoms:
- Thousands of proteins have the same NCBI best hit, `AA049522-RA`, with `Description not found` (TSV) or `Note=Unknown Protein:(...)` / `Dbxref=AA049522-RA:...` (GFF3).
- Subject coordinates in `mmseqs_output/<SP>_ncbi.ofmt6` go up to the millions (e.g. `tstart=1071720`).

**Cause.** This is the glued-header defect of #19. MMseqs only starts a new record at a `>` at the beginning of a line, so in the original `GenomesDatabase.v1.prot` everything after `>AA049522-RA` merged into **one 18,868,593 aa record**. That record holds the 49,542 `AA…` proteins, a whole proteome (described as *Avena sativa* in the IWGC v1 release `SalTr_v01.0.gff`). For a grass genome, that block contains the closest homolog of most genes, so it became the best hit of **26,230 of 44,960 *P. minor* v3 proteins (58 %)**, and of 28,207 in MAKER v2. It hides the real best hit. The IWGC v1 database already had this defect: all 21 `AA` hits in `SalTr_v01.0.gff` are `AA049522-RA`.

**Fix (handled).** Use the repaired FASTA (#19) and rerun the `ncbi` search: delete `mmseqs_output/<SP>_ncbi.ofmt6` and `<SP>_ncbi_best_hit.out`, then rerun the pipeline. Only the NCBI column, the TSV and the functional GFF3 change; the UniRef50, InterPro, Pfam, SignalP/TargetP and MultiLoc2 results are not affected.

After the repair, the 49,542 `AA…` proteins are separate records, but their headers have no name. For grass genomes they still win as the plain best hit (27,214 / 44,960 *P. minor* proteins). Use `ncbi_prefer_named_hits: true`, the default, to report the best **named** hit instead; see [FUNCTIONAL.md §3](FUNCTIONAL.md#3-mmseqs2--best_hits-homology). To check a finished run:
```bash
cut -f2 mmseqs_output/<SP>_ncbi_best_hit.out | sort | uniq -c | sort -rn | head -3   # no single subject should dominate
```

### 24. Fewer proteins in the functional TSV than genes in the GFF

**Cause.** This is expected. The functional pipeline:
1. keeps the longest isoform per gene (AGAT);
2. drops tRNAs and incomplete models with `gffread -J` (no start or stop codon, or an internal stop).

*P. minor* v3: 44,965 genes → 44,960 proteins. MAKER v2: 57,163 → 47,494, because MAKER leaves about 7,600 partial models.
