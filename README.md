# IWGC Annotation v3: BRAKER3 + Iso-Seq

Structural and functional genome annotation for weed genomes in the
International Weed Genomics Consortium (IWGC), first applied to
***Phalaris minor*** (`PhaMi`, allotetraploid, 4.92 Gb, 15 sequences).

This repository contains three things:

1. A **diagnosis** of why an earlier BRAKER run on *P. minor* produced very poor BUSCO scores.
2. A corrected **structural pipeline** (BRAKER3 in ETP mode, trained once on the whole genome, using Iso-Seq reads and a large protein database). It reaches the same BUSCO completeness as the MAKER-based v2 annotation in about **1/10 of the wall-clock time**.
3. The **functional pipeline** (IWGC functional Snakefile) with a fix for a race condition that silently dropped the NCBI homology search.

> **Status:** structural annotation finished and validated. Functional annotation of v3 is still running on MSU HPCC (InterProScan stage); this README will be updated with its results.

---

## Contents

- [Results at a glance](#results-at-a-glance)
- [Why the previous BRAKER run failed](#why-the-previous-braker-run-failed)
- [What v3 changes](#what-v3-changes)
- [Repository layout](#repository-layout)
- [Structural pipeline](#structural-pipeline)
- [Functional pipeline](#functional-pipeline)
- [Output: PhaMi v3 annotation](#output-phami-v3-annotation)
- [Runtime: BRAKER3 vs MAKER](#runtime-braker3-vs-maker)
- [HPCC notes: quotas and containers](#hpcc-notes-quotas-and-containers)
- [Known issues and open questions](#known-issues-and-open-questions)
- [Requirements](#requirements)
- [Contact](#contact)

---

## Results at a glance

BUSCO v5.8.2, `poales_odb10` (n = 4,896), protein mode.

| | Previous BRAKER run | v3 control | **v3 final** | MAKER v2 (reference) |
|---|---|---|---|---|
| Method | BRAKER2, 15 per-chromosome models, *Lolium* proteins only, Iso-Seq only in TSEBRA | BRAKER2, **1 genome-wide model**, Lolium + OrthoDB proteins | **BRAKER3 (ETP), 1 genome-wide model, Lolium + OrthoDB proteins + Iso-Seq FLNC reads** | MAKER est2genome + protein2genome |
| AUGUSTUS accuracy | 0.21 – 0.42 | 0.76 | **0.75** | – |
| Genes | ~55.6 k | 96.8 k | **45.0 k** | 55.1 k protein-coding (+2,039 tRNA) |
| Transcripts | 57.6 k | 102.2 k | **49.2 k** | 55.1 k |
| Single-exon transcripts | 62 % | 38 % | **25 %** | – |
| CDS exons per transcript | 2.7 | 3.2 | **5.3** | – |
| Median protein length (aa) | 221 | 219 | **353** | – |
| **BUSCO C** [S, D] | 67.7 % [23.0, 44.7] | 71.5 % [16.8, 54.7] | **95.9 % [18.6, 77.3]** | 96.1 % [17.6, 78.4] |
| BUSCO F | 2.3 % | 5.4 % | **0.2 %** | 1.5 % |
| **BUSCO M** | 30.1 % | 23.1 % | **3.9 %** | 2.5 % |
| Wall-clock (gene prediction) | – | 16 h | **11 h** | ~5 days |

How to read these numbers:

- **Duplicated BUSCOs are expected.** *P. minor* is an allotetraploid with B and C subgenomes, so most conserved genes are annotated once in each homeolog.
- **The previous run's BUSCO** was measured on all isoforms. The v3 and MAKER numbers use the longest isoform per gene, which is the stricter measure.
- **The 0 % BUSCO** reported for the previous run's merged `braker_all.aa` was not a biological result. `hmmsearch` aborts on proteins longer than 100,000 aa (`Target sequence length > 100K`), and that file contained 216 proteins longer than 20,000 aa.

**Bottom line:** v3 fixes the BRAKER run (C 67.7 % → 95.9 %, M 30.1 % → 3.9 %). It matches the MAKER v2 annotation in completeness, with about 10 k fewer genes and roughly 10× less wall-clock time. The Iso-Seq data was sufficient all along; it simply never reached BRAKER's training step.

---

## Why the previous BRAKER run failed

All points below were checked directly in the previous run's `braker.log` files, scripts and outputs.

1. **BRAKER ran as BRAKER2, not BRAKER3.** The only evidence given to `braker.pl` was `--prot_seq`, and every log says:
   ```
   Only Protein input detected, BRAKER will be executed in EP mode (BRAKER2).
   ```
   The Iso-Seq data was added afterwards as intron hints to **TSEBRA**. TSEBRA only *selects* among gene models that already exist; it cannot create or correct them. So the Iso-Seq never influenced GeneMark or AUGUSTUS training.

2. **The protein database was only *Lolium rigidum*.** `braker_prot.faa` contained 54,112 proteins. The OrthoDB Viridiplantae set (5.3 M proteins) was configured but never concatenated in (`prep: 54112 proteins for BRAKER`). BRAKER2 depends heavily on a broad protein set.

3. **Fifteen independent models, one per chromosome.** Each chromosome trained its own GeneMark-EP and AUGUSTUS on a small gene set. AUGUSTUS accuracy after training was 0.21–0.42; the C-subgenome chromosomes were worst (Chr06C 0.21), and Chr00 had to reuse Chr01C's model. The resulting models show the damage: 62 % single-exon transcripts, a median protein of 221 aa, and "genes" spanning about 500 kb.

4. **The 0 % BUSCO was a crash**, as described in [Results at a glance](#results-at-a-glance).

---

## What v3 changes

| Problem | Fix in v3 |
|---|---|
| Iso-Seq used only in TSEBRA | The 4.1 M Iso-Seq **FLNC reads** are spliced-aligned with minimap2 (`-ax splice:hq -uf`, as the BRAKER docs require) and passed with `--bam` to the `teambraker/braker3:isoseq` container, so BRAKER runs in **ETP mode (BRAKER3)**. The long-read StringTie inside GeneMark-ETP needs many reads per intron, so raw reads are used instead of collapsed transcripts. |
| Lolium-only protein DB | *L. rigidum* + **OrthoDB Viridiplantae**: 5,364,589 proteins. |
| 15 per-chromosome models | **One genome-wide BRAKER run**: a single GeneMark-ETP model and a single AUGUSTUS model (`PhaMi_v3_A_isoseq`). |
| BUSCO crash on giant proteins | Evaluation drops proteins > 10,000 aa and uses the longest isoform per gene. |
| Containers could not read inputs | Inputs symlinked into `/mnt/research` do not resolve inside Singularity, so the genome is copied as a real file and all mount points are bound. |
| GeneMark file explosion (see [quotas](#hpcc-notes-quotas-and-containers)) | `99_quota_guard.sbatch` deletes per-chunk GeneMark files once each GeneMark step has finished. |
| Output not in the lab format | `finalize_v3.py` writes a clean GFF3 with the lab ID scheme, plus protein/CDS FASTA. |

A protein-only **control** (attempt B) uses the same proteins and the same single-model setup, without Iso-Seq. Comparing it with the previous run and with v3 separates the two effects:

- **Single model + more proteins** (previous run → B): AUGUSTUS accuracy 0.4 → 0.76, but BUSCO only 67.7 % → 71.5 %. Because both changes were made together, their individual contributions cannot be separated.
- **Adding Iso-Seq to training** (B → A): BUSCO 71.5 % → 95.9 %, and every gene-structure metric becomes plant-like.

---

## Repository layout

```
IWGC_annotation_v3/
├── structural/                  # BRAKER3 pipeline (SLURM)
│   ├── config.sh                # all paths, containers, modules; edit this first
│   ├── 00_setup.sbatch          # pull braker3:isoseq, build protein DB, AUGUSTUS configs
│   ├── 01_isoseq_align.sbatch   # FLNC BAM -> FASTA -> minimap2 splice:hq -> sorted BAM
│   ├── 02_braker.sbatch         # one genome-wide BRAKER run (A_isoseq | B_protein)
│   ├── 03_eval.sbatch           # gene-set stats + BUSCO (longest isoform, >10k aa dropped)
│   ├── 04_finalize.sbatch       # lab-format GFF3, proteins, CDS, BUSCO v3 vs v2
│   ├── 99_quota_guard.sbatch    # file-quota guard for GeneMark per-chunk files
│   ├── finalize_v3.py           # braker.gtf -> GFF3 with PhaMi IDs
│   └── submit_all.sh            # submits 00 -> 01 -> 02 -> 03 with dependencies
├── functional/                  # IWGC functional annotation (Snakemake)
│   ├── Snakefile                # IWGC functional Snakefile + mmseqs tmp fix
│   ├── config.yaml              # v3 config (GFF = PhaMi.v3.gff3)
│   ├── run_functional_v3.sb     # Snakemake + TSV in assembly order + AGAT functional GFF
│   ├── fix_mmseqs_ncbi.sb       # one-off rescue of the NCBI search (old Snakefile only)
│   └── scripts/                 # functional_merge.py, ml2_one_by_one.sh
└── results/
    ├── PhaMi_v3/                # PhaMi.v3.gff3.gz, release stats, BUSCO v3 & MAKER v2
    └── evaluation/              # stats + BUSCO for attempt A, control B, previous run
```

---

## Structural pipeline

### Inputs

| Input | Description |
|---|---|
| Soft-masked genome | EDTA library + RepeatMasker, 87–92 % masked per chromosome. Reused from the previous run. |
| Protein evidence | *Lolium rigidum* (APGP/CSIRO) + OrthoDB Viridiplantae |
| Iso-Seq | PacBio FLNC BAM (`isoseq3 refine` output), 4,116,976 reads |
| BUSCO lineage | `poales_odb10`, local copy (compute nodes run BUSCO `--offline`) |

### Steps

| Step | Script | Resources | Time on *P. minor* |
|---|---|---|---|
| Setup | `00_setup.sbatch` | 4 CPU / 32 GB | 18 min (mostly the container pull) |
| Iso-Seq alignment | `01_isoseq_align.sbatch` | 32 CPU / 120 GB | 10 min (31 GB peak RAM; 99.84 % of reads mapped) |
| BRAKER3 | `02_braker.sbatch A_isoseq` | 48 CPU / 400 GB | **10 h 50 min** |
| BRAKER2 control | `02_braker.sbatch B_protein` | 48 CPU / 400 GB | 15 h 56 min |
| Evaluation | `03_eval.sbatch` | 16 CPU / 64 GB | 7–15 min |
| Release + BUSCO | `04_finalize.sbatch` | 16 CPU / 64 GB | ~20 min |

The core command (attempt A) is:

```bash
singularity exec braker3_isoseq.sif braker.pl \
  --genome=PhaMi.softmasked.fa \
  --prot_seq=proteins_v3.faa \
  --bam=PhaMi.flnc.minimap2.bam \
  --species=PhaMi_v3_A_isoseq \
  --threads=48 --gff3
```

BRAKER itself warns that more than 48 threads can crash GeneMark.

### Running it

```bash
cd structural/
# 1. edit config.sh (V3, J, input paths, modules)
# 2. submit everything with dependencies
bash submit_all.sh
# 3. in parallel, start the file-quota guard (strongly recommended)
sbatch 99_quota_guard.sbatch
# 4. when attempt A is done
sbatch 04_finalize.sbatch
```

All jobs source `config.sh` from `$SLURM_SUBMIT_DIR`, so submit from `structural/`.

### Lab-format GFF3 (`finalize_v3.py`)

`braker.gtf` mixes AUGUSTUS, GeneMark.hmm3 and GeneMarkS-T records. `finalize_v3.py`:

- rebuilds each transcript from its `CDS` and `stop_codon` lines. BRAKER predicts no UTRs, so exons equal CDS. Stop codons are folded into the CDS where a source leaves them out.
- sorts genes in genome FASTA order (Chr01C…Chr07C, Chr01B…Chr07B, Chr00), then by start position.
- assigns **the same ID scheme as the lab's v2 annotation**, with the counter increasing by 10 across the whole genome and never reset per chromosome:
  ```
  gene  PhaMiChr01Cg000010
  mRNA  PhaMiChr01Cg000010-mRNA-1      (keeps braker_id=g4474.t1 for traceability)
  exon  PhaMiChr01Cg000010-mRNA-1:1
  CDS   PhaMiChr01Cg000010-mRNA-1:cds
  ```
- computes CDS phase from cumulative coding length.

`04_finalize.sbatch` then uses gffread to extract proteins (terminal `*` removed) and CDS, writes a longest-isoform protein set and sanity checks, and runs BUSCO on both v3 and the MAKER v2 proteins.

---

## Functional pipeline

This is the IWGC functional Snakefile, unchanged except for the MMseqs2 fix below. It runs:

1. **AGAT** to keep the longest isoform, then **gffread** `-J` to keep complete models only.
2. **MMseqs2** against NCBI (GenomesDatabase v1), a herbicide-resistance KO set and UniRef50, keeping the best hit per query at e < 1e-5.
3. **InterProScan 5.77-108.0** in 2,000-protein chunks (GO, pathways, IPR), **HMMER/Pfam**, **SignalP + TargetP**, and **MultiLoc2** (which uses the InterProScan GO terms).
4. `functional_merge.py` → `functional_annotation.tsv` with 11 columns.

`run_functional_v3.sb` adds the two post-processing steps used for the v2 release:

- reorders the TSV into **assembly order** (gffread output order) and keeps an alphabetical copy;
- writes a **functional GFF** with `agat_sp_manage_functional_annotation.pl`, using the NCBI best hits, a hit-only subset of the NCBI FASTA, and the InterProScan TSV without the pathways column.

### Bug fix: MMseqs2 race on `/tmp`

The `mmseqs2` rule runs once per database, and Snakemake runs the three searches **in parallel**. All three used `/tmp` as the MMseqs tmp directory. MMseqs creates a `latest` symlink inside its tmp directory, so the parallel runs raced, and one search died almost immediately with:

```
Could not create symlink of /tmp/11310392332275771123!
Error in rule mmseqs2: ... databases/current/GenomesDatabase.v1.prot
```

Because the pipeline runs with `--keep-going`, every other rule kept running for about a day before the missing NCBI search was noticed. The fix gives each database its own tmp directory inside the work volume:

```diff
-  mmseqs easy-search {input.aa} {input.database} {output.ofmt6} /tmp --threads 20 ...
+  rm -rf {volume}/tmp/mmseqs_{database_name} && mkdir -p {volume}/tmp/mmseqs_{database_name}
+  mmseqs easy-search {input.aa} {input.database} {output.ofmt6} /data/tmp/mmseqs_{database_name} --threads 20 ...
+  rm -rf {volume}/tmp/mmseqs_{database_name}
```

This was verified with `snakemake -n`. `fix_mmseqs_ncbi.sb` reruns only the NCBI search for a workdir that was started with the old Snakefile; with the fixed Snakefile it is not needed.

### Running it

```bash
cd /path/to/functional_workdir    # contains Snakefile, config.yaml, images/, databases/, inputs/
# config.yaml: volume_name, fasta_file, final_filtering_gff, species_name
sbatch run_functional_v3.sb
```

To set up a second workdir next to an existing one without using extra space or file quota, hard-link the databases and genome:

```bash
cp -al /old/functional/databases  /new/functional/databases
ln     /old/functional/inputs/genome.fasta /new/functional/inputs/
```

---

## Output: PhaMi v3 annotation

`results/PhaMi_v3/PhaMi.v3.gff3.gz` contains gene / mRNA / exon / CDS records:

| | |
|---|---|
| Genes | 44,965 |
| mRNAs | 49,216 |
| Longest-isoform proteins | 44,965 |
| Complete models used for functional annotation | 44,960 |
| Edge cases (out of 49,216) | 1 internal stop, 3 not starting with Met, 6 CDS not a multiple of 3 |

Genes per chromosome:

| Subgenome C | genes | Subgenome B | genes |
|---|---|---|---|
| Chr01C | 3,993 | Chr01B | 4,421 |
| Chr02C | 3,428 | Chr02B | 3,960 |
| Chr03C | 3,213 | Chr03B | 3,614 |
| Chr04C | 2,959 | Chr04B | 3,416 |
| Chr05C | 2,686 | Chr05B | 2,929 |
| Chr06C | 2,383 | Chr06B | 2,636 |
| Chr07C | 2,475 | Chr07B | 2,800 |
| **Total C** | **21,137** | **Total B** | **23,776** |

Chr00: 52 genes.

The full release lives on MSU HPCC at `/mnt/research/Weed_lab/projects/Phalaris_minor_v3/`. Besides the GFF3, it includes proteins, CDS, inputs, the full BRAKER working directory and the trained AUGUSTUS species `PhaMi_v3_A_isoseq`. The species can be reused by copying it into `$AUGUSTUS_CONFIG_PATH/species/`.

---

## Runtime: BRAKER3 vs MAKER

These figures come from SLURM accounting (`sacct`) for jobs that completed. Repeat masking is excluded from both, because v3 reused an existing soft-masked genome.

| | Wall-clock | CPU-hours (approx.) | Jobs |
|---|---|---|---|
| **BRAKER3 v3** | **~11 h** | **~530** | 1 BRAKER job (48 CPU) + 10 min alignment |
| MAKER v2 | ~5 days | ~12,000–13,000 | ~30 per-chromosome jobs (B subgenome ~2 days; C chromosomes 2.4–3.3 days each) |

MAKER also needed several failed or restarted attempts, including a 7-day timeout and a 2.5-day failed canary run. BRAKER3 ran in a single job without restarts.

Two caveats:

- The MAKER CPU-hour figure is a rough sum; a few Snakemake jobs are identified only by UUID.
- MAKER v2 additionally predicts UTRs and tRNAs, which BRAKER does not.

---

## HPCC notes: quotas and containers

- **Scratch file quota (1,048,576 files).** GeneMark splits a heavily masked 4.9 Gb genome into about 178 k pieces and writes two files per piece for every GeneMark run. That is about 357 k files per run, so two parallel BRAKER runs pushed usage from 138 k to **856 k files within a few hours**. `99_quota_guard.sbatch` removes `output/data` and `output/gmhmm` from `GeneMark-ES`/`GeneMark-EP` once their `genemark.gtf` has been stable for 15 min. Contents of `GeneMark-ETP` are only counted, never deleted. As a last resort, if usage passes 980 k files, the guard cancels the protein-only control run.
- **Full `$HOME`.** All temporary and container cache directories (`TMPDIR`, `APPTAINER_CACHEDIR`, `SINGULARITY_CACHEDIR`) point to scratch.
- **Symlinks and containers.** Symlinks into `/mnt/research/...` (an automounted, real path `/mnt/ufs18/...`) did not resolve inside Singularity, even with bind mounts. Large inputs that containers need are therefore copied as real files.
- **BUSCO** keeps about 10 k small files per run, so its run folders are packed into `.tar.gz` after the summary is copied out.
- **GeneMark license.** The GeneMark inside the braker3 container only *tests* for `~/.gm_key` and runs without it.

---

## Known issues and open questions

- **Fewer genes on the C subgenome than on B**, even though the C chromosomes are larger (21.1 k vs 23.8 k). The C subgenome is also more heavily masked (~92 % vs ~88 %). This may be over-masking by the uncurated EDTA library, which was not cleaned of gene fragments. Worth checking against Iso-Seq loci falling in masked regions.
- **45 k genes may be conservative for a tetraploid.** MAKER v2 has about 10 k more protein-coding genes at the same BUSCO completeness. The functional results (hit rates per gene set) should show whether those extra genes are real.
- **No UTRs or tRNAs** in v3.
- The NCBI protein FASTA `GenomesDatabase.v1.fixed.prot` is missing `>` on its very first header, so one record is ignored when building the AGAT hit subset.
- The BRAKER README now marks BRAKER as deprecated in favour of **BRAKER4**. The Iso-Seq workflow here uses the documented `braker3:isoseq` container (BRAKER 3.0.8).

---

## Requirements

- SLURM cluster (developed on MSU ICER HPCC).
- Singularity/Apptainer.
- Containers: `teambraker/braker3:latest`, `teambraker/braker3:isoseq`, gffread, plus the IWGC functional images (AGAT, MMseqs2, InterProScan 5.77-108.0, HMMER3, SignalP/TargetP, MultiLoc2, samtools).
- Modules: `minimap2/2.26`, `SAMtools/1.18`, `BUSCO/5.8.2`.
- Snakemake (conda env `iwgc`) for the functional pipeline.

Related repositories: [IWGC_Structural](https://github.com/danifilho/IWGC_Structural) (MAKER-based v2) and [IWGC_Funcional](https://github.com/danifilho/IWGC_Funcional).

---

## Contact

Maintained by **@danifilho**. Questions and issues are welcome through the GitHub issue tracker.
