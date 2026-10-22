# IWGC Genome Annotation Pipeline v3

Structural and functional genome annotation for the **International Weed Genomics Consortium (IWGC)**, built with Snakemake and Singularity on SLURM clusters (developed on MSU ICER HPCC).

- **Structural annotation:** repeat soft-masking, then **BRAKER3** with protein evidence plus **PacBio Iso-Seq** long reads, then a clean GFF3 with IWGC gene IDs, protein/CDS FASTA and BUSCO.
- **Functional annotation:** MMseqs2 (NCBI, UniRef50, herbicide-resistance set), InterProScan, Pfam, SignalP/TargetP and MultiLoc2. The results are merged into one TSV and a functional GFF.

v3 replaces MAKER (v2, [IWGC_Structural](https://github.com/danifilho/IWGC_Structural)) with BRAKER3. On *Phalaris minor* (4.9 Gb, allotetraploid) it reaches the same BUSCO completeness as MAKER (**95.9 % vs 96.1 %**) in **~11 h instead of ~5 days**. The functional module is [IWGC_Funcional](https://github.com/danifilho/IWGC_Funcional) plus bug fixes and two new output steps.

> **Validation (2026-10-22)**
>
> - **Structural `Snakefile`:** tested end-to-end on a single chromosome (Chr07B, 216 Mb). It ran with no errors in 3 h (48 CPUs), with BRAKER3 in ETP mode and AUGUSTUS accuracy 0.71. The individual steps also ran on the full 4.9 Gb genome (results below).
> - **Functional `Snakefile`:** completed on the full *P. minor* v3 proteome (44,960 proteins), including the new `order_tsv` and `agat_functional_gff` rules.

---

## Contents

1. [Pipeline overview](#1-pipeline-overview)
2. [Repository layout](#2-repository-layout)
3. [Requirements](#3-requirements)
4. [Container images](#4-container-images)
5. [Databases and evidence](#5-databases-and-evidence)
6. [Setup](#6-setup)
7. [Configuration](#7-configuration)
8. [Running](#8-running)
9. [Outputs](#9-outputs)
10. [Example results: *Phalaris minor*](#10-example-results-phalaris-minor)
11. [What changed from v2](#11-what-changed-from-v2)
12. [Troubleshooting](#12-troubleshooting)
13. [Citations](#13-citations)
14. [Contact](#14-contact)

---

## 1. Pipeline overview

```
STRUCTURAL  (structural/Snakefile)
  genome.fasta ─► RepeatModeler ─► RepeatMasker ─► bedtools soft-mask ─┐  (or a pre-masked genome)
  proteins ─────► protein DB ──────────────────────────────────────────┤
  Iso-Seq reads ─► minimap2 splice:hq ─► sorted BAM ───────────────────┴─► BRAKER3 (ETP, one genome-wide model)
                                                                              │
                         BUSCO ◄── proteins/CDS (gffread) ◄── GFF3 with IWGC IDs
                                                                              │
FUNCTIONAL  (functional/Snakefile)                                            ▼
  GFF3 ─► AGAT longest isoform ─► gffread -J proteins ─┬─► MMseqs2 × {NCBI, herbicide-resistance, UniRef50}
                                                       ├─► InterProScan (chunks) ─► MultiLoc2
                                                       ├─► HMMER / Pfam
                                                       └─► SignalP + TargetP
                         ─► functional_annotation.tsv (assembly order) + <SP>.functional.gff (AGAT)
```

Step-by-step details, with the result of every step on *P. minor*:
- [docs/STRUCTURAL.md](docs/STRUCTURAL.md)
- [docs/FUNCTIONAL.md](docs/FUNCTIONAL.md)

---

## 2. Repository layout

```
IWGC_annotation_v3/
├── structural/
│   ├── Snakefile                 # structural workflow
│   ├── config.yaml               # template: edit before running
│   ├── run_structural.sb         # SLURM wrapper (one node)
│   └── scripts/
│       ├── finalize_gff3.py      # braker.gtf -> GFF3 with IWGC IDs
│       ├── postprocess_proteins.sh
│       ├── gene_set_stats.sh
│       └── genemark_cleanup.sh   # keeps the file count under control during BRAKER
├── functional/
│   ├── Snakefile                 # functional workflow
│   ├── config.yaml               # template: edit before running
│   ├── run_functional.sb         # SLURM wrapper (one node)
│   └── scripts/
│       ├── functional_merge.py   # merges all results into the TSV
│       ├── ml2_one_by_one.sh     # MultiLoc2 driver
│       └── extract_hits.py       # NCBI hit subset for AGAT (repairs glued FASTA headers)
├── docs/
│   ├── STRUCTURAL.md             # every structural step + P. minor results
│   ├── FUNCTIONAL.md             # every functional step + P. minor results
│   ├── TROUBLESHOOTING.md        # errors seen in practice and how to fix them
│   └── CASE_STUDY_PhaMi.md       # why an earlier BRAKER run failed; BRAKER3 vs MAKER
└── results/PhaMi/
    ├── structural/               # PhaMi.v3.gff3.gz, stats, BUSCO, Iso-Seq flagstat
    ├── functional/               # functional TSV + GFF
    └── comparison/               # BUSCO/stats of MAKER v2 and of the BRAKER2 runs
```

---

## 3. Requirements

| Software | Version used | Notes |
|---|---|---|
| Snakemake | 9.13.7 | conda env (`iwgc`); also needs `biopython` 1.86 and `pandas` 2.3.3 |
| SingularityCE / Apptainer | SingularityCE 4.1.2 (Apptainer 1.4.1 also used for pulls) | |
| SLURM | — | optional; the `run_*.sb` wrappers run Snakemake on a single node |
| Python | 3.13 (host) | only the standard library is used by the host-side scripts |

```bash
conda create -n iwgc -c conda-forge -c bioconda snakemake=9.13.7 biopython pandas
```

**Hardware used for a 4.9 Gb genome:**
- Structural: 1 node, 64 CPUs, 400 GB RAM. BRAKER itself uses 48 threads.
- Functional: 1 node, 64 CPUs, 250 GB RAM.

Plan for a scratch **file** quota as well as space; see [TROUBLESHOOTING #2](docs/TROUBLESHOOTING.md#2-scratch-file-quota-explodes-during-braker-genemark).

---

## 4. Container images

All images go in `images/` inside each working directory, **with the file names below**. Library images are unsigned, so pull them with `-U`.

### Structural

| File in `images/` | Pull command | Tools inside (verified versions) |
|---|---|---|
| `braker3_isoseq.sif` | `singularity pull braker3_isoseq.sif docker://teambraker/braker3:isoseq` | BRAKER 3.0.8, GeneMark-ETP (long-read StringTie fork), AUGUSTUS 3.5.0, StringTie 2.2.1, DIAMOND 2.0.15, TSEBRA 1.1.2.4, samtools 1.13 |
| `braker3.sif` (only for protein-only runs) | `singularity pull braker3.sif docker://teambraker/braker3:latest` | BRAKER 3.0.8, AUGUSTUS 3.5.0, GeneMark-ETP/EP, TSEBRA |
| `minimap2.sif` | `singularity pull minimap2.sif docker://staphb/minimap2:2.26` | minimap2 2.26-r1175 |
| `samtools.sif` | `singularity pull samtools.sif docker://staphb/samtools:1.19` | samtools 1.19 |
| `gffread.sif` | `singularity pull -U gffread.sif library://danifilho/functional_annotation_images/gffread:latest` | gffread 0.12.7 |
| `busco.sif` | `singularity pull busco.sif docker://ezlabgva/busco:v5.8.2_cv1` | BUSCO 5.8.2 |
| `dfam_tetools.sif` (masking) | `singularity pull dfam_tetools.sif docker://dfam/tetools` | RepeatModeler 2.0.7, RepeatMasker 4.2.3 |
| `bedtools.sif` (masking) | `singularity pull -U bedtools.sif library://danifilho/structural_annotation/bedtools:latest` | bedtools 2.30.0 |

### Functional

Source: [library://danifilho/functional_annotation_images](https://cloud.sylabs.io/library/danifilho/functional_annotation_images).

| File in `images/` | Pull command | Tools inside (verified versions) |
|---|---|---|
| `agat.sif` | `singularity pull -U agat.sif library://danifilho/functional_annotation_images/agat:latest` | AGAT 1.4.0 |
| `gffread.sif` | `singularity pull -U gffread.sif library://danifilho/functional_annotation_images/gffread:latest` | gffread 0.12.7 |
| `mmseqs2.sif` | `singularity pull -U mmseqs2.sif library://danifilho/functional_annotation_images/mmseqs2:latest` | MMseqs2 (commit `ad5837b3444728411e6c90f8c6ba9370f665c443`) |
| `iprscan.sif` | `singularity pull -U iprscan.sif library://danifilho/functional_annotation_images/iprscan:latest` | InterProScan 5.77-108.0 |
| `hmmer3.sif` | `singularity pull -U hmmer3.sif library://danifilho/functional_annotation_images/hmmer3:latest` | HMMER 3.3 |
| `sigtarp.sif` | `singularity pull -U sigtarp.sif library://danifilho/functional_annotation_images/sigtarp:latest` | SignalP 5.0b, TargetP 2.0 |
| `multiloc2_v3.img` | `singularity pull -U multiloc2_v3.img library://danifilho/functional_annotation_images/multiloc2:latest` | MultiLoc2 (Python 2.7.18) |
| `samtools.sif` | `singularity pull -U samtools.sif library://danifilho/functional_annotation_images/samtools:latest` | samtools 1.19.2 |

> Versions were read from inside each image (`--version`) on 2026-10-21. `:latest` and `dfam/tetools` tags can move; check `--version` after pulling a new copy.

---

## 5. Databases and evidence

### Structural evidence (`structural/inputs/`)

| Input | Example (*P. minor*) | Notes |
|---|---|---|
| Genome FASTA | `phalaris_swapped.fasta`, 4.92 Gb, 15 sequences | Simple headers only (no spaces). |
| Close-relative proteins | *Lolium rigidum* (APGP/CSIRO), 54,112 | |
| Broad protein set | [OrthoDB Viridiplantae](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/), 5.3 M | **Required in practice**: BRAKER's ProtHint needs it. |
| Iso-Seq reads | `*.flnc.bam` from `isoseq refine`, 4.1 M reads | Also accepted: FASTA/FASTQ(.gz) of FLNC or HiFi reads. **Not** collapsed transcripts. |
| BUSCO lineage | `busco_downloads/lineages/poales_odb10` | Download beforehand; compute nodes run BUSCO `--offline`. |

Download OrthoDB Viridiplantae and a BUSCO lineage:
```bash
wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Viridiplantae.fa.gz && gunzip Viridiplantae.fa.gz
singularity exec images/busco.sif busco --download poales_odb10    # creates busco_downloads/
```

### Functional databases (`functional/databases/`)

| Path | Content | Used by |
|---|---|---|
| `current/GenomesDatabase.v1.prot` | IWGC NCBI protein set (5.99 M proteins, FASTA + MMseqs index) | MMseqs2 `ncbi`: descriptions and AGAT product names |
| `herb_res_prots/total_KO_labeled.fasta` | herbicide-resistance KO protein set | MMseqs2 `herbres` |
| `mmseqs_db/uniref50` | UniRef50 MMseqs2 database (`mmseqs databases UniRef50 ...`) | MMseqs2 `uniref` |
| `Pfam-A.hmm` (+ `hmmpress` files) | Pfam-A HMMs | `hmmscan` |
| `interproscan-5.77-108.0/data` | InterProScan data release **matching the image** (5.77-108.0) | InterProScan |

Setting up a second working directory? Hard-link the databases instead of copying them. On the same filesystem this costs neither space nor file quota:
```bash
cp -al /path/old_workdir/databases /path/new_workdir/databases
```

---

## 6. Setup

Each module runs in its **own working directory** (`volume_name`). Everything the containers read must be a **real file or hard link inside it**; see [TROUBLESHOOTING #1](docs/TROUBLESHOOTING.md#1-genome-file--does-not-exist-inside-a-container).

```
structural_workdir/                     functional_workdir/
├── Snakefile                           ├── Snakefile
├── config.yaml                         ├── config.yaml
├── run_structural.sb                   ├── run_functional.sb
├── scripts/                            ├── scripts/
├── images/      (structural images)    ├── images/      (functional images)
├── inputs/      genome, proteins,      ├── inputs/      genome + <SP>.v3.gff3
│                Iso-Seq reads          └── databases/   (see above)
└── busco_downloads/lineages/<lineage>/
```

```bash
git clone https://github.com/danifilho/IWGC_annotation_v3.git
mkdir -p $SCRATCH/PhaMi_v3 && cd $SCRATCH/PhaMi_v3
cp -r IWGC_annotation_v3/structural structural
cp -r IWGC_annotation_v3/functional functional
# then: pull images into */images, place inputs, edit both config.yaml files
```

---

## 7. Configuration

### `structural/config.yaml`

| Key | Meaning | *P. minor* value |
|---|---|---|
| `volume_name` | absolute path of the structural workdir | `/mnt/gs21/scratch/.../structural` |
| `species_abbreviation` | IWGC code, used as the gene-ID prefix | `PhaMi` |
| `fasta_file` | unmasked genome in `inputs/` | `phalaris_swapped.fasta` |
| `softmasked_genome` | optional pre-masked genome in `inputs/`; `""` = run masking | `""` |
| `protein_files` | list of protein FASTA files in `inputs/` | Lolium + Viridiplantae |
| `isoseq_reads` | FLNC/HiFi reads in `inputs/` (`.bam`, `.fa[.gz]`, `.fq[.gz]`); `""` = proteins only | `PhaMi.flnc.bam` |
| `repeatmodeler_threads`, `repeatmasker_threads` | masking threads | 64 |
| `minimap2_threads` / `minimap2_index_batch` | alignment threads / `-I` (must exceed the genome size) | 32 / `16G` |
| `braker_threads` | BRAKER threads (**≤ 48**) | 48 |
| `braker_species` | name of the new AUGUSTUS model | `PhaMi_v3` |
| `drop_internal_stop_models` | drop models whose BRAKER protein has an internal stop codon (before ID numbering) | `true` |
| `busco_lineage` | local lineage directory | `busco_downloads/lineages/poales_odb10` |
| `busco_max_protein_length` | longer proteins are excluded from BUSCO | 10000 |
| `busco_threads` | BUSCO threads | 16 |

### `functional/config.yaml`

| Key | Meaning |
|---|---|
| `volume_name` | absolute path of the functional workdir |
| `fasta_file` | genome in `inputs/` |
| `final_filtering_gff` | structural GFF3 in `inputs/` (`<SP>.v3.gff3`) |
| `species_name` | output prefix |
| `mmseqs_databases` | list of `{name, path}`; the `ncbi` entry is required and also needs `ncbidata` (FASTA with descriptions) |
| `multiloc_script` | `ml2_one_by_one.sh` (do not change) |

---

## 8. Running

### Structural

```bash
cd structural_workdir
snakemake -n                       # dry run: check the job list
sbatch run_structural.sb           # or run directly:
# snakemake --cores 64 --use-singularity --singularity-args "-B $PWD" --rerun-triggers mtime --keep-going --rerun-incomplete -p
```

Check early that BRAKER runs in the expected mode:
```bash
grep "will be executed in" braker_outputs/braker.log
# Both protein and RNA-Seq data in input detected. BRAKER will be executed in ETP mode (BRAKER3).
```

### Functional

```bash
cp structural_workdir/final_outputs/PhaMi.v3.gff3 functional_workdir/inputs/
ln structural_workdir/inputs/<genome>.fasta functional_workdir/inputs/     # hard link
cd functional_workdir
snakemake -n
sbatch run_functional.sb
```

After each run, check for failed rules (`--keep-going` does not stop on errors):
```bash
grep -n "Error in rule" *.err
```

---

## 9. Outputs

### Structural (`structural_workdir/`)

| Path | Content |
|---|---|
| `final_outputs/<SP>.v3.gff3` | **final annotation** (gene / mRNA / exon / CDS; IWGC IDs) |
| `final_outputs/<SP>.v3.proteins.fa`, `<SP>.v3.cds.fa` | all isoforms, same IDs as the GFF3 |
| `final_outputs/<SP>.v3.longest.proteins.fa` | one protein per gene |
| `final_outputs/<SP>.v3.stats.txt` | gene-set statistics and sanity checks |
| `final_outputs/<SP>.v3.excluded_internal_stop.txt` | BRAKER transcripts dropped for internal stop codons |
| `busco_outputs/<SP>.v3.busco_summary.txt` | BUSCO summary (+ `busco_<SP>.tar.gz`) |
| `braker_outputs/` | full BRAKER run (`braker.gtf`, `braker.log`, hints, GeneMark-ETP, AUGUSTUS) |
| `braker_inputs/augustus_config/species/<braker_species>` | **trained AUGUSTUS model** (reusable) |
| `braker_inputs/isoseq.minimap2.bam`, `isoseq.flagstat.txt` | Iso-Seq alignment |
| `masking_outputs/`, `repeat_modeler_outputs/`, `repeat_masker_outputs/` | soft-masked genome, repeat library, repeat GFF |

Gene IDs follow the IWGC scheme, with the counter increasing by 10 across the genome in FASTA order:

```
PhaMiChr01Cg000010                  gene
PhaMiChr01Cg000010-mRNA-1           mRNA   (braker_id=g4474.t1)
PhaMiChr01Cg000010-mRNA-1:1         exon
PhaMiChr01Cg000010-mRNA-1:cds       CDS
```

### Functional (`functional_workdir/functional_outputs/`)

| Path | Content |
|---|---|
| `functional_annotation.tsv` | **one row per protein, in assembly order** (see the columns below) |
| `functional_annotation.alphabetical.tsv` | same content, alphabetical order |
| `<SP>.functional.gff` | longest-isoform GFF with NCBI product names and InterPro/GO/Pfam terms (AGAT) |

TSV columns:
- `ID sequence`
- `NCBI DESCRIPTION`, `NCBI Subject;E-value;Bit`
- `herbres Subject;E-value;Bit`, `uniref Subject;E-value;Bit`
- `SignalP Pos;Pr`, `TargetP ...`, `MULTILOC`
- `IPRSCAN GO`, `IPRSCAN IPR`, `Hmmer Pfam`

Intermediate outputs of every step are listed in [docs/FUNCTIONAL.md](docs/FUNCTIONAL.md).

---

## 10. Example results: *Phalaris minor*

Files in [`results/PhaMi/`](results/PhaMi/):
- [`PhaMi.v3.gff3.gz`](results/PhaMi/structural/PhaMi.v3.gff3.gz): the structural GFF3.
- [`PhaMi.v3.functional_annotation.tsv.gz`](results/PhaMi/functional/PhaMi.v3.functional_annotation.tsv.gz): the functional table.
- [`PhaMi.v3.functional.gff.gz`](results/PhaMi/functional/PhaMi.v3.functional.gff.gz): the functional GFF.

| | **v3 (this pipeline)** | MAKER v2 |
|---|---|---|
| Genes / mRNAs | 44,965 / 49,216 | 55,124 protein-coding (+ 2,039 tRNA) |
| BUSCO `poales_odb10` C [S, D] | **95.9 % [18.6, 77.3]** | 96.1 % [17.6, 78.4] |
| BUSCO F / M | 0.2 % / 3.9 % | 1.5 % / 2.5 % |
| Single-exon transcripts | 25.3 % | – |
| Median protein (longest isoform) | 353 aa | – |
| AUGUSTUS accuracy | 0.75 | – |
| Gene-prediction wall-clock | **~11 h** (1 job, 48 CPUs) | ~5 days (~30 jobs, ~12–13 k CPU-h) |
| Functional: NCBI / UniRef50 hit | 96.8 % / 96.9 % | 97.3 % / 97.7 % |
| Functional: InterPro / GO | 87.5 % / 70.2 % | 89.1 % / 69.6 % |
| Functional: herbicide-resistance hits | 429 | 291 |

Before v3, BRAKER had been run on this genome in protein-only mode with 15 per-chromosome models. That run reached only C:67.7 % / M:30.1 %. [docs/CASE_STUDY_PhaMi.md](docs/CASE_STUDY_PhaMi.md) explains why, and which v3 design choices follow from it.

---

## 11. What changed from v2

**Structural (vs. MAKER v2 and vs. a naive BRAKER setup)**
- MAKER is replaced by **BRAKER3 in ETP mode**, with Iso-Seq **reads** aligned by minimap2 `splice:hq`, as BRAKER requires for long reads.
- **One genome-wide training run** instead of per-chromosome runs (AUGUSTUS accuracy 0.75 vs 0.21–0.42).
- Protein DB = close relative **+ OrthoDB Viridiplantae**, with every file listed explicitly in the config.
- A new IWGC-ID GFF3 writer (`finalize_gff3.py`) that handles BRAKER's mixed AUGUSTUS/GeneMark GTF and computes CDS phases.
- A built-in GeneMark file-count cleanup, BAI→CSI fallback for large chromosomes, and a BUSCO guard against giant proteins.
- Gene models whose BRAKER protein has an internal stop codon (a few GeneMark.hmm3 predictions) are dropped before ID numbering (`drop_internal_stop_models`).

**Functional (vs. IWGC_Funcional)**
- **Bug fix:** the parallel MMseqs2 searches shared `/tmp` and raced on MMseqs' `latest` symlink, which **silently dropped the NCBI search**. Each database now has its own tmp directory.
- **New rule `order_tsv`:** the TSV is written in assembly order, taken from the input GFF3 (AGAT sorts sequence names alphabetically, so its order cannot be used). An alphabetical copy is kept.
- **New rule `agat_functional_gff`:** functional GFF via AGAT, with `extract_hits.py` repairing the glued headers of the NCBI FASTA on the fly.
- The duplicate `configfile:` line was removed, and `egrep \>` became `grep '^>'` (it raised a Python `SyntaxWarning`).

---

## 12. Troubleshooting

[docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) lists every error met while running this pipeline: the message, the cause and the fix. It covers 25 entries, including container paths, file quotas, BRAKER modes, BUSCO 0 %, the MMseqs `/tmp` race, AGAT and InterProScan.

---

## 13. Citations

If you use this pipeline, please cite the underlying tools:
- **BRAKER3:** Gabriel L. et al. (2024) *Genome Research* 34:769–777.
- **GeneMark-ETP:** Brůna T. et al. (2024) *Genome Research* 34:757–768.
- **AUGUSTUS:** Stanke M. et al. (2008) *Bioinformatics* 24:637–644.
- **TSEBRA:** Gabriel L. et al. (2021) *BMC Bioinformatics* 22:566.
- **minimap2:** Li H. (2018) *Bioinformatics* 34:3094–3100.
- **BUSCO:** Manni M. et al. (2021) *Mol Biol Evol* 38:4647–4654.
- **RepeatModeler2:** Flynn J.M. et al. (2020) *PNAS* 117:9451–9457.
- **InterProScan:** Jones P. et al. (2014) *Bioinformatics* 30:1236–1240.
- **MMseqs2:** Steinegger M. & Söding J. (2017) *Nat Biotechnol* 35:1026–1028.
- **AGAT:** Dainat J. *AGAT* (Zenodo).
- **SignalP 5.0:** Almagro Armenteros J.J. et al. (2019) *Nat Biotechnol* 37:420–423.
- **TargetP 2.0:** Almagro Armenteros J.J. et al. (2019) *Life Sci Alliance* 2:e201900429.
- **MultiLoc2:** Blum T. et al. (2009) *BMC Bioinformatics* 10:274.
- **HMMER:** Eddy S.R. (2011) *PLoS Comput Biol* 7:e1002195.
- **Pfam:** Mistry J. et al. (2021) *Nucleic Acids Res* 49:D412–D419.

## 14. Contact

- **Email:** dasilvaf@msu.edu
- **GitHub:** [danifilho](https://github.com/danifilho)
