# Structural pipeline: step by step

Snakefile: [`structural/Snakefile`](../structural/Snakefile) · config: [`structural/config.yaml`](../structural/config.yaml)

```
genome.fasta ──► repeatmodeler ──► repeatmasker ──► softmask ─┐   (skipped if softmasked_genome is set)
                                                              ▼
protein_files ──► protein_db ───────────────────────────► braker ──► finalize_gff3 ──► extract_sequences ──► busco
isoseq_reads ──► isoseq_fasta ──► minimap2 ─────────────────┘  ▲
                                         augustus_config ──────┘
```

Each section below gives what the step does, the tool, its outputs, and **the result on *Phalaris minor*** (allotetraploid, 4.92 Gb, 15 sequences: Chr01C–Chr07C, Chr01B–Chr07B, Chr00).

---

## 1. Repeat masking: `repeatmodeler` → `repeatmasker` → `softmask`

These are the same rules as the IWGC v2 structural pipeline. They are skipped when `softmasked_genome` is set in the config.

| Rule | Tool (image) | Output |
|---|---|---|
| `repeatmodeler` | `BuildDatabase` + `RepeatModeler -LTRStruct` (dfam_tetools: RepeatModeler 2.0.7) | `repeat_modeler_outputs/<genome>-families.fa` (de novo repeat library) |
| `repeatmasker` | `RepeatMasker -gff -a -lib <library>` (dfam_tetools: RepeatMasker 4.2.3) | `repeat_masker_outputs/<genome>.out.gff` |
| `softmask` | `bedtools maskfasta -soft` (bedtools 2.30.0) | `masking_outputs/<genome>.softmask.fasta` (prints the % masked) |
| `faidx_masked` | `samtools faidx` | `.fai` of the masked genome (used for gene ordering) |

**Soft-masking, never hard-masking.** BRAKER uses lowercase as "repeat" information, but still allows genes there.

***P. minor*:**
- **MAKER v2:** the RepeatModeler library had 1,014 classified families.
- **v3 run:** reused an EDTA library + RepeatMasker soft-masked genome (`softmasked_genome`). Masking was 87–88 % per B-subgenome chromosome and 91–92 % per C-subgenome chromosome, **90.6 % overall** (4.46 of 4.92 Gb).

## 2. Protein evidence: `protein_db`

Concatenates every file in `protein_files` and keeps only the first word of each header.

***P. minor*:** *Lolium rigidum* (54,112) + OrthoDB Viridiplantae (5,310,477) = **5,364,589 proteins** (2.4 GB).

> A single related species is not enough for BRAKER's ProtHint. Include a broad set such as OrthoDB.

## 3. Iso-Seq evidence: `isoseq_fasta` → `minimap2`

| Rule | Tool | Output |
|---|---|---|
| `isoseq_fasta` | `samtools fasta` (BAM input) or passthrough (FASTA/FASTQ) | temporary `braker_inputs/isoseq_reads.fa.gz` |
| `minimap2` | `minimap2 -ax splice:hq -uf -I 16G` (2.26) → `samtools sort` (1.19) | `braker_inputs/isoseq.minimap2.bam` (+ `.bai`, or `.csi` for very large chromosomes), `isoseq.flagstat.txt` |

These are exactly the minimap2 options the BRAKER documentation requires for long reads. Use **reads** (FLNC or HiFi), not collapsed transcripts.

***P. minor*:** 4,116,976 FLNC reads → **99.84 % primary mapped**. There were also 2.6 M secondary alignments, which is expected: a read also aligns to the homeolog in the other subgenome. The step took 10 min on 32 threads with a 31 GB peak RAM.

## 4. Gene prediction: `augustus_config` → `braker`

**One genome-wide BRAKER run**, which trains a single GeneMark model and a single AUGUSTUS model for the whole genome:

```bash
braker.pl --genome=<softmasked> --prot_seq=braker_inputs/proteins.faa \
          --bam=braker_inputs/isoseq.minimap2.bam \
          --species=<braker_species> --workingdir=braker_outputs \
          --AUGUSTUS_CONFIG_PATH=braker_inputs/augustus_config --threads=48 --gff3
```

- With `isoseq_reads` set: container `teambraker/braker3:isoseq`, **ETP mode (BRAKER3)**. The GeneMark-ETP fork assembles the long reads with StringTie.
- With `isoseq_reads: ""`: container `teambraker/braker3:latest`, **EP mode (BRAKER2)**, proteins only.
- `augustus_config` copies the container's AUGUSTUS config to a writable location, which BRAKER needs for training.
- `scripts/genemark_cleanup.sh` runs in the background and keeps the file count under control (see [TROUBLESHOOTING #2](TROUBLESHOOTING.md#2-scratch-file-quota-explodes-during-braker-genemark)).

Main outputs in `braker_outputs/`: `braker.gtf`, `braker.aa`, `braker.codingseq`, `braker.gff3`, `braker.log`, `hintsfile.gff`, `GeneMark-ETP/`, `Augustus/`. The trained model is in `braker_inputs/augustus_config/species/<braker_species>`.

***P. minor* (whole genome, 48 threads, 400 GB):**

| | BRAKER3 + Iso-Seq | BRAKER2 (proteins only) |
|---|---|---|
| Mode in `braker.log` | ETP (BRAKER3) | EP (BRAKER2) |
| AUGUSTUS accuracy, initial → optimized | 0.732 → **0.754** | 0.748 → 0.764 |
| Wall-clock | **10 h 50 min** | 15 h 56 min |
| Genes / transcripts | **44,965 / 49,216** | 96,759 / 102,236 |

The earlier per-chromosome BRAKER2 run, with *Lolium* proteins only, had AUGUSTUS accuracy of just 0.21–0.42.

## 5. Lab-format GFF3: `finalize_gff3`

[`scripts/finalize_gff3.py`](../structural/scripts/finalize_gff3.py) turns `braker.gtf` into a clean GFF3 that follows the **IWGC/lab ID scheme**. `braker.gtf` mixes AUGUSTUS, GeneMark.hmm3 and GeneMarkS-T records.

```
Chr01C  BRAKER3  gene  1429105 1430284 . - . ID=PhaMiChr01Cg000010;Name=PhaMiChr01Cg000010
Chr01C  BRAKER3  mRNA  1429105 1430284 . - . ID=PhaMiChr01Cg000010-mRNA-1;Name=PhaMiChr01Cg000010-mRNA-1;Parent=PhaMiChr01Cg000010;braker_id=g4474.t1
Chr01C  BRAKER3  exon  1429893 1430284 . - . ID=PhaMiChr01Cg000010-mRNA-1:1;Parent=PhaMiChr01Cg000010-mRNA-1
Chr01C  BRAKER3  CDS   1429893 1430284 . - 0 ID=PhaMiChr01Cg000010-mRNA-1:cds;Parent=PhaMiChr01Cg000010-mRNA-1
```

- **Gene order:** genome FASTA order, then start coordinate.
- **Gene IDs:** `<species><seqid>g<NNNNNN>`, with the counter increasing **by 10 across the whole genome** and never reset per chromosome.
- **Transcripts:** isoforms are `-mRNA-1..n`, ordered by start coordinate; `braker_id=` keeps the original BRAKER ID.
- **Exons and CDS:** BRAKER predicts no UTRs, so exons equal CDS. Stop codons are folded into the CDS whenever a source omits them. CDS phase is computed from the cumulative coding length.

***P. minor*:** 44,965 genes, 49,216 mRNAs, 0 warnings. Phases were verified against gffread.

## 6. Sequences and statistics: `extract_sequences`

`gffread -S -x -y` (gffread 0.12.7) extracts sequences. Then `postprocess_proteins.sh` strips the terminal `*` and builds the longest-isoform set, and `gene_set_stats.sh` writes the statistics.

| Output | Content |
|---|---|
| `final_outputs/<SP>.v3.gff3` | final annotation |
| `<SP>.v3.proteins.fa` / `<SP>.v3.cds.fa` | all isoforms, same IDs and order as the GFF3 |
| `<SP>.v3.longest.proteins.fa` | one protein per gene |
| `<SP>.v3.stats.txt` | counts, single-exon %, exons per transcript, median protein length, sanity checks, genes per sequence, BRAKER mode and accuracy |

***P. minor*:**

| Metric | v3 | What it tells you |
|---|---|---|
| Single-exon transcripts | **25.3 %** | 20–30 % is typical for grasses; above 50 % suggests fragments or TE ORFs |
| CDS exons per transcript | **5.29** | |
| Median protein length (longest isoform / all isoforms) | **353 / 358 aa** | Below about 250 aa suggests fragmented models |
| Internal stop / not starting with M / CDS not ×3 | 1 / 3 / 6 (of 49,216) | Edge cases |

Genes per chromosome: C subgenome 21,137, B subgenome 23,776, Chr00 52. Full table in [`results/PhaMi/structural/`](../results/PhaMi/structural/).

## 7. Completeness: `busco`

`busco -m proteins --offline` (BUSCO 5.8.2) on the longest-isoform proteins, after dropping proteins longer than `busco_max_protein_length` (see [TROUBLESHOOTING #16](TROUBLESHOOTING.md#16-busco-reports-0-c0-f0-m100)). The run directory is archived as `busco_outputs/busco_<SP>.tar.gz`, because BUSCO leaves about 10 k small files.

***P. minor*, `poales_odb10` (n = 4,896):**

| Annotation | C [S, D] | F | M |
|---|---|---|---|
| **v3 (this pipeline)** | **95.9 % [18.6, 77.3]** | 0.2 % | 3.9 % |
| MAKER v2 (IWGC_Structural) | 96.1 % [17.6, 78.4] | 1.5 % | 2.5 % |
| BRAKER2, 1 model, full proteins (control) | 71.5 % [16.8, 54.7] | 5.4 % | 23.1 % |
| Earlier BRAKER2 run: 15 per-chromosome models, *Lolium* only | 67.7 % [23.0, 44.7] | 2.3 % | 30.1 % |

A high D is expected for an allopolyploid, because both homeologs are annotated.

---

## Resources and runtime (*P. minor*, 4.9 Gb)

| Step | CPUs | Peak RAM | Wall-clock |
|---|---|---|---|
| minimap2 + sort | 32 | 31 GB | 10 min |
| BRAKER3 (ETP) | 48 | < 400 GB requested | 10 h 50 min |
| finalize + gffread | 1 | < 4 GB | ~1 min |
| BUSCO | 16 | < 64 GB | ~10 min |

For comparison, MAKER v2 needed about 5 days of wall-clock and roughly 12–13 k CPU-hours across about 30 per-chromosome jobs (repeat masking excluded from both).
