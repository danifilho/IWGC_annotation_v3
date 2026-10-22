# Case study: *Phalaris minor*, and why the v3 design looks the way it does

*Phalaris minor* (`PhaMi`) is an allotetraploid weed with a 4.92 Gb assembly: 7 C-subgenome chromosomes, 7 B-subgenome chromosomes and an unplaced Chr00. It was the first genome annotated with this pipeline. Before v3, BRAKER had been tried on it and gave poor BUSCO scores. This page documents that diagnosis, because each problem became a design rule of v3.

## The earlier BRAKER run

| Observation (from its `braker.log` files, scripts and outputs) | Consequence |
|---|---|
| `Only Protein input detected, BRAKER will be executed in EP mode (BRAKER2)` in every log. | The Iso-Seq was **never used for training**; it was only given to TSEBRA afterwards, which can only select among existing models. |
| `braker_prot.faa` = 54,112 proteins (*Lolium rigidum* only); the OrthoDB set was never concatenated. | Weak protein evidence for ProtHint. |
| One BRAKER run **per chromosome** (15 independent models; Chr00 reused the Chr01C model). | AUGUSTUS accuracy 0.21–0.42; the C subgenome was worst (Chr06C 0.21). |
| 62 % single-exon transcripts, median protein 221 aa, 2.7 CDS exons per transcript, "genes" spanning ~500 kb, 216 proteins > 20 k aa. | Fragmented gene models plus TE-derived ORFs. |
| BUSCO 0 % on the merged `braker_all.aa`. | A `hmmsearch` crash on proteins > 100 k aa, not a result. The real post-TSEBRA BUSCO was C:67.7 % / M:30.1 %. |

The first explanation offered was that "the Iso-Seq data was not enough". The logs show something different: the Iso-Seq data **never reached BRAKER**.

## Controlled comparison

Same soft-masked genome (EDTA + RepeatMasker, 90.6 % masked). BUSCO `poales_odb10`.

| Run | Evidence | Models | AUGUSTUS acc. | Genes | Single-exon | Median aa | BUSCO C | BUSCO M |
|---|---|---|---|---|---|---|---|---|
| Earlier run | Lolium proteins (BRAKER2) + Iso-Seq in TSEBRA | 15 (per chr.) | 0.21–0.42 | ~55.6 k | 62 % | 221 | 67.7 % | 30.1 % |
| Control | Lolium + OrthoDB (BRAKER2) | 1 (genome) | 0.76 | 96.8 k | 38 % | 219 | 71.5 % | 23.1 % |
| **v3** | Lolium + OrthoDB + **Iso-Seq FLNC reads (BRAKER3)** | 1 (genome) | 0.75 | **45.0 k** | **25 %** | **353** | **95.9 %** | **3.9 %** |
| MAKER v2 | est2genome + protein2genome | per chromosome | – | 55.1 k | – | – | 96.1 % | 2.5 % |

How to read the table:
- **Earlier run → control** (one genome-wide model and the broad protein set): the models train well, with accuracy rising from about 0.4 to 0.76. However, BUSCO barely moves (67.7 % → 71.5 %). Proteins alone are not enough on this genome. Both factors changed at once, so their separate effects cannot be told apart.
- **Control → v3** (adding the Iso-Seq reads to training): BUSCO jumps from 71.5 % to 95.9 %, and every structural metric becomes plant-like.

The conclusion is that the Iso-Seq data was sufficient; it simply has to go into BRAKER3 itself, as reads.

## Runtime: BRAKER3 vs MAKER

Wall-clock and CPU time come from SLURM accounting (`sacct`). Repeat masking is excluded from both.

| | Wall-clock | CPU-hours (approx.) | Jobs |
|---|---|---|---|
| **BRAKER3 (v3)** | **~11 h** | **~530** | 1 BRAKER job (48 CPUs) + 10 min minimap2 |
| MAKER (v2) | ~5 days | ~12,000–13,000 | ~30 per-chromosome jobs: B subgenome ~2 days; C chromosomes 2.4–3.3 days each |

MAKER also needed failed and restarted attempts (a 7-day timeout and a 2.5-day failed canary). Two caveats: MAKER additionally predicts UTRs and tRNAs, and the MAKER CPU-hour figure is a rough sum.

## Open questions

- **Fewer genes on the C subgenome than on B** (21.1 k vs 23.8 k), although the C chromosomes are larger. The C subgenome is also more heavily masked (~92 % vs ~88 %), so possible over-masking by an uncurated EDTA library should be checked.
- **45 k genes vs MAKER's 55 k at the same BUSCO.** The functional pipeline answers this in part. After its filters, MAKER v2 keeps 47,494 complete proteins vs 44,960 for v3, and hit *rates* are almost identical (NCBI 97.3 vs 96.8 %, InterPro 89.1 vs 87.5 %). In absolute numbers, MAKER v2 has about 2,900 more proteins with InterPro domains (42,301 vs 39,362). This suggests that part of MAKER's extra models are real, domain-bearing genes that BRAKER3 missed. MAKER may also split some genes into several models, though, and about 7,600 of its models were partial and removed by `gffread -J`. Comparing the two gene sets locus by locus (e.g. with gffcompare) is the next step.
