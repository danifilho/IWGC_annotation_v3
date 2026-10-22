# Functional pipeline: step by step

Snakefile: [`functional/Snakefile`](../functional/Snakefile) · config: [`functional/config.yaml`](../functional/config.yaml)

```
<SP>.v3.gff3 ─► agat (longest isoform) ─► gff_read (-J: complete models only) ─► <SP>.aa
   <SP>.aa ─┬─► mmseqs2 × {ncbi, herbres, uniref} ─► best_hits (e < 1e-5) ──────────────┐
            ├─► interproscan (2,000-protein chunks) ─► prot_list ─► multiloc ───────────┤
            ├─► hmmer (Pfam) ───────────────────────────────────────────────────────────┤
            └─► sigtarp (SignalP 5 + TargetP 2) ────────────────────────────────────────┤
                                                                                         ▼
                         generate_tsv (functional_merge.py) ─► order_tsv ─► functional_annotation.tsv
   agat GFF + ncbi best hits + InterProScan ─► agat_functional_gff ─► <SP>.functional.gff
```

Each section gives what the step does, the tool, its outputs, and the result on ***Phalaris minor* v3** (input: 44,965 genes / 49,216 mRNAs). All numbers below come from the completed *P. minor* v3 run.

---

## 1. `agat`: longest isoform per gene

`agat_sp_keep_longest_isoform.pl` (AGAT 1.4.0) → `agat_outputs/<SP>.gff`

***P. minor*:** 44,965 genes → 44,965 mRNAs (one per gene).

## 2. `gff_read`: protein and CDS extraction

`gffread -J -S -y / -x` (gffread 0.12.7) → `gff_read_output/<SP>.aa`, `<SP>.cds`.
`-J` drops models without a start or stop codon, or with an internal stop.

***P. minor*:** **44,960 proteins**; 5 incomplete models were dropped. For comparison, MAKER v2 lost about 7,600 partial models at this step.

## 3. `mmseqs2` + `best_hits`: homology

For each database: `mmseqs easy-search ... --start-sens 2 -s 7 --sens-steps 3` (iterative sensitivity), with **a private tmp directory per database**. That tmp directory is the v3 fix for the `/tmp` race; see [TROUBLESHOOTING #18](TROUBLESHOOTING.md#18-could-not-create-symlink-of-tmp-in-rule-mmseqs2). `best_hits` keeps the best-scoring hit per protein with e-value < 1e-5.

Outputs: `mmseqs_output/<SP>_<db>.ofmt6` (all hits) and `<SP>_<db>_best_hit.out`.

| Database | *P. minor* proteins with a best hit | Runtime |
|---|---|---|
| UniRef50 | **43,580 / 44,960 (96.9 %)** | ~25 min |
| herbicide-resistance KO set | **429** | ~1 min |
| NCBI GenomesDatabase v1 | **43,528 (96.8 %)** (11,470 distinct subjects) | ~52 min (20 threads) |

## 4. `interproscan`: domains, GO, pathways

InterProScan 5.77-108.0 with `--goterms --pathways --iprlookup`, run in **chunks of 2,000 proteins** (avoids the LevelDB out-of-memory crash; chunks that are already done are skipped, so the step can resume) → `iprscan_output/<SP>.iprscan.tsv`.

***P. minor*:** 23 chunks at about 1.5 h each (30 CPUs), so **~34 h**. This is the bottleneck of the functional module. Result: **39,362 proteins (87.5 %) with InterPro entries** and **31,545 (70.2 %) with GO terms**.

## 5. `prot_list` + `multiloc`: subcellular localization

`prot_list` collects the InterProScan GO terms per protein. MultiLoc2 (plant model, using those GO terms) then runs in about 20 parallel splits → `multiloc_outputs/<SP>_ml2.txt`.

***P. minor*:** a prediction for 44,949 proteins (100 %), in **~9.5 h** after InterProScan finished.

## 6. `hmmer`: Pfam domains

`hmmscan --domtblout` against Pfam-A (HMMER 3.3) → `hmmer_outputs/<SP>.domtblout`.

***P. minor*:** 41,650 proteins with at least one Pfam hit before thresholds are applied (the merge step keeps the significant ones). Runtime about 7.5 h on 4 CPUs.

## 7. `sigtarp`: signal peptides and targeting

SignalP 5.0b (`-org euk`) and TargetP 2.0 (`-org pl`) → `sigtarp_outputs/<SP>.gff3`, `<SP>_summary.targetp2`.

***P. minor*:**
- SignalP: **4,986** signal peptides.
- TargetP: 34,336 noTP, 5,621 SP (secretory), 3,372 cTP (chloroplast), 1,522 mTP (mitochondrion), 109 luTP (thylakoid lumen).

## 8. `generate_tsv` → `order_tsv`: merged table

`functional_merge.py` merges every result into one row per protein. It writes `functional_annotation.alphabetical.tsv`; `order_tsv` then reorders it into **assembly order** (the order of `<SP>.aa`, which follows the GFF), giving `functional_annotation.tsv`.

Columns: `ID sequence`, `NCBI DESCRIPTION`, `NCBI Subject;E-value;Bit`, `herbres ...`, `uniref ...`, `SignalP Pos;Pr`, `TargetP ...`, `MULTILOC`, `IPRSCAN GO`, `IPRSCAN IPR`, `Hmmer Pfam`.

***P. minor* coverage** (file [`results/PhaMi/functional/coverage.tsv`](../results/PhaMi/functional/coverage.tsv)), compared with the MAKER v2 annotation run through the same pipeline:

| Column | v3 (44,960 proteins) | MAKER v2 (47,494 proteins) |
|---|---|---|
| NCBI hit | 43,528 (96.8 %) | 46,202 (97.3 %) |
| NCBI description (not "Description not found") | 17,298 (38.5 %) | 17,995 (37.9 %) |
| herbicide-resistance hit | 429 | 291 |
| UniRef50 hit | 43,580 (96.9 %) | 46,401 (97.7 %) |
| SignalP signal peptide | 4,986 (11.1 %) | 4,419 (9.3 %) |
| MultiLoc2 | 44,949 (100 %) | 47,494 (100 %) |
| InterProScan GO | 31,545 (70.2 %) | 33,046 (69.6 %) |
| InterProScan IPR | 39,362 (87.5 %) | 42,301 (89.1 %) |

Notes on these columns:
- The `Hmmer Pfam` column is filled for every protein in both versions: `functional_merge.py` records the hmmscan hits without an e-value cut, so use InterPro or the Pfam terms in the functional GFF for strict domain calls.
- Many NCBI subjects are gene-model IDs from the GenomesDatabase without a text description, which is why only about 38 % of proteins have a description.

## 9. `agat_functional_gff`: functional GFF

1. `extract_hits.py` writes the NCBI subjects of the best hits, repairing the glued headers of the NCBI FASTA on the fly.
2. `agat_sp_manage_functional_annotation.pl` adds product names (from the NCBI hits) and InterPro / GO / Pfam / PANTHER / … terms (from InterProScan without the pathways column) to the longest-isoform GFF → `functional_outputs/<SP>.functional.gff`.

The per-database term lists are written to `functional_outputs/agat_<SP>/` (`GO.txt`, `Pfam.txt`, `InterPro.txt`, …), together with AGAT's `report.txt`.

***P. minor*:** 11,470 NCBI subjects extracted (all found), and `PhaMi.functional.gff` written in about 3 min. AGAT prints `Use of uninitialized value $desc` for database headers without a description; this is harmless. AGAT's report: [`results/PhaMi/functional/agat_report.txt`](../results/PhaMi/functional/agat_report.txt).

---

## Runtime (*P. minor*, 44,960 proteins, 1 node, 64 CPUs, 250 GB)

| Step | Wall-clock |
|---|---|
| agat + gffread | ~3 min |
| MMseqs2 (3 databases in parallel) | ~1 h (NCBI is the slowest) |
| HMMER | ~7.5 h (in parallel with InterProScan) |
| SignalP + TargetP | ~1.7 h (starts after HMMER; in parallel with InterProScan) |
| InterProScan | ~34 h (bottleneck) |
| MultiLoc2 | ~9.5 h (after InterProScan) |
| merge + order + AGAT functional GFF | ~5 min |
| **Total** | **~44 h** wall-clock (2026-10-20 12:59 → 10-22 08:24, plus the ~5 min finishing run) |
