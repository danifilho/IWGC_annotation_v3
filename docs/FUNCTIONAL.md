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
   input GFF3 + functional TSV + InterProScan ─► functional_gff ─► <SP>.functional.gff3 (IWGC layout)
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

For each database: `mmseqs easy-search ... --start-sens 2 -s 7 --sens-steps 3` (iterative sensitivity), with **a private tmp directory per database**. That tmp directory is the v3 fix for the `/tmp` race; see [TROUBLESHOOTING #18](TROUBLESHOOTING.md#18-could-not-create-symlink-of-tmp-in-rule-mmseqs2). Before searching a FASTA target, the rule checks that it starts with `>` and has no glued headers ([#28](TROUBLESHOOTING.md#28-most-ncbi-best-hits-are-the-same-id-aa049522-ra)).

`best_hits` keeps one hit per protein with e-value < 1e-5:
- **herbres and UniRef50:** the highest-bitscore hit.
- **NCBI**, with `ncbi_prefer_named_hits: true` (the default): the highest-bitscore hit **whose FASTA header has a description**. The unnamed best hit is used only when no named hit passes the cut.

Why the NCBI rule? `GenomesDatabase.v1` contains 49,542 proteins with bare `AA…-RA` headers and no name or species. They are an oat proteome added from CoGe for IWGC v1; `SalTr_v01.0.gff` labels them `[Avena sativa]`. Oat is the closest relative of grass weeds in the database, so these unnamed proteins are the plain best hit of most grass genes: 27,214 of 44,960 *P. minor* v3 proteins. That leaves only 36 % of proteins with a product name. `ncbi_unnamed_subjects` lists the unnamed subjects directly from the FASTA headers.

Outputs: `mmseqs_output/<SP>_<db>.ofmt6` (all hits), `<SP>_<db>_best_hit.out`, `ncbi_unnamed_subjects.txt`.

| Database | *P. minor* v3 proteins with a best hit | Runtime |
|---|---|---|
| UniRef50 | **43,580 / 44,960 (96.9 %)** | ~25 min |
| herbicide-resistance KO set | **429** | ~1 min |
| NCBI GenomesDatabase v1 (repaired, named hits) | **43,566 (96.9 %)**, with a product name for **43,293 (96.3 %)**; only 273 keep an unnamed `AA` hit | < 10 min (20 threads; ~52 min against the unrepaired DB, probably slowed by aligning against its 18.9 M aa record) |

Effect of the NCBI rule on *P. minor* v3:
- 26,941 proteins moved from an unnamed oat hit to a named RefSeq hit.
- Identity to the hit dropped by 5.5 points on average, and the bitscore kept 92 % of its value.
- The named hits are mostly *Triticum dicoccoides* (20,077), *Aegilops tauschii* (11,339) and *Brachypodium distachyon* (8,664).

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
| NCBI hit | 43,566 (96.9 %) | 46,238 (97.4 %) |
| NCBI description (not "Description not found") | 43,293 (96.3 %) | 45,938 (96.7 %) |
| herbicide-resistance hit | 429 | 291 |
| UniRef50 hit | 43,580 (96.9 %) | 46,401 (97.7 %) |
| SignalP signal peptide | 4,986 (11.1 %) | 4,419 (9.3 %) |
| MultiLoc2 | 44,949 (100 %) | 47,494 (100 %) |
| InterProScan GO | 31,545 (70.2 %) | 33,046 (69.6 %) |
| InterProScan IPR | 39,362 (87.5 %) | 42,301 (89.1 %) |

Notes on these columns:
- The `Hmmer Pfam` column is filled for every protein in both versions: `functional_merge.py` records the hmmscan hits without an e-value cut, so use the InterPro column for strict domain calls.
- Before the database repair and the named-hit rule (§3), only about 38 % of proteins had a description, because most best hits were the unnamed `AA` record.

## 9. `functional_gff`: IWGC-format functional GFF3

`scripts/iwgc_functional_gff.py` writes the release GFF3 in the **same layout as the IWGC v1 annotations** (e.g. `Weed_lab/SalTr_v01.0.gff`):
- It reads the input structural GFF3 (assembly order, every feature kept), `functional_annotation.tsv` and the InterProScan TSV.
- Each **gene** gets:
  - `Dbxref=<NCBI hit>:<GeneID>:<bitscore>`
  - `Ontology_term=<GO list>`
  - `Note=<NCBI description>:(<db>|<accession>|<name>)`, or `Note=Unknown Protein` when there is no NCBI hit. tRNA genes and models dropped by `gffread -J` therefore show `Unknown Protein`, as in the v1 releases.
- mRNA, exon, CDS, UTR and tRNA lines keep only `ID` and `Parent`.
- Column 2 is set to `gff_source` (e.g. `IWGC_v3`), and the header records the species, version and tools.

The original IWGC v1 pipeline ([PattersonWeedLab/IWGC_annotation_pipeline](https://github.com/PattersonWeedLab/IWGC_annotation_pipeline)) does not publish the script that wrote these files, so the layout was reproduced from `SalTr_v01.0.gff` itself.

**GeneID.** In the v1 files, the middle `Dbxref` field is the NCBI GeneID; the v1 NCBI set was built with NCBI Datasets, which carry GeneIDs. The current GenomesDatabase FASTA has no GeneIDs, so `scripts/make_geneid_map.sh` builds them. It downloads NCBI `gene2refseq` (~2.4 GB), keeps only the best-hit accessions and deletes the download:
```bash
bash scripts/make_geneid_map.sh inputs/ncbi_protein2geneid.tsv mmseqs_output/<SP>_ncbi_best_hit.out
# then set  ncbi_geneid_map: inputs/ncbi_protein2geneid.tsv  in config.yaml
```
Without a map, the script falls back to the `LOC<number>` in the description when there is one, and otherwise writes `NA`. Two kinds of hits never get a GeneID:
- **IWGC gene models** (`AA…` IDs), which are part of the GenomesDatabase but not NCBI records;
- **retired accession versions** (e.g. `.1` after NCBI moved to `.2`).

***P. minor* (runtime ~6 s):**

| | v3 (BRAKER3) | v2 (MAKER) |
|---|---|---|
| Genes | 44,965 | 57,163 (incl. 2,039 tRNA) |
| With NCBI hit (`Dbxref`) | 43,566 (96.9 %) | 46,238 (80.9 %; tRNA and partial models have none) |
| `Note=Unknown Protein` | 1,399 | 10,925 |
| With GO (`Ontology_term`) | 31,545 | 33,046 |
| `Note` with an InterProScan match | 42,707 | 45,608 |
| `Dbxref` with GeneID (gene2refseq, 2026-10-22) | 42,164 | 44,662 |
| `Dbxref` GeneID = `NA` (unnamed `AA` hits + retired accession versions) | 1,402 | 1,576 |

> This step replaces `agat_sp_manage_functional_annotation.pl`. AGAT expects UniProt-style FASTA headers and could not read the NCBI descriptions, so every gene came out as `product=hypothetical protein`, in v2 and v3 alike. See [TROUBLESHOOTING #27](TROUBLESHOOTING.md#27-agat-functional-gff-every-gene-is-hypothetical-protein).

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
| merge + order + functional GFF3 | ~2 min |
| **Total** | **~44 h** wall-clock (2026-10-20 12:59 → 10-22 08:24, plus the ~5 min finishing run) |
