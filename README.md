# IWGC Genome Annotation Pipeline v3

Structural and functional genome annotation for the **International Weed Genomics Consortium (IWGC)**, built with Snakemake and Singularity and developed on MSU ICER HPCC (SLURM).

| Module | What it does | Main outputs |
|---|---|---|
| **Structural** | repeat soft-masking → **BRAKER3** with proteins + **PacBio Iso-Seq** reads (one genome-wide model) → GFF3 with IWGC gene IDs → BUSCO | `<SP>.v3.gff3`, proteins, CDS, BUSCO |
| **Functional** | MMseqs2 (NCBI, UniRef50, herbicide resistance), InterProScan, Pfam, SignalP/TargetP, MultiLoc2 | `functional_annotation.tsv`, IWGC-format `<SP>.functional.gff3` |

On *Phalaris minor* (4.9 Gb, allotetraploid), BRAKER3 reached the same BUSCO completeness as the earlier MAKER pipeline (**95.9 % vs 96.1 %**) in **~11 h instead of ~5 days**.

📖 **Full documentation: [Wiki](https://github.com/danifilho/IWGC_annotation_v3/wiki)**. It covers every step, every config key, the [functional GFF3 format](https://github.com/danifilho/IWGC_annotation_v3/wiki/Functional-GFF3), [troubleshooting](https://github.com/danifilho/IWGC_annotation_v3/wiki/Troubleshooting) and the [*P. minor* results](https://github.com/danifilho/IWGC_annotation_v3/wiki/Phalaris-minor-Case-Study).

---

## 1. What you need

- **Software:** Snakemake 9.13.7 and SingularityCE 4.1.2 (or Apptainer). SLURM is optional.
  ```bash
  conda create -n iwgc -c conda-forge -c bioconda snakemake=9.13.7 biopython pandas
  ```
- **Hardware** (4.9 Gb genome):
  - structural: 1 node, 64 CPUs, 400 GB RAM;
  - functional: 1 node, 64 CPUs, 250 GB RAM;
  - plenty of scratch **file** quota, since GeneMark can create ~350 k temporary files.
- **Structural inputs:**
  - genome FASTA (simple headers, no spaces);
  - proteins of a close relative **and** OrthoDB Viridiplantae;
  - PacBio Iso-Seq **reads** (FLNC `.bam`, or FASTA/FASTQ), not collapsed transcripts;
  - a local BUSCO lineage.
- **Functional databases:** the IWGC NCBI protein set (`GenomesDatabase.v1`), UniRef50 (MMseqs2), the herbicide-resistance KO set, Pfam-A and InterProScan data 5.77-108.0.

## 2. Get the code and the images

```bash
git clone https://github.com/danifilho/IWGC_annotation_v3.git
mkdir -p $SCRATCH/MySp && cd $SCRATCH/MySp
cp -r ~/IWGC_annotation_v3/structural structural
cp -r ~/IWGC_annotation_v3/functional functional
```

Container images go in each workdir's `images/` **with these file names**:

```bash
# structural/images
cd structural && mkdir -p images && cd images
singularity pull braker3_isoseq.sif docker://teambraker/braker3:isoseq      # BRAKER 3.0.8 (Iso-Seq variant)
singularity pull braker3.sif        docker://teambraker/braker3:latest      # only for protein-only runs
singularity pull minimap2.sif       docker://staphb/minimap2:2.26
singularity pull samtools.sif       docker://staphb/samtools:1.19
singularity pull busco.sif          docker://ezlabgva/busco:v5.8.2_cv1
singularity pull dfam_tetools.sif   docker://dfam/tetools                   # RepeatModeler/RepeatMasker
singularity pull -U bedtools.sif    library://danifilho/structural_annotation/bedtools:latest
singularity pull -U gffread.sif     library://danifilho/functional_annotation_images/gffread:latest

# functional/images
cd ../../functional && mkdir -p images && cd images
for i in agat gffread hmmer3 iprscan mmseqs2 samtools sigtarp; do
  singularity pull -U $i.sif library://danifilho/functional_annotation_images/$i:latest
done
singularity pull -U multiloc2_v3.img library://danifilho/functional_annotation_images/multiloc2:latest
```

Verified tool versions and database downloads are listed in the [Installation and setup](https://github.com/danifilho/IWGC_annotation_v3/wiki/Installation-and-Setup) page.

## 3. Prepare inputs

```
structural/inputs/    genome.fasta, protein FASTA files, Iso-Seq reads
structural/busco_downloads/lineages/<lineage>/
functional/inputs/    genome.fasta + the structural GFF3 (after step 5)
functional/databases/ current/ herb_res_prots/ mmseqs_db/ Pfam-A.hmm* interproscan-5.77-108.0/
```

> ⚠️ Every file the containers read must be a **real file or a hard link inside the workdir**. A symlink pointing outside it fails inside the container.

> ⚠️ **Repair the NCBI database once** before the functional run. The original `GenomesDatabase.v1.prot` merges 49,542 proteins into one record ([why](https://github.com/danifilho/IWGC_annotation_v3/wiki/Troubleshooting#28-most-ncbi-best-hits-are-the-same-id-aa049522-ra)):
> ```bash
> python3 functional/scripts/fix_fasta_headers.py \
>   functional/databases/current/GenomesDatabase.v1.prot functional/databases/current/GenomesDatabase.v1.fixed.prot
> ```

## 4. Configure

Edit `structural/config.yaml`. The keys you must set:

```yaml
volume_name: /abs/path/to/structural          # this workdir
species_abbreviation: MySp                    # gene-ID prefix (MySpChr01g000010)
fasta_file: genome.fasta
softmasked_genome: ""                         # or a pre-masked genome to skip masking
protein_files: [relative_proteins.fa, Viridiplantae.fa]
isoseq_reads: sample.flnc.bam                 # "" = protein-only BRAKER
braker_species: MySp_v3                       # new AUGUSTUS model name
busco_lineage: busco_downloads/lineages/poales_odb10
```

Edit `functional/config.yaml`:

```yaml
volume_name: /abs/path/to/functional
fasta_file: genome.fasta
final_filtering_gff: MySp.v3.gff3             # from the structural run
species_name: MySp
# the 'ncbi' entry must point to the REPAIRED FASTA (path and ncbidata)
```

All keys and their defaults: [Configuration](https://github.com/danifilho/IWGC_annotation_v3/wiki/Configuration).

## 5. Run

```bash
# structural
cd structural
snakemake -n                    # dry run
sbatch run_structural.sb        # ~11 h for 4.9 Gb
grep "will be executed in" braker_outputs/braker.log   # must say: ETP mode (BRAKER3)

# functional
cp final_outputs/MySp.v3.gff3 ../functional/inputs/
ln inputs/genome.fasta ../functional/inputs/
cd ../functional
snakemake -n
sbatch run_functional.sb        # ~2 days (InterProScan is the bottleneck)

# after each run: --keep-going does not stop on errors, so check
grep -n "Error in rule" *.err
```

Optional, needs internet: fill NCBI GeneIDs in the functional GFF3 with `scripts/make_geneid_map.sh` (see [Installation and setup](https://github.com/danifilho/IWGC_annotation_v3/wiki/Installation-and-Setup#ncbi-geneids-optional)).

## 6. Main outputs

| File | Content |
|---|---|
| `structural/final_outputs/<SP>.v3.gff3` | gene models with IWGC IDs |
| `structural/final_outputs/<SP>.v3.{proteins,cds,longest.proteins}.fa` | sequences |
| `structural/busco_outputs/<SP>.v3.busco_summary.txt` | BUSCO |
| `functional/functional_outputs/functional_annotation.tsv` | one row per protein: position, NCBI / UniRef50 / herbicide-resistance hits, InterPro, GO, Pfam, SignalP/TargetP, MultiLoc2 |
| `functional/functional_outputs/<SP>.functional.gff3` | **release GFF3**: structure + `Dbxref`, `Ontology_term`, `Note` per gene ([format](https://github.com/danifilho/IWGC_annotation_v3/wiki/Functional-GFF3)) |

Every output: [Running and outputs](https://github.com/danifilho/IWGC_annotation_v3/wiki/Running-and-Outputs). Example results for *P. minor*: [`results/PhaMi/`](results/PhaMi/).

## 7. Repository layout

```
structural/   Snakefile, config.yaml, run_structural.sb, scripts/
functional/   Snakefile, config.yaml, run_functional.sb, scripts/
results/PhaMi/   example outputs (GFF3, TSV, BUSCO, comparisons)
```

## Citation and contact

Please cite the underlying tools ([list](https://github.com/danifilho/IWGC_annotation_v3/wiki#citations)).
Contact: **dasilvaf@msu.edu** · [@danifilho](https://github.com/danifilho)
