#!/usr/bin/env python3
"""Write samples.csv + config.ini for one BRAKER4 run (pre-masked genome,
proteins, paired-end RNA-seq, optional Iso-Seq reads) into OUTDIR.

samples.csv: the genome goes in both `genome` and `genome_masked`, so BRAKER4
does not mask it again. RNA-seq pairs are joined with ':' (BRAKER4 aligns
them with HISAT2 and names each library after its R1 file without _R1).
Iso-Seq reads (FASTA/FASTQ) go in isoseq_fastq: with RNA-seq this is BRAKER4's
dual mode (minimap2 inside BRAKER4; one GeneMark-ETP run with both,
StringTie --mix of short + long reads).

config.ini: BRAKER4's config.ini.example with SECTION.key=value overrides
(comments are not kept).

Usage: braker4_setup.py --outdir D --sample SP --genome G --proteins P
         --busco-lineage L --template config.ini.example
         --r1 a_R1.fq.gz ... --r2 a_R2.fq.gz ... [--isoseq reads.fa.gz ...]
         [--set SECTION.key=value ...]
"""
import argparse, configparser, csv, os, sys

COLUMNS = ["sample_name", "genome", "genome_masked", "protein_fasta", "bam_files",
           "fastq_r1", "fastq_r2", "sra_ids", "varus_genus", "varus_species",
           "isoseq_bam", "isoseq_fastq", "busco_lineage", "reference_gtf"]

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("--outdir", "--sample", "--genome", "--proteins", "--busco-lineage", "--template"):
        ap.add_argument(a, required=True)
    ap.add_argument("--r1", nargs="+", required=True)
    ap.add_argument("--r2", nargs="+", required=True)
    ap.add_argument("--isoseq", nargs="*", default=[])
    ap.add_argument("--set", nargs="*", default=[])
    a = ap.parse_args()
    if len(a.r1) != len(a.r2):
        sys.exit(f"{len(a.r1)} R1 files but {len(a.r2)} R2 files")

    row = dict.fromkeys(COLUMNS, "")
    row.update(sample_name=a.sample, genome=a.genome, genome_masked=a.genome,
               protein_fasta=a.proteins, fastq_r1=":".join(a.r1), fastq_r2=":".join(a.r2),
               isoseq_fastq=":".join(a.isoseq), busco_lineage=a.busco_lineage)
    with open(os.path.join(a.outdir, "samples.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerow(row)

    ini = configparser.ConfigParser(interpolation=None)
    ini.optionxform = str
    ini.read(a.template)
    for s in a.set:
        key, value = s.split("=", 1)
        section, option = key.split(".", 1)
        if not ini.has_section(section):
            ini.add_section(section)
        ini.set(section, option, value)
    with open(os.path.join(a.outdir, "config.ini"), "w") as fh:
        ini.write(fh)
    print(f"BRAKER4 inputs: {len(a.r1)} RNA-seq libraries, {len(a.isoseq)} Iso-Seq files, overrides: {' '.join(a.set) or 'none'}")

if __name__ == "__main__":
    main()
