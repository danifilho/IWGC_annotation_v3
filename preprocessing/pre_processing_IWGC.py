#!/usr/bin/env python3
"""IWGC pre-processing: prepares the RNA-seq input of the structural pipeline.
(Genome masking, EDTA + bedtools, is part of the structural Snakefile.)

  trim      Illumina RNA-seq pairs: FastQC (raw) -> fastp (adapter, poly-G,
            minimum length) -> FastQC (trimmed) -> MultiQC, plus a per-library
            summary table. Output <out>/trimmed/<lib>_R1.fq.gz / _R2.fq.gz is the
            RNA-seq input of predictor: braker4 (BRAKER4 aligns it with HISAT2).

Each tool runs in a Singularity/Apptainer image from --images (fastqc.sif,
fastp.sif, multiqc.sif). Paths are resolved to their real
location and those directories are bound, so symlinked paths such as
/mnt/research/... also work inside the containers. Run it inside a SLURM job
(see pre_processing_IWGC.sb); trim uses --jobs x --threads CPUs.

Examples:
  pre_processing_IWGC.py trim --reads-dir /path/illumina_data \\
      --r1 '{lib}-pairsR1.fq.gz' --r2 '{lib}-pairsR2.fq.gz' \\
      --adapter CTGTCTCTTATACACATCT --out rnaseq --name AmaPa \\
      --images images --jobs 5 --threads 16
"""
import argparse, json, os, re, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

def real(p):
    return os.path.realpath(p)

def container_cmd(a, image, binds):
    exe = a.container or shutil.which("singularity") or shutil.which("apptainer")
    if not exe:
        sys.exit("singularity/apptainer not found on PATH")
    img = real(os.path.join(a.images, image))
    if not os.path.isfile(img):
        sys.exit(f"missing image {img}")
    b = []
    for d in sorted(set(binds)):
        b += ["-B", d]
    return [exe, "exec"] + b + [img]

def run(cmd, log):
    log.write("$ " + " ".join(cmd) + "\n")
    log.flush()
    subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)


# ------------------------------------------------------------------ trim
def find_libraries(reads_dir, r1):
    pat = re.compile("^" + re.escape(r1).replace(re.escape("{lib}"), "(.+)") + "$")
    return sorted(m.group(1) for f in os.listdir(reads_dir) if (m := pat.match(f)))

def trim_one(lib, a, out, binds):
    r1 = real(os.path.join(a.reads_dir, a.r1.format(lib=lib)))
    r2 = real(os.path.join(a.reads_dir, a.r2.format(lib=lib)))
    done = f"{out}/logs/{lib}.done"
    if os.path.exists(done):
        return lib, "already done"
    t0 = time.time()
    # R1/R2 names: BRAKER4 names each library after its R1 file without "_R1"
    t1, t2 = f"{out}/trimmed/{lib}_R1.fq.gz", f"{out}/trimmed/{lib}_R2.fq.gz"
    adapter = (["--adapter_sequence", a.adapter, "--adapter_sequence_r2", a.adapter]
               if a.adapter else ["--detect_adapter_for_pe"])
    fastqc = container_cmd(a, "fastqc.sif", binds)
    fastp = container_cmd(a, "fastp.sif", binds)
    with open(f"{out}/logs/{lib}.log", "w") as log:
        run(fastqc + ["fastqc", "-t", "2", "-o", f"{out}/fastqc_raw", r1, r2], log)
        run(fastp + ["fastp", "-w", str(a.threads), "-i", r1, "-I", r2, "-o", t1, "-O", t2]
            + adapter + a.fastp_extra.split()
            + ["-j", f"{out}/fastp/{lib}.fastp.json", "-h", f"{out}/fastp/{lib}.fastp.html"], log)
        run(fastqc + ["fastqc", "-t", "2", "-o", f"{out}/fastqc_trimmed", t1, t2], log)
    open(done, "w").write(f"{time.time() - t0:.0f} s\n")
    return lib, f"{(time.time() - t0) / 60:.1f} min"

def fastp_summary(out, libs):
    rows = ["library\tpairs_raw\tpairs_trimmed\tpct_kept\treads_with_adapter\tpct_with_adapter"
            "\tq30_trimmed\tduplication"]
    tot_in = tot_out = 0
    for lib in libs:
        j = json.load(open(f"{out}/fastp/{lib}.fastp.json"))
        s = j["summary"]
        pin = s["before_filtering"]["total_reads"] // 2
        pout = s["after_filtering"]["total_reads"] // 2
        ad = j.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0)
        tot_in += pin; tot_out += pout
        rows.append(f"{lib}\t{pin}\t{pout}\t{100 * pout / pin:.1f}\t{ad}\t{100 * ad / (2 * pin):.1f}"
                    f"\t{100 * s['after_filtering']['q30_rate']:.1f}\t{100 * j['duplication']['rate']:.1f}")
    rows.append(f"TOTAL\t{tot_in}\t{tot_out}\t{100 * tot_out / tot_in:.1f}\t\t\t\t")
    open(f"{out}/trimming_summary.tsv", "w").write("\n".join(rows) + "\n")
    return tot_in, tot_out

def cmd_trim(a):
    out = real(a.out)
    for d in ("trimmed", "fastp", "fastqc_raw", "fastqc_trimmed", "multiqc", "logs"):
        os.makedirs(f"{out}/{d}", exist_ok=True)
    libs = a.libraries or find_libraries(a.reads_dir, a.r1)
    if not libs:
        sys.exit(f"no file in {a.reads_dir} matches {a.r1}")
    binds = [real(a.reads_dir), out]
    print(f"{len(libs)} libraries, {a.jobs} at a time x {a.threads} fastp threads: {' '.join(libs)}", flush=True)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for lib, msg in ex.map(lambda l: trim_one(l, a, out, binds), libs):
            print(f"  {lib}: {msg}", flush=True)
    multiqc = container_cmd(a, "multiqc.sif", binds)
    with open(f"{out}/logs/multiqc.log", "w") as log:
        run(multiqc + ["multiqc", "-f", "-o", f"{out}/multiqc", "-n", f"{a.name}_raw", f"{out}/fastqc_raw"], log)
        run(multiqc + ["multiqc", "-f", "-o", f"{out}/multiqc", "-n", f"{a.name}_trimmed",
                       f"{out}/fastqc_trimmed", f"{out}/fastp"], log)
    pin, pout = fastp_summary(out, libs)
    print(f"pairs: {pin:,} raw -> {pout:,} trimmed ({100 * pout / pin:.1f} %)")
    print(f"summary: {out}/trimming_summary.tsv  reports: {out}/multiqc/{a.name}_raw.html, {a.name}_trimmed.html")
    print(f"BRAKER4 input: rnaseq_dir: {out}/trimmed  rnaseq_r1: \"{{lib}}_R1.fq.gz\"  rnaseq_r2: \"{{lib}}_R2.fq.gz\"")
    print(f"elapsed: {(time.time() - t0) / 60:.1f} min")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--images", required=True, help="directory with the .sif images")
    common.add_argument("--container", help="singularity or apptainer executable (default: found on PATH)")

    t = sub.add_parser("trim", parents=[common], help="FastQC + fastp + MultiQC on paired RNA-seq")
    t.add_argument("--reads-dir", required=True)
    t.add_argument("--r1", default="{lib}_R1.fastq.gz", help="R1 file name, {lib} = library (default %(default)s)")
    t.add_argument("--r2", default="{lib}_R2.fastq.gz", help="R2 file name (default %(default)s)")
    t.add_argument("--libraries", nargs="*", help="library names (default: every file matching --r1)")
    t.add_argument("--adapter", default="", help="adapter for both reads, e.g. Nextera CTGTCTCTTATACACATCT, "
                                                 "TruSeq AGATCGGAAGAGC (default: fastp detects it)")
    t.add_argument("--fastp-extra", default="--trim_poly_g --length_required 36", help="default: %(default)s")
    t.add_argument("--out", required=True)
    t.add_argument("--name", default="rnaseq", help="MultiQC report prefix (default %(default)s)")
    t.add_argument("--jobs", type=int, default=4, help="libraries in parallel (default %(default)s)")
    t.add_argument("--threads", type=int, default=16, help="fastp threads per library, max 16 (default %(default)s)")
    t.set_defaults(func=cmd_trim)

    a = ap.parse_args()
    a.func(a)

if __name__ == "__main__":
    main()
