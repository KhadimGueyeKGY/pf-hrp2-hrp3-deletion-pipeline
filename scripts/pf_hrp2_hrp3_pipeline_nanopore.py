#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import gzip
import json
import logging
import os
import re
import shutil
import statistics
import subprocess
import sys
import traceback
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam


AUTHOR = "Khadim Gueye"

RESET = "\033[0m"
BOLD = "\033[1m"
RED = "\033[31m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
BLUE = "\033[34m"
MAGENTA = "\033[35m"
CYAN = "\033[36m"


def color(tag: str, c: str, msg: str) -> str:
    return f"{c}{BOLD}{tag}{RESET} {msg}"


def info(msg: str) -> str:
    return color("[INFO]", CYAN, msg)


def step(msg: str) -> str:
    return color("[STEP]", BLUE, msg)


def done(msg: str) -> str:
    return color("[DONE]", GREEN, msg)


def warn(msg: str) -> str:
    return color("[WARN]", YELLOW, msg)


def fail(msg: str) -> str:
    return color("[FAIL]", RED, msg)


class PipelineError(RuntimeError):
    pass


@dataclass(frozen=True)
class BarcodeSample:
    state: str
    barcode: str
    barcode_dir: Path
    sample: str


@dataclass(frozen=True)
class GeneRegion:
    chrom: str
    start0: int
    end: int
    name: str
    strand: str

    @property
    def start1(self) -> int:
        return self.start0 + 1

    @property
    def length(self) -> int:
        return self.end - self.start0


@dataclass(frozen=True)
class Config:
    demux_dir: Path
    output_dir: Path
    reference: Path
    genes_bed: Path
    metadata: Optional[Path]
    threads: int
    sample_workers: int
    min_length: int
    min_qscore: int
    min_mapq: int
    min_baseq: int
    min_site_dp: int
    min_site_qual: int
    gene_presence_min_mean_depth: float
    gene_presence_min_breadth_1x: float
    gene_presence_min_breadth_5x: float
    gene_presence_min_breadth_10x: float
    fastqc_bin: str
    multiqc_bin: str
    chopper_bin: str
    minimap2_bin: str
    samtools_bin: str
    bcftools_bin: str
    freebayes_bin: str
    bgzip_bin: str
    tabix_bin: str
    clair3_bin: str
    clair3_model_path: Optional[Path]
    clair3_platform: str
    clair3_include_all_ctgs: bool
    run_clair3: bool
    run_freebayes: bool
    run_bcftools: bool
    freebayes_ploidy: int
    bcftools_ploidy: int
    snpeff_bin: Path
    snpeff_db: str
    snpeff_config: Optional[Path]
    run_deploid: bool
    deploid_bin: str
    population_vcf: Optional[Path]
    genetic_map: Optional[Path]
    exclude_vcf: Optional[Path]
    hap_min_reads: int
    hap_min_informative_sites: int
    hap_min_fraction: float
    resume: bool
    keep_tmp: bool
    log_file: Path


def setup_logger(log_file: Path) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("pf_hrp_ont")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fh = logging.FileHandler(log_file, mode="a")
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s"))
    logger.addHandler(fh)
    logger.propagate = False
    return logger


def get_logger(name: str) -> logging.Logger:
    lg = logging.getLogger(f"pf_hrp_ont.{name}")
    lg.setLevel(logging.INFO)
    lg.propagate = True
    return lg


def q(x: Path | str) -> str:
    s = str(x)
    return "'" + s.replace("'", "'\"'\"'") + "'"


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def ensure_parent(p: Path) -> None:
    p.parent.mkdir(parents=True, exist_ok=True)


def nonempty(p: Path) -> bool:
    return p.exists() and p.stat().st_size > 0


def which_or_raise(tool: str) -> str:
    p = shutil.which(tool)
    if not p:
        raise PipelineError(f"Missing dependency in PATH: {tool}")
    return p


def run_cmd(cmd: List[str], logger: logging.Logger, cwd: Optional[Path] = None, stdout_path: Optional[Path] = None) -> None:
    logger.info("CMD: %s", " ".join(map(str, cmd)))
    if stdout_path:
        ensure_parent(stdout_path)
        with open(stdout_path, "wb") as out:
            p = subprocess.run(cmd, cwd=str(cwd) if cwd else None, stdout=out, stderr=subprocess.PIPE, check=False)
        if p.stderr:
            logger.info("STDERR:\n%s", p.stderr.decode(errors="replace"))
    else:
        p = subprocess.run(cmd, cwd=str(cwd) if cwd else None, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        if p.stdout:
            logger.info("STDOUT:\n%s", p.stdout)
        if p.stderr:
            logger.info("STDERR:\n%s", p.stderr)
    if p.returncode != 0:
        raise PipelineError(f"Command failed: {' '.join(map(str, cmd))}")


def run_bash(script: str, logger: logging.Logger, cwd: Optional[Path] = None) -> None:
    logger.info("BASH: %s", script)
    p = subprocess.run(["bash", "-lc", script], cwd=str(cwd) if cwd else None, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.stdout:
        logger.info("STDOUT:\n%s", p.stdout)
    if p.stderr:
        logger.info("STDERR:\n%s", p.stderr)
    if p.returncode != 0:
        raise PipelineError(f"Bash failed: {script}")


def safe_unlink(p: Path) -> None:
    try:
        p.unlink(missing_ok=True)
    except Exception:
        pass


def check_dependencies(cfg: Config) -> None:
    required = [
        cfg.fastqc_bin,
        cfg.multiqc_bin,
        cfg.chopper_bin,
        cfg.minimap2_bin,
        cfg.samtools_bin,
        cfg.bcftools_bin,
        cfg.bgzip_bin,
        cfg.tabix_bin,
    ]
    if cfg.run_freebayes:
        required.append(cfg.freebayes_bin)
    if cfg.run_clair3:
        required.append(cfg.clair3_bin)
    if cfg.run_deploid:
        required.append(cfg.deploid_bin)
    for tool in required:
        which_or_raise(tool)
    if not cfg.snpeff_bin.exists() and not shutil.which(str(cfg.snpeff_bin)):
        raise PipelineError(f"snpEff not found: {cfg.snpeff_bin}")


def index_reference(ref: Path, cfg: Config, logger: logging.Logger) -> None:
    if not Path(str(ref) + ".fai").exists():
        run_cmd([cfg.samtools_bin, "faidx", str(ref)], logger)


def fai_contigs(fai: Path) -> Dict[str, int]:
    out: Dict[str, int] = {}
    with open(fai) as fh:
        for line in fh:
            if line.strip():
                t = line.rstrip("\n").split("\t")
                out[t[0]] = int(t[1])
    return out


def parse_bed(bed: Path) -> List[GeneRegion]:
    genes: List[GeneRegion] = []
    with open(bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            t = re.split(r"\s+", line.strip())
            if len(t) < 3:
                continue
            genes.append(GeneRegion(t[0], int(t[1]), int(t[2]), t[3] if len(t) >= 4 else f"{t[0]}:{int(t[1])+1}-{t[2]}", t[5] if len(t) >= 6 else "."))
    if not genes:
        raise PipelineError(f"No gene regions parsed from BED: {bed}")
    return genes


def validate_bed(genes: List[GeneRegion], contigs: Dict[str, int]) -> None:
    for g in genes:
        if g.chrom not in contigs:
            raise PipelineError(f"BED contig not found in reference: {g.chrom}")
        if g.start0 < 0 or g.end > contigs[g.chrom] or g.start0 >= g.end:
            raise PipelineError(f"Invalid BED interval: {g.chrom}:{g.start1}-{g.end}")


def normalize_name(x: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", x.replace("|", "_").replace(" ", "_"))


def load_metadata(metadata: Optional[Path]) -> Dict[Tuple[str, str], str]:
    if not metadata or not metadata.exists():
        return {}
    lines = metadata.read_text(errors="replace").splitlines()
    if not lines:
        return {}
    header = [x.strip().lower() for x in lines[0].split("\t")]

    def idx(names: List[str]) -> Optional[int]:
        for n in names:
            if n in header:
                return header.index(n)
        return None

    i_sample = idx(["sample", "sample_id", "sample_name", "id"])
    i_state = idx(["state", "site", "location"])
    i_barcode = idx(["barcode", "ont_barcode", "barcodes"])
    if i_sample is None or i_state is None or i_barcode is None:
        return {}

    out: Dict[Tuple[str, str], str] = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        t = line.split("\t")
        if len(t) <= max(i_sample, i_state, i_barcode):
            continue
        out[(t[i_state].strip(), t[i_barcode].strip())] = t[i_sample].strip()
    return out


def discover_barcode_samples(cfg: Config) -> List[BarcodeSample]:
    meta = load_metadata(cfg.metadata)
    samples: List[BarcodeSample] = []
    for state_dir in sorted([p for p in cfg.demux_dir.iterdir() if p.is_dir()]):
        fastq_pass = state_dir / "fastq_pass"
        if not fastq_pass.exists():
            continue
        for bc_dir in sorted([p for p in fastq_pass.glob("barcode*") if p.is_dir()]):
            state = state_dir.name
            barcode = bc_dir.name
            sample = meta.get((state, barcode), barcode)
            samples.append(BarcodeSample(state, barcode, bc_dir, sample))
    if not samples:
        raise PipelineError(f"No barcode directories found under {cfg.demux_dir}/*/fastq_pass/barcode*")
    return samples


def list_read_files(barcode_dir: Path) -> Tuple[List[Path], str]:
    pats = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz", "*.fasta", "*.fasta.gz", "*.fa", "*.fa.gz"]
    files: List[Path] = []
    for pat in pats:
        files.extend(sorted(barcode_dir.glob(pat)))
    files = sorted(set(files))
    if not files:
        return [], "none"
    if any(str(f).endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")) for f in files):
        return [f for f in files if str(f).endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))], "fastq"
    return [f for f in files if str(f).endswith((".fasta", ".fa", ".fasta.gz", ".fa.gz"))], "fasta"


def concat_reads(inputs: List[Path], out_gz: Path, logger: logging.Logger) -> None:
    ensure_parent(out_gz)
    with gzip.open(out_gz, "wb") as out:
        for fp in inputs:
            logger.info("Concatenating %s", fp)
            if str(fp).endswith(".gz"):
                with gzip.open(fp, "rb") as inp:
                    shutil.copyfileobj(inp, out)
            else:
                with open(fp, "rb") as inp:
                    shutil.copyfileobj(inp, out)
    if not nonempty(out_gz):
        raise PipelineError(f"Concatenated read file is empty: {out_gz}")


def fastqc(cfg: Config, files: List[Path], out_dir: Path, logger: logging.Logger) -> None:
    ensure_dir(out_dir)
    good = [f for f in files if nonempty(f)]
    if good:
        run_cmd([cfg.fastqc_bin, "-t", str(max(1, min(4, cfg.threads))), "-o", str(out_dir)] + [str(x) for x in good], logger)


def multiqc(cfg: Config, search_dir: Path, out_dir: Path, name: str, logger: logging.Logger) -> None:
    ensure_dir(out_dir)
    run_cmd([cfg.multiqc_bin, str(search_dir), "-o", str(out_dir), "-n", name, "--force"], logger)


def filter_reads(cfg: Config, input_gz: Path, out_fastq_gz: Path, logger: logging.Logger) -> Path:
    ensure_parent(out_fastq_gz)
    script = (
        f"set -euo pipefail; "
        f"{cfg.chopper_bin} --input {q(input_gz)} "
        f"--quality {cfg.min_qscore} --minlength {cfg.min_length} "
        f"--threads {max(1, cfg.threads)} | {cfg.bgzip_bin} -c > {q(out_fastq_gz)}"
    )
    run_bash(script, logger)
    return out_fastq_gz if nonempty(out_fastq_gz) else input_gz


def map_ont(cfg: Config, reads: Path, bam: Path, logger: logging.Logger) -> None:
    ensure_parent(bam)
    script = (
        f"set -euo pipefail; "
        f"{cfg.minimap2_bin} -t {max(1, cfg.threads)} -ax map-ont {q(cfg.reference)} {q(reads)} "
        f"| {cfg.samtools_bin} sort -@ {max(1, cfg.threads // 2)} -o {q(bam)} -"
    )
    run_bash(script, logger)
    run_cmd([cfg.samtools_bin, "index", str(bam)], logger)


def bam_stats(cfg: Config, bam: Path, out_dir: Path, prefix: str, logger: logging.Logger) -> None:
    ensure_dir(out_dir)
    run_cmd([cfg.samtools_bin, "flagstat", str(bam)], logger, stdout_path=out_dir / f"{prefix}.flagstat.txt")
    run_cmd([cfg.samtools_bin, "idxstats", str(bam)], logger, stdout_path=out_dir / f"{prefix}.idxstats.txt")
    run_cmd([cfg.samtools_bin, "stats", str(bam)], logger, stdout_path=out_dir / f"{prefix}.stats.txt")


def tabix_vcf(cfg: Config, vcf: Path, logger: logging.Logger) -> None:
    if not nonempty(vcf):
        raise PipelineError(f"Cannot index empty VCF: {vcf}")
    safe_unlink(Path(str(vcf) + ".tbi"))
    safe_unlink(Path(str(vcf) + ".csi"))
    run_cmd([cfg.tabix_bin, "-f", "-p", "vcf", str(vcf)], logger)


def call_clair3(cfg: Config, bam: Path, out_vcf: Path, sample: str, work_dir: Path, logger: logging.Logger) -> Optional[Path]:
    if not cfg.run_clair3:
        return None
    if cfg.resume and nonempty(out_vcf):
        tabix_vcf(cfg, out_vcf, logger)
        return out_vcf
    if not cfg.clair3_model_path:
        raise PipelineError("--clair3_model_path is required when --run_clair3 is enabled")
    ensure_dir(work_dir)
    call_dir = work_dir / "clair3"
    if call_dir.exists():
        shutil.rmtree(call_dir)
    ensure_dir(call_dir)
    cmd = [
        cfg.clair3_bin,
        f"--bam_fn={bam}",
        f"--ref_fn={cfg.reference}",
        f"--output={call_dir}",
        f"--threads={max(1, cfg.threads)}",
        f"--platform={cfg.clair3_platform}",
        f"--model_path={cfg.clair3_model_path}",
        f"--sample_name={sample}",
    ]
    if cfg.clair3_include_all_ctgs:
        cmd.append("--include_all_ctgs")
    run_cmd(cmd, logger)
    produced = call_dir / "merge_output.vcf.gz"
    if not produced.exists():
        hits = list(call_dir.rglob("*.vcf.gz"))
        if not hits:
            raise PipelineError(f"Clair3 produced no VCF in {call_dir}")
        produced = hits[0]
    shutil.copy2(produced, out_vcf)
    tabix_vcf(cfg, out_vcf, logger)
    return out_vcf


def call_freebayes(cfg: Config, bam: Path, out_vcf: Path, logger: logging.Logger) -> Optional[Path]:
    if not cfg.run_freebayes:
        return None
    ensure_parent(out_vcf)
    script = (
        f"set -euo pipefail; "
        f"{cfg.freebayes_bin} -f {q(cfg.reference)} --ploidy {cfg.freebayes_ploidy} "
        f"--min-base-quality {cfg.min_baseq} --min-mapping-quality {cfg.min_mapq} "
        f"--min-alternate-count 3 --min-alternate-fraction 0.02 --use-best-n-alleles 4 {q(bam)} "
        f"| {cfg.bcftools_bin} view -Oz -o {q(out_vcf)}"
    )
    run_bash(script, logger)
    tabix_vcf(cfg, out_vcf, logger)
    return out_vcf


def call_bcftools(cfg: Config, bam: Path, out_vcf: Path, logger: logging.Logger) -> Optional[Path]:
    if not cfg.run_bcftools:
        return None
    ensure_parent(out_vcf)
    script = (
        f"set -euo pipefail; "
        f"{cfg.bcftools_bin} mpileup -Ou -f {q(cfg.reference)} "
        f"-q {cfg.min_mapq} -Q {cfg.min_baseq} {q(bam)} "
        f"| {cfg.bcftools_bin} call -mv --ploidy {cfg.bcftools_ploidy} -Ou "
        f"| {cfg.bcftools_bin} view -Oz -o {q(out_vcf)}"
    )
    run_bash(script, logger)
    tabix_vcf(cfg, out_vcf, logger)
    return out_vcf


def sample_name_from_vcf(cfg: Config, vcf: Path, logger: logging.Logger) -> str:
    p = subprocess.run([cfg.bcftools_bin, "query", "-l", str(vcf)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        logger.info("bcftools query -l stderr: %s", p.stderr)
        raise PipelineError(f"Cannot read sample name from {vcf}")
    names = [x for x in p.stdout.splitlines() if x.strip()]
    return names[0] if names else "SAMPLE"


def reheader_vcf(cfg: Config, in_vcf: Path, out_vcf: Path, sample: str, logger: logging.Logger) -> None:
    ensure_parent(out_vcf)
    p = subprocess.run(
        [cfg.bcftools_bin, "reheader", "-s", "/dev/stdin", "-o", str(out_vcf), str(in_vcf)],
        input=sample + "\n",
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if p.returncode != 0:
        logger.info("reheader stderr: %s", p.stderr)
        raise PipelineError(f"bcftools reheader failed: {in_vcf}")
    tabix_vcf(cfg, out_vcf, logger)


def sanitize_vcf(cfg: Config, in_vcf: Path, out_vcf: Path, logger: logging.Logger) -> None:
    ensure_parent(out_vcf)
    script = (
        f"set -euo pipefail; "
        f"{cfg.bcftools_bin} annotate -x FORMAT/GQ,INFO/MQ,FORMAT/MQ {q(in_vcf)} -Ou "
        f"| {cfg.bcftools_bin} view -Oz -o {q(out_vcf)}"
    )
    run_bash(script, logger)
    tabix_vcf(cfg, out_vcf, logger)


def merge_callers(cfg: Config, vcfs: Dict[str, Path], out_vcf: Path, work_dir: Path, logger: logging.Logger) -> Path:
    if not vcfs:
        raise PipelineError("No caller VCFs available")
    ensure_dir(work_dir)
    canonical = sample_name_from_vcf(cfg, next(iter(vcfs.values())), logger)
    prepared: List[Path] = []
    for caller, vcf in vcfs.items():
        rh = work_dir / f"{caller}.rehead.vcf.gz"
        san = work_dir / f"{caller}.sanitized.vcf.gz"
        reheader_vcf(cfg, vcf, rh, canonical, logger)
        sanitize_vcf(cfg, rh, san, logger)
        prepared.append(san)
    ordered = [work_dir / f"{x}.sanitized.vcf.gz" for x in ["clair3", "freebayes", "bcftools"] if (work_dir / f"{x}.sanitized.vcf.gz").exists()]
    tmp = work_dir / "ensemble.concat.vcf.gz"
    script = (
        f"set -euo pipefail; "
        f"{cfg.bcftools_bin} concat -a -D -Oz -o {q(tmp)} "
        + " ".join(q(x) for x in ordered)
        + f"; {cfg.bcftools_bin} sort -Oz -o {q(out_vcf)} {q(tmp)}"
    )
    run_bash(script, logger)
    tabix_vcf(cfg, out_vcf, logger)
    return out_vcf


def filter_good_bad(cfg: Config, in_vcf: Path, good_vcf: Path, bad_vcf: Path, logger: logging.Logger) -> None:
    expr = f"QUAL>={cfg.min_site_qual} && (INFO/DP>={cfg.min_site_dp} || FORMAT/DP[0]>={cfg.min_site_dp})"
    script_good = f"set -euo pipefail; {cfg.bcftools_bin} filter -i '{expr}' {q(in_vcf)} -Oz -o {q(good_vcf)}"
    script_bad = f"set -euo pipefail; {cfg.bcftools_bin} filter -e '{expr}' {q(in_vcf)} -Oz -o {q(bad_vcf)}"
    run_bash(script_good, logger)
    run_bash(script_bad, logger)
    tabix_vcf(cfg, good_vcf, logger)
    tabix_vcf(cfg, bad_vcf, logger)


def snpeff(cfg: Config, in_vcf: Path, out_vcf: Path, logger: logging.Logger) -> None:
    ensure_parent(out_vcf)
    cmd = [str(cfg.snpeff_bin)]
    if cfg.snpeff_config:
        cmd += ["-c", str(cfg.snpeff_config)]
    cmd += [cfg.snpeff_db, str(in_vcf)]
    script = "set -euo pipefail; " + " ".join(q(x) for x in cmd) + f" | {cfg.bgzip_bin} -c > {q(out_vcf)}"
    run_bash(script, logger)
    tabix_vcf(cfg, out_vcf, logger)


def depth_bed(cfg: Config, bam: Path, bed: Path, out_tsv: Path, logger: logging.Logger) -> None:
    script = f"set -euo pipefail; {cfg.samtools_bin} depth -aa -b {q(bed)} {q(bam)} > {q(out_tsv)}"
    run_bash(script, logger)


def compute_gene_coverage(depth_tsv: Path, genes: List[GeneRegion], out_tsv: Path) -> List[Dict[str, object]]:
    per_gene = {g.name: [] for g in genes}
    idx = {(g.chrom, p): g.name for g in genes for p in range(g.start1, g.end + 1)}
    with open(depth_tsv) as fh:
        for line in fh:
            if not line.strip():
                continue
            c, p, d = line.rstrip("\n").split("\t")
            g = idx.get((c, int(p)))
            if g:
                per_gene[g].append(int(d))
    rows: List[Dict[str, object]] = []
    for g in genes:
        vals = per_gene[g.name]
        if len(vals) < g.length:
            vals += [0] * (g.length - len(vals))
        mean_dp = sum(vals) / g.length if g.length else 0
        rows.append({
            "gene": g.name,
            "chrom": g.chrom,
            "start_1based": g.start1,
            "end_1based": g.end,
            "length_bp": g.length,
            "mean_depth": round(mean_dp, 4),
            "median_depth": round(float(statistics.median(vals)) if vals else 0, 4),
            "max_depth": max(vals) if vals else 0,
            "breadth_1x_pct": round(sum(x >= 1 for x in vals) / g.length * 100 if g.length else 0, 4),
            "breadth_5x_pct": round(sum(x >= 5 for x in vals) / g.length * 100 if g.length else 0, 4),
            "breadth_10x_pct": round(sum(x >= 10 for x in vals) / g.length * 100 if g.length else 0, 4),
        })
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return rows


def gene_calls(cfg: Config, sample: str, rows: List[Dict[str, object]], out_tsv: Path) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    for r in rows:
        present = (
            float(r["mean_depth"]) >= cfg.gene_presence_min_mean_depth
            and float(r["breadth_1x_pct"]) >= cfg.gene_presence_min_breadth_1x
            and float(r["breadth_5x_pct"]) >= cfg.gene_presence_min_breadth_5x
            and float(r["breadth_10x_pct"]) >= cfg.gene_presence_min_breadth_10x
        )
        row = {"sample": sample, **r, "gene_status": "present" if present else "deleted_or_not_detected"}
        out.append(row)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    return out


def region_snps(vcf: Path, gene: GeneRegion) -> List[Tuple[str, int, str, str, str]]:
    snps = []
    with pysam.VariantFile(str(vcf)) as vf:
        for rec in vf.fetch(gene.chrom, gene.start0, gene.end):
            if rec.alts is None or len(rec.alts) != 1:
                continue
            ref = str(rec.ref).upper()
            alt = str(rec.alts[0]).upper()
            if len(ref) != 1 or len(alt) != 1:
                continue
            ann = ""
            try:
                x = rec.info.get("ANN", None)
                ann = ",".join(map(str, x)) if isinstance(x, (tuple, list)) else (str(x) if x else "")
            except Exception:
                pass
            snps.append((str(rec.contig), int(rec.pos), ref, alt, ann))
    return snps


def extract_gene_haplotypes(cfg: Config, sample: str, bam: Path, gene: GeneRegion, snps: List[Tuple[str, int, str, str, str]], out_prefix: Path) -> List[Dict[str, object]]:
    reads_tsv = Path(str(out_prefix) + ".read_haplotypes.tsv")
    haps_tsv = Path(str(out_prefix) + ".haplotypes.tsv")
    fasta = Path(str(out_prefix) + ".haplotypes.fasta")
    ensure_parent(reads_tsv)
    if not snps:
        reads_tsv.write_text("sample\tgene\tread_name\thaplotype_string\tinformative_sites\n")
        haps_tsv.write_text("sample\tgene\thaplotype_id\thaplotype_string\tsupporting_reads\tfrequency\tinformative_sites\n")
        fasta.write_text("")
        return []
    snp_index = {(c, p): i for i, (c, p, _r, _a, _ann) in enumerate(snps)}
    refalt = {(c, p): (r, a) for c, p, r, a, _ann in snps}
    read_calls: Dict[str, List[str]] = {}
    with pysam.AlignmentFile(str(bam), "rb") as bf:
        for read in bf.fetch(gene.chrom, gene.start0, gene.end):
            if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < cfg.min_mapq or read.query_sequence is None:
                continue
            read_calls.setdefault(read.query_name, ["N"] * len(snps))
            seq = read.query_sequence.upper()
            quals = read.query_qualities or []
            for qpos, rpos in read.get_aligned_pairs(matches_only=False):
                if qpos is None or rpos is None:
                    continue
                key = (gene.chrom, rpos + 1)
                i = snp_index.get(key)
                if i is None or qpos >= len(seq) or qpos >= len(quals):
                    continue
                if int(quals[qpos]) < cfg.min_baseq:
                    continue
                base = seq[qpos]
                ref, alt = refalt[key]
                if base == ref or base == alt:
                    read_calls[read.query_name][i] = base
    filtered = []
    with open(reads_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "gene", "read_name", "haplotype_string", "informative_sites"])
        for rn, h in read_calls.items():
            hs = "".join(h)
            ninfo = sum(x != "N" for x in h)
            if ninfo >= cfg.hap_min_informative_sites:
                filtered.append(hs)
                w.writerow([sample, gene.name, rn, hs, ninfo])
    counts = Counter(filtered)
    total = sum(counts.values())
    min_reads = max(cfg.hap_min_reads, int(total * cfg.hap_min_fraction))
    rows = []
    with open(haps_tsv, "w", newline="") as fh, open(fasta, "w") as fa:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "gene", "haplotype_id", "haplotype_string", "supporting_reads", "frequency", "informative_sites"])
        n = 0
        for hap, cnt in counts.most_common():
            if cnt < min_reads:
                continue
            n += 1
            hid = f"{normalize_name(gene.name)}_HAP{n}"
            freq = cnt / total if total else 0
            row = {
                "sample": sample,
                "gene": gene.name,
                "haplotype_id": hid,
                "haplotype_string": hap,
                "supporting_reads": cnt,
                "frequency": round(freq, 6),
                "informative_sites": sum(x != "N" for x in hap),
            }
            rows.append(row)
            w.writerow(row.values())
            fa.write(f">{sample}|{gene.name}|{hid}|reads={cnt}|freq={freq:.6f}\n{hap}\n")
    return rows


def write_region_variant_table(sample: str, vcf: Path, genes: List[GeneRegion], out_tsv: Path) -> None:
    rows = []
    with pysam.VariantFile(str(vcf)) as vf:
        for g in genes:
            for rec in vf.fetch(g.chrom, g.start0, g.end):
                ann = ""
                try:
                    x = rec.info.get("ANN", None)
                    ann = ",".join(map(str, x)) if isinstance(x, (tuple, list)) else (str(x) if x else "")
                except Exception:
                    pass
                if rec.alts:
                    for alt in rec.alts:
                        rows.append({"sample": sample, "gene": g.name, "chrom": rec.contig, "pos": int(rec.pos), "ref": rec.ref, "alt": str(alt), "qual": rec.qual, "annotation": ann})
    ensure_parent(out_tsv)
    fields = ["sample", "gene", "chrom", "pos", "ref", "alt", "qual", "annotation"]
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def run_deploid(cfg: Config, sample: str, vcf: Path, out_dir: Path, logger: logging.Logger) -> Optional[Path]:
    if not cfg.run_deploid:
        return None
    ensure_dir(out_dir)
    prefix = out_dir / sample
    cmd = [cfg.deploid_bin, "-vcf", str(vcf), "-o", str(prefix)]
    if cfg.population_vcf:
        cmd += ["-panel", str(cfg.population_vcf)]
    if cfg.genetic_map:
        cmd += ["-map", str(cfg.genetic_map)]
    if cfg.exclude_vcf:
        cmd += ["-exclude", str(cfg.exclude_vcf)]
    run_cmd(cmd, logger)
    return prefix


def write_json(obj: object, path: Path) -> None:
    ensure_parent(path)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2)


def process_sample(cfg: Config, bs: BarcodeSample, genes: List[GeneRegion]) -> Optional[Dict[str, str]]:
    logger = get_logger(f"{bs.state}.{bs.barcode}")
    sample_id = f"{bs.state}_{bs.barcode}"
    sample = bs.sample or sample_id
    root = cfg.output_dir / bs.state / bs.barcode
    dirs = {
        "concat": root / "concat",
        "qc_raw": root / "qc_raw",
        "filtered": root / "filtered_fastq",
        "qc_filtered": root / "qc_after_filtering",
        "bam": root / "bam",
        "stats": root / "stats",
        "variants": root / "variants",
        "coverage": root / "coverage",
        "haplotypes": root / "haplotypes",
        "deploid": root / "deploid",
        "multiqc": root / "multiqc",
        "summary": root / "summary",
        "tmp": root / "tmp",
    }
    for d in dirs.values():
        ensure_dir(d)
    try:
        print(step(f"{sample_id} :: concatenate Nanopore barcode reads"))
        read_files, mode = list_read_files(bs.barcode_dir)
        if not read_files:
            raise PipelineError(f"No FASTQ/FASTA found in {bs.barcode_dir}")
        concat = dirs["concat"] / f"{bs.barcode}.{mode}.gz"
        if not (cfg.resume and nonempty(concat)):
            concat_reads(read_files, concat, logger)

        print(step(f"{sample_id} :: QC raw reads"))
        fastqc(cfg, [concat], dirs["qc_raw"], logger)

        print(step(f"{sample_id} :: length and quality filtering"))
        filtered = dirs["filtered"] / f"{bs.barcode}.filtered.fastq.gz"
        reads_for_map = filter_reads(cfg, concat, filtered, logger) if mode == "fastq" else concat

        print(step(f"{sample_id} :: QC after filtering"))
        fastqc(cfg, [reads_for_map], dirs["qc_filtered"], logger)

        print(step(f"{sample_id} :: minimap2 mapping"))
        bam = dirs["bam"] / f"{bs.barcode}.Pf3D7.bam"
        if not (cfg.resume and nonempty(bam) and nonempty(Path(str(bam) + ".bai"))):
            map_ont(cfg, reads_for_map, bam, logger)
        bam_stats(cfg, bam, dirs["stats"], bs.barcode, logger)

        print(step(f"{sample_id} :: variant calling with Clair3, FreeBayes and bcftools"))
        vcfs: Dict[str, Path] = {}
        clair3_vcf = dirs["variants"] / f"{bs.barcode}.clair3.raw.vcf.gz"
        fb_vcf = dirs["variants"] / f"{bs.barcode}.freebayes.raw.vcf.gz"
        bc_vcf = dirs["variants"] / f"{bs.barcode}.bcftools.raw.vcf.gz"

        x = call_clair3(cfg, bam, clair3_vcf, sample, dirs["tmp"], logger)
        if x:
            vcfs["clair3"] = x
        x = call_freebayes(cfg, bam, fb_vcf, logger)
        if x:
            vcfs["freebayes"] = x
        x = call_bcftools(cfg, bam, bc_vcf, logger)
        if x:
            vcfs["bcftools"] = x

        ensemble = dirs["variants"] / f"{bs.barcode}.ensemble.vcf.gz"
        annotated = dirs["variants"] / f"{bs.barcode}.ensemble.ann.vcf.gz"
        good_vcf = dirs["variants"] / f"{bs.barcode}.ensemble.good.vcf.gz"
        bad_vcf = dirs["variants"] / f"{bs.barcode}.ensemble.bad.vcf.gz"

        merge_callers(cfg, vcfs, ensemble, dirs["tmp"] / "ensemble", logger)
        snpeff(cfg, ensemble, annotated, logger)
        filter_good_bad(cfg, annotated, good_vcf, bad_vcf, logger)

        print(step(f"{sample_id} :: HRP2/HRP3 deletion by gene coverage"))
        depth_tsv = dirs["coverage"] / f"{bs.barcode}.hrp2_hrp3.depth.tsv"
        cov_tsv = dirs["coverage"] / f"{bs.barcode}.hrp2_hrp3.coverage.tsv"
        calls_tsv = dirs["coverage"] / f"{bs.barcode}.hrp2_hrp3.gene_calls.tsv"
        depth_bed(cfg, bam, cfg.genes_bed, depth_tsv, logger)
        cov_rows = compute_gene_coverage(depth_tsv, genes, cov_tsv)
        call_rows = gene_calls(cfg, sample, cov_rows, calls_tsv)

        print(step(f"{sample_id} :: HRP2/HRP3 variants and clonal haplotypes"))
        variant_tsv = dirs["variants"] / f"{bs.barcode}.hrp2_hrp3.variants.tsv"
        write_region_variant_table(sample, good_vcf, genes, variant_tsv)

        all_haps: List[Dict[str, object]] = []
        hap_summary = []
        for g in genes:
            snps = region_snps(good_vcf, g)
            rows = extract_gene_haplotypes(cfg, sample, bam, g, snps, dirs["haplotypes"] / f"{bs.barcode}.{normalize_name(g.name)}")
            all_haps.extend(rows)
            status = next((x["gene_status"] for x in call_rows if x["gene"] == g.name), "deleted_or_not_detected")
            hap_summary.append({
                "sample": sample,
                "gene": g.name,
                "gene_status": status,
                "n_gene_snps": len(snps),
                "n_clonal_haplotypes": len(rows),
                "parasite_level_interpretation": "multi-clonal" if len(rows) > 1 else ("single-clone" if len(rows) == 1 else "not_resolved"),
            })

        haps_tsv = dirs["haplotypes"] / f"{bs.barcode}.all_genes.clonal_haplotypes.tsv"
        with open(haps_tsv, "w", newline="") as fh:
            fields = ["sample", "gene", "haplotype_id", "haplotype_string", "supporting_reads", "frequency", "informative_sites"]
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(all_haps)

        print(step(f"{sample_id} :: optional DEploid2"))
        deploid_prefix = run_deploid(cfg, sample, good_vcf, dirs["deploid"], logger)

        print(step(f"{sample_id} :: MultiQC"))
        multiqc(cfg, root, dirs["multiqc"], f"{bs.barcode}.multiqc_report.html", logger)

        summary_json = dirs["summary"] / f"{bs.barcode}.summary.json"
        write_json({
            "author": AUTHOR,
            "sample": sample,
            "state": bs.state,
            "barcode": bs.barcode,
            "n_input_fastq_files": len(read_files),
            "concat_reads": str(concat),
            "filtered_reads": str(reads_for_map),
            "bam": str(bam),
            "ensemble_vcf": str(ensemble),
            "annotated_vcf": str(annotated),
            "good_vcf": str(good_vcf),
            "bad_vcf": str(bad_vcf),
            "gene_calls": call_rows,
            "clonal_haplotype_summary": hap_summary,
            "deploid_prefix": str(deploid_prefix) if deploid_prefix else None,
            "multiqc_report": str(dirs["multiqc"] / f"{bs.barcode}.multiqc_report.html"),
        }, summary_json)

        if not cfg.keep_tmp:
            shutil.rmtree(dirs["tmp"], ignore_errors=True)

        print(done(f"{sample_id} :: complete"))
        return {
            "sample": sample,
            "state": bs.state,
            "barcode": bs.barcode,
            "gene_calls": str(calls_tsv),
            "coverage": str(cov_tsv),
            "haplotypes": str(haps_tsv),
            "variants": str(variant_tsv),
            "summary": str(summary_json),
        }

    except Exception as e:
        logger.error("Sample failed: %s", e)
        logger.error("TRACEBACK:\n%s", traceback.format_exc())
        print(fail(f"{sample_id} :: {e}"))
        return None


def aggregate_tsv(files: List[Path], out_tsv: Path, empty_header: List[str]) -> None:
    rows = []
    fields = None
    for fp in files:
        if not fp.exists():
            continue
        with open(fp) as fh:
            dr = csv.DictReader(fh, delimiter="\t")
            if fields is None:
                fields = dr.fieldnames
            rows.extend(list(dr))
    ensure_parent(out_tsv)
    with open(out_tsv, "w", newline="") as fh:
        if rows and fields:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        else:
            csv.writer(fh, delimiter="\t").writerow(empty_header)


def aggregate_prevalence(gene_call_files: List[Path], out_tsv: Path) -> None:
    rows = []
    for fp in gene_call_files:
        with open(fp) as fh:
            rows.extend(list(csv.DictReader(fh, delimiter="\t")))
    by_gene = defaultdict(list)
    for r in rows:
        by_gene[r["gene"]].append(r)
    out = []
    for gene, vals in sorted(by_gene.items()):
        n = len(vals)
        present = sum(v["gene_status"] == "present" for v in vals)
        deleted = n - present
        out.append({
            "gene": gene,
            "n_samples": n,
            "n_present": present,
            "n_deleted_or_not_detected": deleted,
            "prevalence_deleted_or_not_detected_pct": round(deleted / n * 100 if n else 0, 4),
        })
    ensure_parent(out_tsv)
    with open(out_tsv, "w", newline="") as fh:
        fields = ["gene", "n_samples", "n_present", "n_deleted_or_not_detected", "prevalence_deleted_or_not_detected_pct"]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out)


def build_layout(out: Path) -> None:
    for p in [out / "cohort" / "coverage", out / "cohort" / "variants", out / "cohort" / "haplotypes", out / "cohort" / "multiqc", out / "logs"]:
        ensure_dir(p)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Nanopore PF HRP2/HRP3 deletion pipeline with three callers and optional DEploid2.")
    p.add_argument("--demux_dir", required=True, type=Path)
    p.add_argument("--output_dir", required=True, type=Path)
    p.add_argument("--reference", required=True, type=Path)
    p.add_argument("--genes_bed", required=True, type=Path)
    p.add_argument("--metadata", type=Path, default=None)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--sample_workers", type=int, default=2)
    p.add_argument("--min_length", type=int, default=150)
    p.add_argument("--min_qscore", type=int, default=10)
    p.add_argument("--min_mapq", type=int, default=20)
    p.add_argument("--min_baseq", type=int, default=15)
    p.add_argument("--min_site_dp", type=int, default=20)
    p.add_argument("--min_site_qual", type=int, default=20)
    p.add_argument("--gene_presence_min_mean_depth", type=float, default=5.0)
    p.add_argument("--gene_presence_min_breadth_1x", type=float, default=90.0)
    p.add_argument("--gene_presence_min_breadth_5x", type=float, default=80.0)
    p.add_argument("--gene_presence_min_breadth_10x", type=float, default=70.0)
    p.add_argument("--fastqc_bin", default="fastqc")
    p.add_argument("--multiqc_bin", default="multiqc")
    p.add_argument("--chopper_bin", default="chopper")
    p.add_argument("--minimap2_bin", default="minimap2")
    p.add_argument("--samtools_bin", default="samtools")
    p.add_argument("--bcftools_bin", default="bcftools")
    p.add_argument("--freebayes_bin", default="freebayes")
    p.add_argument("--bgzip_bin", default="bgzip")
    p.add_argument("--tabix_bin", default="tabix")
    p.add_argument("--clair3_bin", default="run_clair3.sh")
    p.add_argument("--clair3_model_path", type=Path, default=None)
    p.add_argument("--clair3_platform", choices=["ont", "hifi", "ilmn"], default="ont")
    p.add_argument("--clair3_include_all_ctgs", action="store_true")
    p.add_argument("--run_clair3", type=int, default=1)
    p.add_argument("--run_freebayes", type=int, default=1)
    p.add_argument("--run_bcftools", type=int, default=1)
    p.add_argument("--freebayes_ploidy", type=int, default=2)
    p.add_argument("--bcftools_ploidy", type=int, default=2)
    p.add_argument("--snpeff", required=True, type=Path)
    p.add_argument("--snpeff_db", default="Pf3D7_v2")
    p.add_argument("--snpeff_config", type=Path, default=None)
    p.add_argument("--run_deploid", action="store_true")
    p.add_argument("--deploid_bin", default="dEploid")
    p.add_argument("--population_vcf", type=Path, default=None)
    p.add_argument("--genetic_map", type=Path, default=None)
    p.add_argument("--exclude_vcf", type=Path, default=None)
    p.add_argument("--hap_min_reads", type=int, default=3)
    p.add_argument("--hap_min_informative_sites", type=int, default=2)
    p.add_argument("--hap_min_fraction", type=float, default=0.005)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--keep_tmp", action="store_true")
    p.add_argument("--log_file", type=Path, default=Path("pf_hrp_ont.log"))
    return p.parse_args()


def banner(cfg: Config, n: int) -> None:
    print(f"{MAGENTA}{BOLD}======================================================================{RESET}")
    print(f"{MAGENTA}{BOLD} PF-HRP2/HRP3 Nanopore Pipeline{RESET}")
    print(f"{MAGENTA}{BOLD} Author: {AUTHOR}{RESET}")
    print(f"{MAGENTA}{BOLD}======================================================================{RESET}")
    print(info(f"Demux directory        : {cfg.demux_dir}"))
    print(info(f"Output directory       : {cfg.output_dir}"))
    print(info(f"Reference              : {cfg.reference}"))
    print(info(f"Genes BED              : {cfg.genes_bed}"))
    print(info(f"Barcode samples        : {n}"))
    print(info("Callers                : Clair3, FreeBayes, bcftools"))
    print(info(f"DEploid2               : {'enabled' if cfg.run_deploid else 'disabled'}"))
    print(info("Goal                   : HRP2/HRP3 deletion at sample and clonal/parasite level"))
    print(f"{MAGENTA}{BOLD}======================================================================{RESET}")


def main() -> int:
    args = parse_args()
    cfg = Config(
        demux_dir=args.demux_dir,
        output_dir=args.output_dir,
        reference=args.reference,
        genes_bed=args.genes_bed,
        metadata=args.metadata,
        threads=args.threads,
        sample_workers=args.sample_workers,
        min_length=args.min_length,
        min_qscore=args.min_qscore,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        min_site_dp=args.min_site_dp,
        min_site_qual=args.min_site_qual,
        gene_presence_min_mean_depth=args.gene_presence_min_mean_depth,
        gene_presence_min_breadth_1x=args.gene_presence_min_breadth_1x,
        gene_presence_min_breadth_5x=args.gene_presence_min_breadth_5x,
        gene_presence_min_breadth_10x=args.gene_presence_min_breadth_10x,
        fastqc_bin=args.fastqc_bin,
        multiqc_bin=args.multiqc_bin,
        chopper_bin=args.chopper_bin,
        minimap2_bin=args.minimap2_bin,
        samtools_bin=args.samtools_bin,
        bcftools_bin=args.bcftools_bin,
        freebayes_bin=args.freebayes_bin,
        bgzip_bin=args.bgzip_bin,
        tabix_bin=args.tabix_bin,
        clair3_bin=args.clair3_bin,
        clair3_model_path=args.clair3_model_path,
        clair3_platform=args.clair3_platform,
        clair3_include_all_ctgs=args.clair3_include_all_ctgs,
        run_clair3=bool(args.run_clair3),
        run_freebayes=bool(args.run_freebayes),
        run_bcftools=bool(args.run_bcftools),
        freebayes_ploidy=args.freebayes_ploidy,
        bcftools_ploidy=args.bcftools_ploidy,
        snpeff_bin=args.snpeff,
        snpeff_db=args.snpeff_db,
        snpeff_config=args.snpeff_config,
        run_deploid=args.run_deploid,
        deploid_bin=args.deploid_bin,
        population_vcf=args.population_vcf,
        genetic_map=args.genetic_map,
        exclude_vcf=args.exclude_vcf,
        hap_min_reads=args.hap_min_reads,
        hap_min_informative_sites=args.hap_min_informative_sites,
        hap_min_fraction=args.hap_min_fraction,
        resume=args.resume,
        keep_tmp=args.keep_tmp,
        log_file=args.log_file,
    )

    build_layout(cfg.output_dir)
    setup_logger(cfg.output_dir / "logs" / cfg.log_file.name)
    logger = get_logger("main")

    try:
        check_dependencies(cfg)
        index_reference(cfg.reference, cfg, logger)
        genes = parse_bed(cfg.genes_bed)
        validate_bed(genes, fai_contigs(Path(str(cfg.reference) + ".fai")))
        samples = discover_barcode_samples(cfg)
        banner(cfg, len(samples))

        results: List[Dict[str, str]] = []
        workers = max(1, min(cfg.sample_workers, len(samples)))
        if workers == 1:
            for s in samples:
                r = process_sample(cfg, s, genes)
                if r:
                    results.append(r)
        else:
            with cf.ThreadPoolExecutor(max_workers=workers) as ex:
                futs = {ex.submit(process_sample, cfg, s, genes): s for s in samples}
                for fut in cf.as_completed(futs):
                    r = fut.result()
                    if r:
                        results.append(r)

        gene_files = [Path(x["gene_calls"]) for x in results]
        cov_files = [Path(x["coverage"]) for x in results]
        var_files = [Path(x["variants"]) for x in results]
        hap_files = [Path(x["haplotypes"]) for x in results]

        print(step("Cohort :: aggregate coverage and gene deletion calls"))
        aggregate_tsv(gene_files, cfg.output_dir / "cohort" / "coverage" / "all_samples.hrp2_hrp3.gene_calls.tsv", [])
        aggregate_tsv(cov_files, cfg.output_dir / "cohort" / "coverage" / "all_samples.hrp2_hrp3.coverage.tsv", [])
        aggregate_prevalence(gene_files, cfg.output_dir / "cohort" / "coverage" / "hrp2_hrp3.deletion_prevalence.tsv")

        print(step("Cohort :: aggregate variants and clonal haplotypes"))
        aggregate_tsv(var_files, cfg.output_dir / "cohort" / "variants" / "all_samples.hrp2_hrp3.variants.tsv", [])
        aggregate_tsv(hap_files, cfg.output_dir / "cohort" / "haplotypes" / "all_samples.clonal_haplotypes.tsv", [])

        print(step("Cohort :: final MultiQC"))
        multiqc(cfg, cfg.output_dir, cfg.output_dir / "cohort" / "multiqc", "cohort.multiqc_report.html", logger)

        write_json({
            "author": AUTHOR,
            "n_barcode_samples_detected": len(samples),
            "n_barcode_samples_completed": len(results),
            "reference": str(cfg.reference),
            "genes_bed": str(cfg.genes_bed),
            "callers": {
                "clair3": cfg.run_clair3,
                "freebayes": cfg.run_freebayes,
                "bcftools": cfg.run_bcftools,
            },
            "deploid2": cfg.run_deploid,
            "multiqc_report": str(cfg.output_dir / "cohort" / "multiqc" / "cohort.multiqc_report.html"),
        }, cfg.output_dir / "cohort" / "cohort_summary.json")

        print(done("Pipeline completed"))
        print(info(f"Log: {cfg.output_dir / 'logs' / cfg.log_file.name}"))
        return 0

    except Exception as e:
        logger.error("Pipeline failed: %s", e)
        logger.error("TRACEBACK:\n%s", traceback.format_exc())
        print(fail(str(e)))
        return 1


if __name__ == "__main__":
    sys.exit(main())