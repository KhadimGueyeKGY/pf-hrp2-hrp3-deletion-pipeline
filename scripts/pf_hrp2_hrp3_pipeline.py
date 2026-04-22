#!/usr/bin/env python3
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import gzip
import json
import logging
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import traceback
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pysam


RESET = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"
RED = "\033[31m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
BLUE = "\033[34m"
MAGENTA = "\033[35m"
CYAN = "\033[36m"


def c(tag: str, color: str, msg: str) -> str:
    return f"{color}{BOLD}{tag}{RESET} {msg}"


def info(msg: str) -> str:
    return c("[INFO]", CYAN, msg)


def step(msg: str) -> str:
    return c("[STEP]", BLUE, msg)


def done(msg: str) -> str:
    return c("[DONE]", GREEN, msg)


def warn(msg: str) -> str:
    return c("[WARN]", YELLOW, msg)


def err(msg: str) -> str:
    return c("[FAIL]", RED, msg)


class PipelineError(RuntimeError):
    pass


@dataclass(frozen=True)
class SamplePair:
    sample: str
    r1: Path
    r2: Path


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
    input_dir: Path
    output_dir: Path
    pf_ref: Path
    human_ref: Path
    genes_bed: Path
    threads: int
    sample_workers: int
    fastp_bin: str
    fastqc_bin: str
    bwa_bin: str
    samtools_bin: str
    bedtools_bin: str
    bcftools_bin: str
    bgzip_bin: str
    tabix_bin: str
    snpeff_jar: Path
    snpeff_db: str
    snpeff_config: Optional[Path]
    min_read_len: int
    min_base_qual: int
    min_mapq: int
    variant_min_dp: int
    variant_min_qual: int
    gene_presence_min_mean_depth: float
    gene_presence_min_breadth_1x: float
    gene_presence_min_breadth_5x: float
    gene_presence_min_breadth_10x: float
    hap_min_reads: int
    hap_min_informative_sites: int
    hap_min_mapq: int
    hap_min_baseq: int
    resume: bool
    keep_tmp: bool
    log_file: Path


def setup_logger(log_file: Path) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("pf_hrp_pipeline")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    fh = logging.FileHandler(log_file, mode="a")
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s %(message)s"))
    logger.addHandler(fh)
    logger.propagate = False
    return logger


def get_logger(name: str) -> logging.Logger:
    lg = logging.getLogger(f"pf_hrp_pipeline.{name}")
    lg.setLevel(logging.INFO)
    lg.propagate = True
    return lg


def run_cmd(cmd: List[str], logger: logging.Logger, cwd: Optional[Path] = None, stdout_path: Optional[Path] = None) -> None:
    logger.info("CMD: %s", " ".join(map(str, cmd)))
    if stdout_path is None:
        p = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        if p.stdout:
            logger.info("STDOUT:\n%s", p.stdout)
        if p.stderr:
            logger.info("STDERR:\n%s", p.stderr)
        if p.returncode != 0:
            raise PipelineError(f"Command failed: {' '.join(map(str, cmd))}")
    else:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        with open(stdout_path, "wb") as out_fh:
            p = subprocess.run(
                cmd,
                cwd=str(cwd) if cwd else None,
                stdout=out_fh,
                stderr=subprocess.PIPE,
                check=False,
            )
        if p.stderr:
            logger.info("STDERR:\n%s", p.stderr.decode(errors="replace"))
        if p.returncode != 0:
            raise PipelineError(f"Command failed: {' '.join(map(str, cmd))}")


def run_bash(script: str, logger: logging.Logger, cwd: Optional[Path] = None) -> None:
    logger.info("BASH: %s", script)
    p = subprocess.run(
        ["bash", "-lc", script],
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if p.stdout:
        logger.info("STDOUT:\n%s", p.stdout)
    if p.stderr:
        logger.info("STDERR:\n%s", p.stderr)
    if p.returncode != 0:
        raise PipelineError(f"Bash failed: {script}")


def which_or_raise(x: str) -> str:
    p = shutil.which(x)
    if not p:
        raise PipelineError(f"Missing dependency in PATH: {x}")
    return p


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def safe_unlink(p: Path) -> None:
    try:
        p.unlink(missing_ok=True)
    except Exception:
        pass


def file_nonempty(p: Path) -> bool:
    return p.exists() and p.stat().st_size > 0


def fai_contigs(fai: Path) -> Dict[str, int]:
    out: Dict[str, int] = {}
    with open(fai, "r") as fh:
        for line in fh:
            if not line.strip():
                continue
            t = line.rstrip("\n").split("\t")
            out[t[0]] = int(t[1])
    return out


def parse_bed(bed: Path) -> List[GeneRegion]:
    genes: List[GeneRegion] = []
    with open(bed, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            t = re.split(r"\s+", line.strip())
            if len(t) < 3:
                continue
            chrom = t[0]
            start0 = int(t[1])
            end = int(t[2])
            name = t[3] if len(t) >= 4 else f"{chrom}:{start0 + 1}-{end}"
            strand = t[5] if len(t) >= 6 else "."
            genes.append(GeneRegion(chrom, start0, end, name, strand))
    if not genes:
        raise PipelineError(f"No regions parsed from BED: {bed}")
    return genes


def normalize_gene_name(x: str) -> str:
    x2 = x.replace("|", "_").replace(" ", "_")
    x2 = re.sub(r"[^A-Za-z0-9._-]+", "_", x2)
    return x2


def discover_pairs(input_dir: Path) -> List[SamplePair]:
    all_files = []
    for ext in ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]:
        all_files.extend(input_dir.rglob(ext))
    all_files = sorted(set(all_files))
    if not all_files:
        raise PipelineError(f"No FASTQ files found under: {input_dir}")
    r1_map: Dict[str, Path] = {}
    r2_map: Dict[str, Path] = {}

    pats = [
        (re.compile(r"(.+?)(?:_R?1(?:_001)?)(\.f(?:ast)?q(?:\.gz)?)$", re.IGNORECASE), 1),
        (re.compile(r"(.+?)(?:_R?2(?:_001)?)(\.f(?:ast)?q(?:\.gz)?)$", re.IGNORECASE), 2),
        (re.compile(r"(.+?)(?:\.R1)(\.f(?:ast)?q(?:\.gz)?)$", re.IGNORECASE), 1),
        (re.compile(r"(.+?)(?:\.R2)(\.f(?:ast)?q(?:\.gz)?)$", re.IGNORECASE), 2),
        (re.compile(r"(.+?)(?:[_\.]1)(\.f(?:ast)?q(?:\.gz)?)$", re.IGNORECASE), 1),
        (re.compile(r"(.+?)(?:[_\.]2)(\.f(?:ast)?q(?:\.gz)?)$", re.IGNORECASE), 2),
    ]

    for fp in all_files:
        name = fp.name
        matched = False
        for rx, mate in pats:
            m = rx.match(name)
            if m:
                sample = m.group(1)
                if mate == 1:
                    r1_map[sample] = fp
                else:
                    r2_map[sample] = fp
                matched = True
                break
        if not matched:
            continue

    samples = sorted(set(r1_map) & set(r2_map))
    if not samples:
        raise PipelineError("No paired FASTQ files discovered.")
    pairs = [SamplePair(sample=s, r1=r1_map[s], r2=r2_map[s]) for s in samples]
    return pairs


def check_dependencies(cfg: Config) -> None:
    which_or_raise(cfg.fastp_bin)
    which_or_raise(cfg.fastqc_bin)
    which_or_raise(cfg.bwa_bin)
    which_or_raise(cfg.samtools_bin)
    which_or_raise(cfg.bedtools_bin)
    which_or_raise(cfg.bcftools_bin)
    which_or_raise(cfg.bgzip_bin)
    which_or_raise(cfg.tabix_bin)
    if not cfg.snpeff_jar.exists():
        raise PipelineError(f"snpEff jar not found: {cfg.snpeff_jar}")


def index_reference_bwa(ref: Path, bwa_bin: str, samtools_bin: str, logger: logging.Logger) -> None:
    expected = [Path(str(ref) + x) for x in [".amb", ".ann", ".bwt", ".pac", ".sa"]]
    if not all(p.exists() for p in expected):
        run_cmd([bwa_bin, "index", str(ref)], logger)
    fai = Path(str(ref) + ".fai")
    if not fai.exists():
        run_cmd([samtools_bin, "faidx", str(ref)], logger)


def ensure_samtools_fastq_extract_nonhuman(
    samtools_bin: str,
    bam_in: Path,
    r1_out: Path,
    r2_out: Path,
    singleton_out: Path,
    logger: logging.Logger,
) -> None:
    script = (
        f"set -euo pipefail; "
        f"{samtools_bin} fastq "
        f"-1 {shlex_q(r1_out)} "
        f"-2 {shlex_q(r2_out)} "
        f"-s {shlex_q(singleton_out)} "
        f"-0 /dev/null "
        f"-f 12 -F 256 "
        f"{shlex_q(bam_in)}"
    )
    run_bash(script, logger)


def shlex_q(p: Path | str) -> str:
    s = str(p)
    return "'" + s.replace("'", "'\"'\"'") + "'"


def fastqc_on_files(cfg: Config, files: List[Path], out_dir: Path, logger: logging.Logger) -> None:
    ensure_dir(out_dir)
    cmd = [cfg.fastqc_bin, "-t", str(max(1, min(cfg.threads, 4))), "-o", str(out_dir)] + [str(x) for x in files]
    run_cmd(cmd, logger)


def run_fastp(cfg: Config, r1: Path, r2: Path, out_r1: Path, out_r2: Path, html: Path, json_out: Path, logger: logging.Logger) -> None:
    ensure_dir(out_r1.parent)
    cmd = [
        cfg.fastp_bin,
        "-i", str(r1),
        "-I", str(r2),
        "-o", str(out_r1),
        "-O", str(out_r2),
        "--thread", str(max(1, cfg.threads)),
        "--detect_adapter_for_pe",
        "--length_required", str(cfg.min_read_len),
        "--json", str(json_out),
        "--html", str(html),
    ]
    run_cmd(cmd, logger)


def map_paired_bwa(cfg: Config, ref: Path, r1: Path, r2: Path, out_bam: Path, logger: logging.Logger) -> None:
    ensure_dir(out_bam.parent)

    bwa_threads = max(1, cfg.threads // 2)
    sort_threads = max(1, cfg.threads // 2)

    script = (
        f"{cfg.bwa_bin} mem -t {bwa_threads} {shlex_q(ref)} {shlex_q(r1)} {shlex_q(r2)} "
        f"| {cfg.samtools_bin} sort -@ {sort_threads} -o {shlex_q(out_bam)} -"
    )
    #print (script)
    #run_bash(script, logger)
    os.system(script)
    run_cmd([cfg.samtools_bin, "index", str(out_bam)], logger)


def bam_stats(cfg: Config, bam: Path, flagstat_txt: Path, idxstats_txt: Path, stats_txt: Path, logger: logging.Logger) -> None:
    run_cmd([cfg.samtools_bin, "flagstat", str(bam)], logger, stdout_path=flagstat_txt)
    run_cmd([cfg.samtools_bin, "idxstats", str(bam)], logger, stdout_path=idxstats_txt)
    run_cmd([cfg.samtools_bin, "stats", str(bam)], logger, stdout_path=stats_txt)


def depth_bed(cfg: Config, bam: Path, bed: Path, out_tsv: Path, logger: logging.Logger) -> None:
    script = (
        f"set -euo pipefail; "
        f"{cfg.samtools_bin} depth -aa -b {shlex_q(bed)} {shlex_q(bam)} > {shlex_q(out_tsv)}"
    )
    run_bash(script, logger)


def compute_gene_coverage(depth_tsv: Path, genes: List[GeneRegion], out_tsv: Path) -> List[Dict[str, object]]:
    per_gene_depths: Dict[str, List[int]] = {g.name: [] for g in genes}
    gene_by_key = {(g.chrom, pos): g.name for g in genes for pos in range(g.start1, g.end + 1)}
    with open(depth_tsv, "r") as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, pos, dep = line.rstrip("\n").split("\t")
            key = (chrom, int(pos))
            gname = gene_by_key.get(key)
            if gname is not None:
                per_gene_depths[gname].append(int(dep))
    rows: List[Dict[str, object]] = []
    for g in genes:
        vals = per_gene_depths[g.name]
        length = g.length
        if len(vals) < length:
            vals = vals + [0] * (length - len(vals))
        mean_dp = sum(vals) / length if length > 0 else 0.0
        med_dp = statistics.median(vals) if vals else 0.0
        max_dp = max(vals) if vals else 0
        breadth1 = (sum(1 for x in vals if x >= 1) / length * 100.0) if length > 0 else 0.0
        breadth5 = (sum(1 for x in vals if x >= 5) / length * 100.0) if length > 0 else 0.0
        breadth10 = (sum(1 for x in vals if x >= 10) / length * 100.0) if length > 0 else 0.0
        rows.append({
            "gene": g.name,
            "chrom": g.chrom,
            "start_1based": g.start1,
            "end_1based": g.end,
            "length_bp": length,
            "mean_depth": round(mean_dp, 4),
            "median_depth": round(float(med_dp), 4),
            "max_depth": int(max_dp),
            "breadth_1x_pct": round(breadth1, 4),
            "breadth_5x_pct": round(breadth5, 4),
            "breadth_10x_pct": round(breadth10, 4),
        })
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return rows


def gene_presence_calls(rows: List[Dict[str, object]], cfg: Config, out_tsv: Path, sample: str) -> List[Dict[str, object]]:
    out_rows: List[Dict[str, object]] = []
    for r in rows:
        present = (
            float(r["mean_depth"]) >= cfg.gene_presence_min_mean_depth and
            float(r["breadth_1x_pct"]) >= cfg.gene_presence_min_breadth_1x and
            float(r["breadth_5x_pct"]) >= cfg.gene_presence_min_breadth_5x and
            float(r["breadth_10x_pct"]) >= cfg.gene_presence_min_breadth_10x
        )
        out_rows.append({
            "sample": sample,
            "gene": r["gene"],
            "chrom": r["chrom"],
            "start_1based": r["start_1based"],
            "end_1based": r["end_1based"],
            "length_bp": r["length_bp"],
            "mean_depth": r["mean_depth"],
            "median_depth": r["median_depth"],
            "max_depth": r["max_depth"],
            "breadth_1x_pct": r["breadth_1x_pct"],
            "breadth_5x_pct": r["breadth_5x_pct"],
            "breadth_10x_pct": r["breadth_10x_pct"],
            "gene_status": "present" if present else "deleted_or_not_detected",
        })
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    return out_rows


def call_variants(cfg: Config, bam: Path, vcf_gz: Path, logger: logging.Logger) -> None:
    ensure_dir(vcf_gz.parent)
    script = (
        f"set -euo pipefail; "
        f"{cfg.bcftools_bin} mpileup -Ou -f {shlex_q(cfg.pf_ref)} "
        f"-q {cfg.min_mapq} -Q {cfg.min_base_qual} "
        f"{shlex_q(bam)} "
        f"| {cfg.bcftools_bin} call -mv -Ou "
        f"| {cfg.bcftools_bin} filter -Oz "
        f"-i 'QUAL>={cfg.variant_min_qual} && INFO/DP>={cfg.variant_min_dp}' "
        f"-o {shlex_q(vcf_gz)}"
    )
    run_bash(script, logger)
    run_cmd([cfg.tabix_bin, "-f", "-p", "vcf", str(vcf_gz)], logger)


def snpeff_annotate(cfg: Config, vcf_in: Path, vcf_out: Path, logger: logging.Logger) -> None:
    ensure_dir(vcf_out.parent)
    base = [
        "java",
        "-Xmx8g",
        "-jar",
        str(cfg.snpeff_jar),
    ]
    if cfg.snpeff_config:
        base += ["-c", str(cfg.snpeff_config)]
    base += [cfg.snpeff_db, str(vcf_in)]
    script = (
        f"set -euo pipefail; "
        + " ".join(shlex_q(x) if i >= 3 else x for i, x in enumerate(base))
        + f" | {cfg.bgzip_bin} -c > {shlex_q(vcf_out)}"
    )
    run_bash(script, logger)
    run_cmd([cfg.tabix_bin, "-f", "-p", "vcf", str(vcf_out)], logger)


def parse_fastp_json(p: Path) -> Dict[str, object]:
    with open(p, "r") as fh:
        return json.load(fh)


def safe_num(d: Dict[str, object], *keys: str) -> Optional[object]:
    cur = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return None
        cur = cur[k]
    return cur


def summarize_qc(sample: str, fastp_json: Path, out_tsv: Path) -> Dict[str, object]:
    data = parse_fastp_json(fastp_json)
    rows = {
        "sample": sample,
        "raw_total_reads": safe_num(data, "summary", "before_filtering", "total_reads"),
        "raw_total_bases": safe_num(data, "summary", "before_filtering", "total_bases"),
        "raw_q20_bases": safe_num(data, "summary", "before_filtering", "q20_bases"),
        "raw_q30_bases": safe_num(data, "summary", "before_filtering", "q30_bases"),
        "raw_gc_content": safe_num(data, "summary", "before_filtering", "gc_content"),
        "clean_total_reads": safe_num(data, "summary", "after_filtering", "total_reads"),
        "clean_total_bases": safe_num(data, "summary", "after_filtering", "total_bases"),
        "clean_q20_bases": safe_num(data, "summary", "after_filtering", "q20_bases"),
        "clean_q30_bases": safe_num(data, "summary", "after_filtering", "q30_bases"),
        "clean_gc_content": safe_num(data, "summary", "after_filtering", "gc_content"),
        "reads_passed_filtering": safe_num(data, "filtering_result", "passed_filter_reads"),
        "reads_low_quality": safe_num(data, "filtering_result", "low_quality_reads"),
        "reads_too_many_N": safe_num(data, "filtering_result", "too_many_N_reads"),
        "reads_too_short": safe_num(data, "filtering_result", "too_short_reads"),
        "reads_with_too_long_polyX": safe_num(data, "filtering_result", "too_long_polyx_reads"),
    }
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows.keys()), delimiter="\t")
        w.writeheader()
        w.writerow(rows)
    return rows


def count_primary_reads_in_bam(bam: Path) -> int:
    total = 0
    with pysam.AlignmentFile(str(bam), "rb") as bf:
        for rec in bf.fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary or rec.is_duplicate:
                continue
            total += 1
    return total


def region_variants(vcf_gz: Path, region: GeneRegion) -> List[Tuple[str, int, str, str, str]]:
    out = []
    with pysam.VariantFile(str(vcf_gz)) as vf:
        for rec in vf.fetch(region.chrom, region.start0, region.end):
            if rec.alts is None or len(rec.alts) != 1:
                continue
            ref = str(rec.ref).upper()
            alt = str(rec.alts[0]).upper()
            if len(ref) != 1 or len(alt) != 1:
                continue
            ann = ""
            try:
                x = rec.info.get("ANN", None)
                if x is not None:
                    ann = ",".join(map(str, x)) if isinstance(x, (list, tuple)) else str(x)
            except Exception:
                ann = ""
            out.append((region.chrom, int(rec.pos), ref, alt, ann))
    return out


def extract_gene_haplotypes(
    bam_path: Path,
    region: GeneRegion,
    snps: List[Tuple[str, int, str, str, str]],
    sample: str,
    gene_name: str,
    out_reads_tsv: Path,
    out_haps_tsv: Path,
    out_fasta: Path,
    cfg: Config,
) -> List[Dict[str, object]]:
    out_rows: List[Dict[str, object]] = []
    if not snps:
        with open(out_reads_tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample", "gene", "read_name", "haplotype_string", "informative_sites"])
        with open(out_haps_tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample", "gene", "haplotype_id", "haplotype_string", "supporting_reads", "frequency", "informative_sites"])
        with open(out_fasta, "w") as fh:
            pass
        return out_rows

    snp_index = {(chrom, pos): idx for idx, (chrom, pos, ref, alt, ann) in enumerate(snps)}
    ref_alt = {(chrom, pos): (ref, alt) for chrom, pos, ref, alt, ann in snps}
    read_calls: Dict[str, List[str]] = {}
    read_info_count: Dict[str, int] = defaultdict(int)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for aln in bam.fetch(region.chrom, region.start0, region.end):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary or aln.is_duplicate:
                continue
            if aln.mapping_quality < cfg.hap_min_mapq:
                continue
            if aln.query_sequence is None:
                continue
            qname = aln.query_name
            if qname not in read_calls:
                read_calls[qname] = ["N"] * len(snps)
            pairs = aln.get_aligned_pairs(matches_only=False)
            seq = aln.query_sequence.upper()
            quals = aln.query_qualities or []
            for qpos, rpos in pairs:
                if qpos is None or rpos is None:
                    continue
                key = (region.chrom, rpos + 1)
                idx = snp_index.get(key)
                if idx is None:
                    continue
                if qpos >= len(seq):
                    continue
                if qpos >= len(quals):
                    continue
                b = seq[qpos]
                q = int(quals[qpos])
                if q < cfg.hap_min_baseq:
                    continue
                ref, alt = ref_alt[key]
                if b == ref or b == alt:
                    if read_calls[qname][idx] == "N":
                        read_info_count[qname] += 1
                    read_calls[qname][idx] = b

    with open(out_reads_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "gene", "read_name", "haplotype_string", "informative_sites"])

        filtered_haps: List[str] = []
        for qname, hap in read_calls.items():
            h = "".join(hap)
            ninfo = sum(1 for x in hap if x != "N")
            if ninfo < cfg.hap_min_informative_sites:
                continue
            filtered_haps.append(h)
            w.writerow([sample, gene_name, qname, h, ninfo])

    counts = Counter(filtered_haps)
    total = sum(counts.values())

    with open(out_haps_tsv, "w", newline="") as fh, open(out_fasta, "w") as fa:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "gene", "haplotype_id", "haplotype_string", "supporting_reads", "frequency", "informative_sites"])
        idx = 0
        for hap, cnt in counts.most_common():
            if cnt < cfg.hap_min_reads:
                continue
            idx += 1
            hid = f"{normalize_gene_name(gene_name)}_HAP{idx}"
            freq = cnt / total if total > 0 else 0.0
            ninfo = sum(1 for x in hap if x != "N")
            row = {
                "sample": sample,
                "gene": gene_name,
                "haplotype_id": hid,
                "haplotype_string": hap,
                "supporting_reads": cnt,
                "frequency": round(freq, 6),
                "informative_sites": ninfo,
            }
            out_rows.append(row)
            w.writerow(row)
            fa.write(f">{sample}|{gene_name}|{hid}|reads={cnt}|freq={freq:.6f}|informative={ninfo}\n")
            fa.write(hap + "\n")

    return out_rows


def write_json(obj: object, path: Path) -> None:
    ensure_dir(path.parent)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2)


def aggregate_sample_gene_calls(sample_call_files: List[Path], out_tsv: Path, out_prev_tsv: Path) -> None:
    rows: List[Dict[str, str]] = []
    for fp in sample_call_files:
        with open(fp, "r") as fh:
            dr = csv.DictReader(fh, delimiter="\t")
            rows.extend(list(dr))
    if not rows:
        return
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    by_gene: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for r in rows:
        by_gene[r["gene"]].append(r)

    prev_rows = []
    for gene, vals in sorted(by_gene.items()):
        n = len(vals)
        present = sum(1 for x in vals if x["gene_status"] == "present")
        deleted = n - present
        prev_rows.append({
            "gene": gene,
            "n_samples": n,
            "n_present": present,
            "n_deleted_or_not_detected": deleted,
            "prevalence_present_pct": round((present / n * 100.0) if n > 0 else 0.0, 4),
            "prevalence_deleted_or_not_detected_pct": round((deleted / n * 100.0) if n > 0 else 0.0, 4),
        })
    with open(out_prev_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(prev_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(prev_rows)


def aggregate_qc(sample_qc_files: List[Path], out_tsv: Path) -> None:
    rows: List[Dict[str, str]] = []
    for fp in sample_qc_files:
        with open(fp, "r") as fh:
            dr = csv.DictReader(fh, delimiter="\t")
            rows.extend(list(dr))
    if not rows:
        return
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def aggregate_haplotypes(sample_hap_files: List[Path], out_tsv: Path, out_prev_tsv: Path) -> None:
    rows: List[Dict[str, str]] = []
    for fp in sample_hap_files:
        if not fp.exists():
            continue
        with open(fp, "r") as fh:
            dr = csv.DictReader(fh, delimiter="\t")
            rows.extend(list(dr))
    if not rows:
        with open(out_tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample", "gene", "haplotype_id", "haplotype_string", "supporting_reads", "frequency", "informative_sites"])
        with open(out_prev_tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["gene", "haplotype_string", "n_samples", "total_supporting_reads", "sample_prevalence_pct", "read_weighted_fraction_within_gene_pct"])
        return

    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    grouped: Dict[Tuple[str, str], Dict[str, object]] = {}
    total_reads_by_gene: Dict[str, int] = defaultdict(int)
    samples_by_gene: Dict[str, set] = defaultdict(set)

    for r in rows:
        gene = r["gene"]
        hap = r["haplotype_string"]
        samp = r["sample"]
        support = int(float(r["supporting_reads"]))
        total_reads_by_gene[gene] += support
        samples_by_gene[gene].add(samp)
        key = (gene, hap)
        if key not in grouped:
            grouped[key] = {
                "gene": gene,
                "haplotype_string": hap,
                "samples": set(),
                "total_supporting_reads": 0,
            }
        grouped[key]["samples"].add(samp)
        grouped[key]["total_supporting_reads"] += support

    prev_rows = []
    for (gene, hap), d in sorted(grouped.items()):
        n_gene_samples = len(samples_by_gene[gene])
        n_hap_samples = len(d["samples"])
        total_gene_reads = total_reads_by_gene[gene]
        prev_rows.append({
            "gene": gene,
            "haplotype_string": hap,
            "n_samples": n_hap_samples,
            "total_supporting_reads": d["total_supporting_reads"],
            "sample_prevalence_pct": round((n_hap_samples / n_gene_samples * 100.0) if n_gene_samples > 0 else 0.0, 4),
            "read_weighted_fraction_within_gene_pct": round((d["total_supporting_reads"] / total_gene_reads * 100.0) if total_gene_reads > 0 else 0.0, 4),
        })

    with open(out_prev_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(prev_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(prev_rows)


def aggregate_variant_counts(sample_variant_files: List[Path], out_tsv: Path) -> None:
    rows = []
    for fp in sample_variant_files:
        if not fp.exists():
            continue
        with open(fp, "r") as fh:
            dr = csv.DictReader(fh, delimiter="\t")
            rows.extend(list(dr))
    if not rows:
        return
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def write_region_variant_table(sample: str, vcf_gz: Path, genes: List[GeneRegion], out_tsv: Path) -> None:
    rows = []
    with pysam.VariantFile(str(vcf_gz)) as vf:
        for g in genes:
            for rec in vf.fetch(g.chrom, g.start0, g.end):
                if rec.alts is None:
                    continue
                ann = ""
                try:
                    x = rec.info.get("ANN", None)
                    if x is not None:
                        ann = ",".join(map(str, x)) if isinstance(x, (list, tuple)) else str(x)
                except Exception:
                    ann = ""
                for alt in rec.alts:
                    rows.append({
                        "sample": sample,
                        "gene": g.name,
                        "chrom": rec.contig,
                        "pos": int(rec.pos),
                        "ref": str(rec.ref),
                        "alt": str(alt),
                        "qual": "." if rec.qual is None else rec.qual,
                        "annotation": ann,
                    })
    if not rows:
        with open(out_tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample", "gene", "chrom", "pos", "ref", "alt", "qual", "annotation"])
        return
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def process_sample(cfg: Config, sample_pair: SamplePair, genes: List[GeneRegion]) -> Dict[str, str]:
    sample = sample_pair.sample
    logger = get_logger(sample)
    sample_dir = cfg.output_dir / "samples" / sample
    raw_qc_dir = sample_dir / "qc_raw"
    clean_qc_dir = sample_dir / "qc_clean"
    trim_dir = sample_dir / "trimmed"
    host_dir = sample_dir / "host_depletion"
    pf_dir = sample_dir / "pf_mapping"
    cov_dir = sample_dir / "coverage"
    var_dir = sample_dir / "variants"
    hap_dir = sample_dir / "haplotypes"
    sum_dir = sample_dir / "summary"
    tmp_dir = sample_dir / "tmp"

    for d in [raw_qc_dir, clean_qc_dir, trim_dir, host_dir, pf_dir, cov_dir, var_dir, hap_dir, sum_dir, tmp_dir]:
        ensure_dir(d)

    print(step(f"{sample} :: QC raw"))
    fastqc_on_files(cfg, [sample_pair.r1, sample_pair.r2], raw_qc_dir, logger)

    trim_r1 = trim_dir / f"{sample}.trimmed.R1.fastq.gz"
    trim_r2 = trim_dir / f"{sample}.trimmed.R2.fastq.gz"
    fastp_html = trim_dir / f"{sample}.fastp.html"
    fastp_json = trim_dir / f"{sample}.fastp.json"

    print(step(f"{sample} :: Trim/filter"))
    run_fastp(cfg, sample_pair.r1, sample_pair.r2, trim_r1, trim_r2, fastp_html, fastp_json, logger)

    print(step(f"{sample} :: QC clean"))
    fastqc_on_files(cfg, [trim_r1, trim_r2], clean_qc_dir, logger)

    human_bam = host_dir / f"{sample}.human.bam"
    human_flagstat = host_dir / f"{sample}.human.flagstat.txt"
    human_idxstats = host_dir / f"{sample}.human.idxstats.txt"
    human_stats = host_dir / f"{sample}.human.stats.txt"
    nonhuman_r1 = host_dir / f"{sample}.nonhuman.R1.fastq"
    nonhuman_r2 = host_dir / f"{sample}.nonhuman.R2.fastq"
    nonhuman_single = host_dir / f"{sample}.nonhuman.singletons.fastq"

    print(step(f"{sample} :: Human decontamination"))
    map_paired_bwa(cfg, cfg.human_ref, trim_r1, trim_r2, human_bam, logger)
    bam_stats(cfg, human_bam, human_flagstat, human_idxstats, human_stats, logger)
    ensure_samtools_fastq_extract_nonhuman(cfg.samtools_bin, human_bam, nonhuman_r1, nonhuman_r2, nonhuman_single, logger)

    pf_bam = pf_dir / f"{sample}.Pf3D7.bam"
    pf_flagstat = pf_dir / f"{sample}.Pf3D7.flagstat.txt"
    pf_idxstats = pf_dir / f"{sample}.Pf3D7.idxstats.txt"
    pf_stats = pf_dir / f"{sample}.Pf3D7.stats.txt"

    print(step(f"{sample} :: Map to Pf3D7"))
    map_paired_bwa(cfg, cfg.pf_ref, nonhuman_r1, nonhuman_r2, pf_bam, logger)
    bam_stats(cfg, pf_bam, pf_flagstat, pf_idxstats, pf_stats, logger)

    print(step(f"{sample} :: Gene coverage"))
    gene_depth_tsv = cov_dir / f"{sample}.gene_depth.tsv"
    gene_cov_tsv = cov_dir / f"{sample}.gene_coverage.tsv"
    gene_calls_tsv = cov_dir / f"{sample}.gene_calls.tsv"
    depth_bed(cfg, pf_bam, cfg.genes_bed, gene_depth_tsv, logger)
    cov_rows = compute_gene_coverage(gene_depth_tsv, genes, gene_cov_tsv)
    call_rows = gene_presence_calls(cov_rows, cfg, gene_calls_tsv, sample)

    print(step(f"{sample} :: Variant calling"))
    raw_vcf = var_dir / f"{sample}.raw.filtered.vcf.gz"
    anno_vcf = var_dir / f"{sample}.raw.filtered.ann.vcf.gz"
    region_variants_tsv = var_dir / f"{sample}.hrp2_hrp3.variants.tsv"
    call_variants(cfg, pf_bam, raw_vcf, logger)
    snpeff_annotate(cfg, raw_vcf, anno_vcf, logger)
    write_region_variant_table(sample, anno_vcf, genes, region_variants_tsv)

    print(step(f"{sample} :: Haplotype extraction"))
    all_hap_rows: List[Dict[str, object]] = []
    hap_summary_rows: List[Dict[str, object]] = []
    for g in genes:
        gsafe = normalize_gene_name(g.name)
        gene_reads_tsv = hap_dir / f"{sample}.{gsafe}.read_haplotypes.tsv"
        gene_haps_tsv = hap_dir / f"{sample}.{gsafe}.haplotypes.tsv"
        gene_fasta = hap_dir / f"{sample}.{gsafe}.haplotypes.fasta"
        snps = region_variants(anno_vcf, g)
        hap_rows = extract_gene_haplotypes(
            bam_path=pf_bam,
            region=g,
            snps=snps,
            sample=sample,
            gene_name=g.name,
            out_reads_tsv=gene_reads_tsv,
            out_haps_tsv=gene_haps_tsv,
            out_fasta=gene_fasta,
            cfg=cfg,
        )
        all_hap_rows.extend(hap_rows)
        status = next((x["gene_status"] for x in call_rows if x["gene"] == g.name), "deleted_or_not_detected")
        hap_summary_rows.append({
            "sample": sample,
            "gene": g.name,
            "gene_status": status,
            "n_gene_variants": len(snps),
            "n_haplotypes_retained": len(hap_rows),
            "total_haplotype_supporting_reads": sum(int(x["supporting_reads"]) for x in hap_rows) if hap_rows else 0,
        })

    sample_haps_tsv = hap_dir / f"{sample}.all_genes.haplotypes.tsv"
    if all_hap_rows:
        with open(sample_haps_tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(all_hap_rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(all_hap_rows)
    else:
        with open(sample_haps_tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["sample", "gene", "haplotype_id", "haplotype_string", "supporting_reads", "frequency", "informative_sites"])

    hap_summary_tsv = hap_dir / f"{sample}.haplotype_summary.tsv"
    with open(hap_summary_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(hap_summary_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(hap_summary_rows)

    print(step(f"{sample} :: Summaries"))
    qc_summary_tsv = sum_dir / f"{sample}.qc_summary.tsv"
    qc_summary = summarize_qc(sample, fastp_json, qc_summary_tsv)

    pf_primary = count_primary_reads_in_bam(pf_bam)
    human_primary = count_primary_reads_in_bam(human_bam)
    sample_json = sum_dir / f"{sample}.summary.json"
    write_json({
        "sample": sample,
        "inputs": {"r1": str(sample_pair.r1), "r2": str(sample_pair.r2)},
        "files": {
            "trimmed_r1": str(trim_r1),
            "trimmed_r2": str(trim_r2),
            "nonhuman_r1": str(nonhuman_r1),
            "nonhuman_r2": str(nonhuman_r2),
            "human_bam": str(human_bam),
            "pf_bam": str(pf_bam),
            "gene_calls_tsv": str(gene_calls_tsv),
            "gene_coverage_tsv": str(gene_cov_tsv),
            "region_variants_tsv": str(region_variants_tsv),
            "annotated_vcf": str(anno_vcf),
            "sample_haplotypes_tsv": str(sample_haps_tsv),
            "sample_haplotype_summary_tsv": str(hap_summary_tsv),
            "qc_summary_tsv": str(qc_summary_tsv),
        },
        "counts": {
            "human_primary_alignments": human_primary,
            "pf_primary_alignments": pf_primary,
        },
        "qc": qc_summary,
        "gene_calls": call_rows,
        "haplotype_summary": hap_summary_rows,
    }, sample_json)

    if not cfg.keep_tmp:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    print(done(f"{sample} :: complete"))
    return {
        "sample": sample,
        "gene_calls": str(gene_calls_tsv),
        "qc_summary": str(qc_summary_tsv),
        "sample_haplotypes": str(sample_haps_tsv),
        "region_variants": str(region_variants_tsv),
    }


def build_output_layout(root: Path) -> None:
    for p in [
        root / "samples",
        root / "cohort" / "coverage",
        root / "cohort" / "haplotypes",
        root / "cohort" / "qc",
        root / "cohort" / "variants",
        root / "logs",
    ]:
        ensure_dir(p)


def validate_bed_against_reference(genes: List[GeneRegion], fai_dict: Dict[str, int]) -> None:
    for g in genes:
        if g.chrom not in fai_dict:
            raise PipelineError(f"BED contig not found in reference: {g.chrom}")
        if g.start0 < 0 or g.end > fai_dict[g.chrom] or g.start0 >= g.end:
            raise PipelineError(f"Invalid BED interval: {g.chrom}:{g.start1}-{g.end}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--input_dir", required=True, type=Path)
    p.add_argument("--output_dir", required=True, type=Path)
    p.add_argument("--pf_ref", required=True, type=Path)
    p.add_argument("--human_ref", required=True, type=Path)
    p.add_argument("--genes_bed", required=True, type=Path)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--sample_workers", type=int, default=2)
    p.add_argument("--fastp_bin", default="fastp")
    p.add_argument("--fastqc_bin", default="fastqc")
    p.add_argument("--bwa_bin", default="bwa")
    p.add_argument("--samtools_bin", default="samtools")
    p.add_argument("--bedtools_bin", default="bedtools")
    p.add_argument("--bcftools_bin", default="bcftools")
    p.add_argument("--bgzip_bin", default="bgzip")
    p.add_argument("--tabix_bin", default="tabix")
    p.add_argument("--snpeff_jar", required=True, type=Path)
    p.add_argument("--snpeff_db", default="Pf3D7_v2")
    p.add_argument("--snpeff_config", type=Path, default=None)
    p.add_argument("--min_read_len", type=int, default=50)
    p.add_argument("--min_base_qual", type=int, default=20)
    p.add_argument("--min_mapq", type=int, default=20)
    p.add_argument("--variant_min_dp", type=int, default=10)
    p.add_argument("--variant_min_qual", type=int, default=30)
    p.add_argument("--gene_presence_min_mean_depth", type=float, default=5.0)
    p.add_argument("--gene_presence_min_breadth_1x", type=float, default=90.0)
    p.add_argument("--gene_presence_min_breadth_5x", type=float, default=80.0)
    p.add_argument("--gene_presence_min_breadth_10x", type=float, default=70.0)
    p.add_argument("--hap_min_reads", type=int, default=3)
    p.add_argument("--hap_min_informative_sites", type=int, default=2)
    p.add_argument("--hap_min_mapq", type=int, default=20)
    p.add_argument("--hap_min_baseq", type=int, default=20)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--keep_tmp", action="store_true")
    p.add_argument("--log_file", type=Path, default=Path("pipeline.log"))
    return p.parse_args()


def banner(cfg: Config, n_samples: int) -> None:
    print(f"{MAGENTA}{BOLD}======================================================================{RESET}")
    print(f"{MAGENTA}{BOLD} PF-HRP2/HRP3 Illumina Pipeline{RESET}")
    print(f"{MAGENTA}{BOLD}======================================================================{RESET}")
    print(info(f"Input directory        : {cfg.input_dir}"))
    print(info(f"Output directory       : {cfg.output_dir}"))
    print(info(f"P. falciparum ref      : {cfg.pf_ref}"))
    print(info(f"Human ref              : {cfg.human_ref}"))
    print(info(f"Genes BED              : {cfg.genes_bed}"))
    print(info(f"snpEff DB              : {cfg.snpeff_db}"))
    print(info(f"Threads                : {cfg.threads}"))
    print(info(f"Sample workers         : {cfg.sample_workers}"))
    print(info(f"Samples detected       : {n_samples}"))
    print(f"{MAGENTA}{BOLD}======================================================================{RESET}")


def main() -> int:
    args = parse_args()
    cfg = Config(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        pf_ref=args.pf_ref,
        human_ref=args.human_ref,
        genes_bed=args.genes_bed,
        threads=args.threads,
        sample_workers=args.sample_workers,
        fastp_bin=args.fastp_bin,
        fastqc_bin=args.fastqc_bin,
        bwa_bin=args.bwa_bin,
        samtools_bin=args.samtools_bin,
        bedtools_bin=args.bedtools_bin,
        bcftools_bin=args.bcftools_bin,
        bgzip_bin=args.bgzip_bin,
        tabix_bin=args.tabix_bin,
        snpeff_jar=args.snpeff_jar,
        snpeff_db=args.snpeff_db,
        snpeff_config=args.snpeff_config,
        min_read_len=args.min_read_len,
        min_base_qual=args.min_base_qual,
        min_mapq=args.min_mapq,
        variant_min_dp=args.variant_min_dp,
        variant_min_qual=args.variant_min_qual,
        gene_presence_min_mean_depth=args.gene_presence_min_mean_depth,
        gene_presence_min_breadth_1x=args.gene_presence_min_breadth_1x,
        gene_presence_min_breadth_5x=args.gene_presence_min_breadth_5x,
        gene_presence_min_breadth_10x=args.gene_presence_min_breadth_10x,
        hap_min_reads=args.hap_min_reads,
        hap_min_informative_sites=args.hap_min_informative_sites,
        hap_min_mapq=args.hap_min_mapq,
        hap_min_baseq=args.hap_min_baseq,
        resume=args.resume,
        keep_tmp=args.keep_tmp,
        log_file=args.log_file,
    )

    build_output_layout(cfg.output_dir)
    root_logger = setup_logger(cfg.output_dir / "logs" / cfg.log_file.name)
    logger = get_logger("main")

    try:
        check_dependencies(cfg)
        index_reference_bwa(cfg.pf_ref, cfg.bwa_bin, cfg.samtools_bin, logger)
        index_reference_bwa(cfg.human_ref, cfg.bwa_bin, cfg.samtools_bin, logger)

        genes = parse_bed(cfg.genes_bed)
        validate_bed_against_reference(genes, fai_contigs(Path(str(cfg.pf_ref) + ".fai")))

        pairs = discover_pairs(cfg.input_dir)
        banner(cfg, len(pairs))

        results: List[Dict[str, str]] = []
        workers = max(1, min(cfg.sample_workers, len(pairs)))

        if workers == 1:
            for sp in pairs:
                results.append(process_sample(cfg, sp, genes))
        else:
            with cf.ThreadPoolExecutor(max_workers=workers) as ex:
                futs = {ex.submit(process_sample, cfg, sp, genes): sp.sample for sp in pairs}
                for fut in cf.as_completed(futs):
                    sample = futs[fut]
                    try:
                        results.append(fut.result())
                    except Exception as e:
                        logger.error("Sample failed: %s", sample)
                        logger.error("TRACEBACK:\n%s", traceback.format_exc())
                        print(err(f"{sample} :: {e}"))

        gene_call_files = [Path(x["gene_calls"]) for x in results if "gene_calls" in x]
        qc_files = [Path(x["qc_summary"]) for x in results if "qc_summary" in x]
        hap_files = [Path(x["sample_haplotypes"]) for x in results if "sample_haplotypes" in x]
        var_files = [Path(x["region_variants"]) for x in results if "region_variants" in x]

        print(step("Cohort :: aggregate sample-level gene calls"))
        aggregate_sample_gene_calls(
            gene_call_files,
            cfg.output_dir / "cohort" / "coverage" / "all_samples.gene_calls.tsv",
            cfg.output_dir / "cohort" / "coverage" / "gene_prevalence.tsv",
        )

        print(step("Cohort :: aggregate QC"))
        aggregate_qc(
            qc_files,
            cfg.output_dir / "cohort" / "qc" / "all_samples.qc_summary.tsv",
        )

        print(step("Cohort :: aggregate haplotypes"))
        aggregate_haplotypes(
            hap_files,
            cfg.output_dir / "cohort" / "haplotypes" / "all_samples.haplotypes.tsv",
            cfg.output_dir / "cohort" / "haplotypes" / "haplotype_prevalence.tsv",
        )

        print(step("Cohort :: aggregate gene-region variants"))
        aggregate_variant_counts(
            var_files,
            cfg.output_dir / "cohort" / "variants" / "all_samples.hrp2_hrp3.variants.tsv",
        )

        cohort_summary = {
            "n_samples_detected": len(pairs),
            "n_samples_completed": len(results),
            "genes_bed": str(cfg.genes_bed),
            "pf_ref": str(cfg.pf_ref),
            "human_ref": str(cfg.human_ref),
            "snpeff_db": cfg.snpeff_db,
        }
        write_json(cohort_summary, cfg.output_dir / "cohort" / "cohort_summary.json")

        print(done("Pipeline completed"))
        print(info(f"Log: {cfg.output_dir / 'logs' / cfg.log_file.name}"))
        return 0

    except Exception as e:
        logger.error("Pipeline failed: %s", e)
        logger.error("TRACEBACK:\n%s", traceback.format_exc())
        print(err(str(e)))
        return 1


if __name__ == "__main__":
    sys.exit(main())