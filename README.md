# pf-hrp2-hrp3-deletion-pipeline

## Overview

`pf-hrp2-hrp3-deletion-pipeline` is an Illumina workflow for the genomic analysis of *Plasmodium falciparum* HRP2 and HRP3 targets.

It supports:

- raw read quality control
- read trimming and filtering
- host read depletion against the human reference genome
- mapping to the *P. falciparum* 3D7 reference genome
- target gene coverage assessment
- HRP2 and HRP3 presence or deletion inference
- variant calling in target regions
- functional annotation with snpEff
- read-based haplotype extraction
- sample-level and cohort-level summary tables

The pipeline is designed for paired-end Illumina FASTQ data and accepts `.fastq`, `.fq`, `.fastq.gz`, and `.fq.gz` files.

---

## Project Structure

```text
pf-hrp2-hrp3-deletion-pipeline/
├── assets/
│   └── hrp2_hrp3_genes.bed
├── scripts/
│   └── pf_hrp2_hrp3_pipeline.py
├── scripts_plot/
│   └── pf_hrp2_hrp3_pipeline_analysis.ipynb
├── results/
├── main.sh
└── README.md
````

---

## Input Requirements

### Sequencing data

The pipeline expects paired-end Illumina reads.

Accepted file patterns include:

* `sample_R1.fastq.gz` and `sample_R2.fastq.gz`
* `sample_R1_001.fastq.gz` and `sample_R2_001.fastq.gz`
* `sample.R1.fastq.gz` and `sample.R2.fastq.gz`
* `sample_1.fq.gz` and `sample_2.fq.gz`

The pipeline automatically detects pairs from the input directory.

### Required reference files

You must provide:

* a *P. falciparum* reference FASTA
* a human reference FASTA for host depletion
* a BED file containing HRP2 and HRP3 target intervals
* a valid `snpEff.jar`
* a valid `snpEff.config`

---

## Target BED File

Example BED file used in this project:

```bed
Pf3D7_08_v3	1373211	1376988	PF3D7_0831800|HRP2	0	-
Pf3D7_13_v3	2840235	2842840	PF3D7_1372200|HRP3	0	-
```

---

## Main Analysis Steps

### 1. Raw read quality control

The pipeline runs FastQC on raw paired-end reads.

Outputs:

* raw FastQC reports

### 2. Read trimming and filtering

The pipeline runs fastp to remove adapters and low-quality sequence.

Outputs:

* trimmed R1 and R2 FASTQ files
* fastp HTML report
* fastp JSON report

### 3. Post-trimming quality control

The pipeline runs FastQC again on cleaned reads.

Outputs:

* cleaned FastQC reports

### 4. Human read depletion

Trimmed reads are mapped to the human reference genome with BWA.
Reads not assigned to the host are extracted for downstream analysis.

Outputs:

* host BAM file
* host alignment statistics
* non-human R1 and R2 FASTQ files
* singleton FASTQ file

### 5. Mapping to *P. falciparum* 3D7

Non-human reads are mapped to the *P. falciparum* reference genome.

Outputs:

* parasite BAM file
* BAM index
* mapping statistics

### 6. HRP2 and HRP3 coverage analysis

The pipeline computes per-base depth across the BED targets and summarizes coverage by gene.

Outputs:

* per-base gene depth table
* per-gene coverage summary
* gene presence or deletion call table

### 7. Variant calling

Variants are called from the parasite BAM with bcftools.

Outputs:

* filtered VCF
* indexed VCF

### 8. Variant annotation

Filtered variants are annotated with snpEff.

Outputs:

* annotated VCF
* indexed annotated VCF
* HRP2 and HRP3 regional variant table

### 9. Haplotype extraction

Read-supported haplotypes are reconstructed from informative SNP positions within each target gene.

Outputs:

* read haplotype table
* haplotype summary table
* haplotype FASTA file
* pooled sample haplotype table

### 10. Cohort-level aggregation

The pipeline merges per-sample results into cohort-wide summaries.

Outputs:

* all-sample gene call table
* gene prevalence table
* all-sample QC summary
* all-sample haplotype table
* haplotype prevalence table
* all-sample regional variant table
* cohort summary JSON

---

## Output Layout

```text
results/
├── logs/
├── samples/
│   └── SAMPLE_ID/
│       ├── qc_raw/
│       ├── qc_clean/
│       ├── trimmed/
│       ├── host_depletion/
│       ├── pf_mapping/
│       ├── coverage/
│       ├── variants/
│       ├── haplotypes/
│       ├── summary/
│       └── tmp/
└── cohort/
    ├── coverage/
    ├── haplotypes/
    ├── qc/
    ├── variants/
    └── cohort_summary.json
```

---

## Key Output Files

### Per-sample files

For each sample, the pipeline generates:

* `coverage/<sample>.gene_depth.tsv`
* `coverage/<sample>.gene_coverage.tsv`
* `coverage/<sample>.gene_calls.tsv`
* `variants/<sample>.raw.filtered.vcf.gz`
* `variants/<sample>.raw.filtered.ann.vcf.gz`
* `variants/<sample>.hrp2_hrp3.variants.tsv`
* `haplotypes/<sample>.*.haplotypes.tsv`
* `haplotypes/<sample>.all_genes.haplotypes.tsv`
* `haplotypes/<sample>.haplotype_summary.tsv`
* `summary/<sample>.qc_summary.tsv`
* `summary/<sample>.summary.json`

### Cohort files

The cohort directory contains:

* `coverage/all_samples.gene_calls.tsv`
* `coverage/gene_prevalence.tsv`
* `qc/all_samples.qc_summary.tsv`
* `haplotypes/all_samples.haplotypes.tsv`
* `haplotypes/haplotype_prevalence.tsv`
* `variants/all_samples.hrp2_hrp3.variants.tsv`
* `cohort_summary.json`

---

## Interpretation of Gene Calls

Gene calls are based on coverage thresholds across each target region.

A gene is classified as `present` when all selected coverage criteria are met.

A gene is classified as `deleted_or_not_detected` when one or more criteria are not met.

This label reflects lack of sufficient target support. It should be interpreted with caution in low-parasitemia or low-coverage samples.

---

## Software Requirements

The pipeline requires the following tools:

* Python 3.11
* fastp
* FastQC
* BWA
* samtools
* bcftools
* bedtools
* htslib
* pysam
* Java
* snpEff

---

## Conda Environment

Example minimal environment:

```yaml
name: pf-hrp2-hrp3
channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - python=3.11
  - fastp
  - fastqc
  - bwa
  - samtools
  - bcftools
  - bedtools
  - htslib
  - pysam
  - snpeff
  - openjdk=21
```

Create and activate the environment:

```bash
conda env create -f environment.yml
conda activate pf-hrp2-hrp3
```

---

## snpEff Setup

The current pipeline uses a JAR-based snpEff call.

You must provide:

* `--snpeff_jar`
* `--snpeff_config`
* `--snpeff_db`

Example:

```bash
SNPEFF_JAR="/path/to/snpEff.jar"
SNPEFF_CONFIG="/path/to/snpEff.config"
SNPEFF_DB="Pf3D7_v2"
```

---

## Running the Pipeline

### Main launcher

The recommended entry point is:

```bash
bash main.sh
```

### Direct Python execution

You can also run the script directly:

```bash
python3 scripts/pf_hrp2_hrp3_pipeline.py \
  --input_dir /path/to/fastq \
  --output_dir /path/to/results \
  --pf_ref /path/to/Pf3D7_v2.fasta \
  --human_ref /path/to/hg38.fa \
  --genes_bed assets/hrp2_hrp3_genes.bed \
  --snpeff_jar /path/to/snpEff.jar \
  --snpeff_db Pf3D7_v2 \
  --snpeff_config /path/to/snpEff.config \
  --threads 8 \
  --sample_workers 1 \
  --bwa_threads_cap 4 \
  --sort_threads_cap 4 \
  --fastp_threads_cap 4 \
  --bcftools_threads_cap 2 \
  --resume
```

---

## Recommended Runtime Settings

Human read depletion against hg38 is the most resource-intensive step.

Recommended starting parameters:

* `threads = 8`
* `sample_workers = 1`
* `bwa_threads_cap = 4`
* `sort_threads_cap = 4`

These settings are safer on shared compute systems and reduce failure caused by thread oversubscription.

---

## Command-Line Arguments

### Required

* `--input_dir`
* `--output_dir`
* `--pf_ref`
* `--human_ref`
* `--genes_bed`
* `--snpeff_jar`

### Common optional arguments

* `--snpeff_db`
* `--snpeff_config`
* `--threads`
* `--sample_workers`
* `--resume`
* `--keep_tmp`
* `--log_file`

### Performance controls

* `--fastqc_threads_cap`
* `--fastp_threads_cap`
* `--bwa_threads_cap`
* `--sort_threads_cap`
* `--bcftools_threads_cap`

### Analysis thresholds

* `--min_read_len`
* `--min_base_qual`
* `--min_mapq`
* `--variant_min_dp`
* `--variant_min_qual`
* `--gene_presence_min_mean_depth`
* `--gene_presence_min_breadth_1x`
* `--gene_presence_min_breadth_5x`
* `--gene_presence_min_breadth_10x`
* `--hap_min_reads`
* `--hap_min_informative_sites`
* `--hap_min_mapq`
* `--hap_min_baseq`

---

## Logging

Pipeline logs are written to:

```text
results/logs/pipeline.log
```

This file captures:

* executed commands
* stdout and stderr
* pipeline failures
* traceback information

---

## Resume Mode

If `--resume` is enabled, the pipeline skips samples that already have a completed summary JSON file.

This is useful after interruption or partial failure.

---

## Limitations

* The pipeline is designed for paired-end Illumina data.
* Gene deletion inference is coverage-based.
* Low coverage can reduce confidence in deletion calls.
* Haplotype reconstruction is limited to informative SNP positions present in the annotated VCF.
* Structural deletion breakpoints are not explicitly resolved.

---

## Recommended Validation

For high-confidence HRP2 and HRP3 deletion studies, it is good practice to complement this workflow with:

* parasite density assessment
* broader genome coverage review
* flanking region analysis
* independent molecular confirmation when required

---

## Citation and Use

If you use this workflow in a report, manuscript, or surveillance analysis, document:

* reference genome build
* BED target definitions
* coverage thresholds
* variant filters
* snpEff database version
* software versions
* runtime parameters

This improves reproducibility and interpretability.

---

## Author

**Khadim GUEYE**

IGH / ACEGID (Institute of Genomics and Global Health)


