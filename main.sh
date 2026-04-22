#!/usr/bin/env bash

###############################################################################
# Project: pf-hrp2-hrp3-deletion-pipeline
# Script:  main.sh
# Author:  Khadim IGH / ACEGID
# Date:    2026-04-22
# Purpose: Run complete Illumina pipeline for HRP2/HRP3 deletion detection
#          - QC (before/after)
#          - Human decontamination
#          - Mapping to Pf3D7
#          - Variant calling + annotation (snpEff)
#          - Gene coverage + deletion detection
#          - Haplotype extraction
#          - Cohort-level prevalence
###############################################################################

set -euo pipefail

############################################
# CONFIGURATION
############################################

# Input FASTQ directory (auto-detects fq/fastq(.gz))
INPUT_DIR="/mnt/hpc_acegid/home/khadmig/work/data/malaria/ILLUMINA_DATA/260420_VH00635_5_AAHL3MNM5/fastq"

# Output directory
OUTPUT_DIR="/mnt/hpc_acegid/home/khadmig/work/pf-hrp2-hrp3-deletion-pipeline/results"

# References
PF_REF="/mnt/hpc_acegid/home/khadmig/work/ref_3d7/Pf3D7_v2.fasta"
HUMAN_REF="/mnt/hpc_acegid/nfsscratch/DATABASE/hg38/hg38.fa"   

# Target genes (HRP2 / HRP3 BED)
GENES_BED="assets/hrp2_hrp3_genes.bed"

# snpEff configuration
SNPEFF_JAR="/mnt/hpc_acegid/home/khadmig/miniconda/envs/pf_ont_pipeline/share/snpeff-5.3.0a-0/snpEff.jar"
SNPEFF_DB="Pf3D7_v2"
SNPEFF_CONFIG="/mnt/hpc_acegid/nfsscratch/DATABASE/snpeff/snpEff.config"

# Resources
THREADS=16
SAMPLE_WORKERS=4

# Pipeline script
SCRIPT="scripts/pf_hrp2_hrp3_pipeline.py"

############################################
# TERMINAL COLORS
############################################

GREEN="\033[32m"
RED="\033[31m"
BLUE="\033[34m"
YELLOW="\033[33m"
RESET="\033[0m"

log()  { echo -e "${BLUE}[RUN]${RESET} $1"; }
ok()   { echo -e "${GREEN}[OK]${RESET} $1"; }
warn() { echo -e "${YELLOW}[WARN]${RESET} $1"; }
fail() { echo -e "${RED}[FAIL]${RESET} $1"; }

############################################
# FUNCTIONS
############################################

#--------------------------------------------------
# Check required inputs before running pipeline
#--------------------------------------------------
check_inputs() {
    log "Checking input files and directories..."

    [ -d "$INPUT_DIR" ] || { fail "Missing INPUT_DIR"; exit 1; }
    [ -f "$PF_REF" ] || { fail "Missing PF_REF"; exit 1; }
    [ -f "$GENES_BED" ] || { fail "Missing GENES_BED"; exit 1; }
    [ -f "$SNPEFF_JAR" ] || { fail "Missing SNPEFF_JAR"; exit 1; }
    [ -f "$SNPEFF_CONFIG" ] || { fail "Missing SNPEFF_CONFIG"; exit 1; }

    if [ ! -f "$HUMAN_REF" ]; then
        warn "HUMAN_REF not found → host decontamination will FAIL"
    fi

    mkdir -p "$OUTPUT_DIR"

    ok "All inputs validated"
}

#--------------------------------------------------
# Run full pipeline
#--------------------------------------------------
run_pipeline() {
    log "Launching HRP2/HRP3 pipeline..."

    python3 "$SCRIPT" \
      --input_dir "$INPUT_DIR" \
      --output_dir "$OUTPUT_DIR" \
      --pf_ref "$PF_REF" \
      --human_ref "$HUMAN_REF" \
      --genes_bed "$GENES_BED" \
      --snpeff_jar "$SNPEFF_JAR" \
      --snpeff_db "$SNPEFF_DB" \
      --snpeff_config "$SNPEFF_CONFIG" \
      --threads "$THREADS" \
      --sample_workers "$SAMPLE_WORKERS"

    ok "Pipeline execution completed"
}

############################################
# MAIN EXECUTION
############################################

main() {
    echo -e "${BLUE}====================================================${RESET}"
    echo -e "${BLUE}  PF HRP2/HRP3 Deletion Pipeline - Illumina${RESET}"
    echo -e "${BLUE}====================================================${RESET}"

    check_inputs
    run_pipeline

    echo -e "${GREEN}====================================================${RESET}"
    echo -e "${GREEN}  ALL DONE ${RESET}"
    echo -e "${GREEN}====================================================${RESET}"
}

main "$@"