#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Project: pf-hrp2-hrp3-deletion-pipeline
# Script:  main.sh
# Author:  Khadim Gueye / ACEGID
# Purpose: Run HRP2/HRP3 deletion pipeline for Illumina or Nanopore data
###############################################################################

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
# USER CONFIGURATION
############################################

# Choose: illumina or nanopore
PLATFORM="${1:-illumina}"

# Common references
PF_REF="/mnt/hpc_acegid/home/khadmig/work/ref_3d7/Pf3D7_v2.fasta"
GENES_BED="assets/hrp2_hrp3_genes.bed"

# snpEff
SNPEFF_JAR="/mnt/hpc_acegid/home/khadmig/miniconda/envs/pf_ont_pipeline/share/snpeff-5.3.0a-0/snpEff.jar"
SNPEFF_DB="Pf3D7_v2"
SNPEFF_CONFIG="/mnt/hpc_acegid/nfsscratch/DATABASE/snpeff/snpEff.config"

# Resources
THREADS=16
SAMPLE_WORKERS=4

############################################
# ILLUMINA CONFIGURATION
############################################

ILLUMINA_INPUT_DIR="/mnt/hpc_acegid/home/khadmig/work/data/malaria/ILLUMINA_DATA/260420_VH00635_5_AAHL3MNM5/fastq"
ILLUMINA_OUTPUT_DIR="/mnt/hpc_acegid/home/khadmig/work/data/malaria/ILLUMINA_DATA/260420_VH00635_5_AAHL3MNM5/results"
HUMAN_REF="/mnt/hpc_acegid/nfsscratch/DATABASE/hg38/hg38.fa"
ILLUMINA_SCRIPT="scripts/pf_hrp2_hrp3_pipeline_illumina.py"

############################################
# NANOPORE CONFIGURATION
############################################

NANOPORE_DEMUX_DIR="/mnt/hpc_acegid/home/khadmig/work/data/malaria/NANOPORE_DATA/demux"
NANOPORE_OUTPUT_DIR="/mnt/hpc_acegid/home/khadmig/work/data/malaria/NANOPORE_DATA/results"
NANOPORE_METADATA="/mnt/hpc_acegid/nfsscratch/khadmig/data/malaria/demux_input/metadata_with_sample.tsv"
NANOPORE_SCRIPT="scripts/pf_hrp2_hrp3_pipeline_nanopore.py"

CLAIR3_MODEL_PATH="$HOME/miniconda/envs/pf_ont_pipeline/bin/models/ont_guppy5"
CLAIR3_PLATFORM="ont"

# Optional DEploid2
RUN_DEPLOID=0
POPULATION_VCF=""
GENETIC_MAP=""
EXCLUDE_VCF=""

############################################
# INPUT VALIDATION
############################################

check_common_inputs() {
    [ -f "$PF_REF" ] || { fail "Missing PF_REF: $PF_REF"; exit 1; }
    [ -f "$GENES_BED" ] || { fail "Missing GENES_BED: $GENES_BED"; exit 1; }
    [ -f "$SNPEFF_JAR" ] || { fail "Missing SNPEFF_JAR: $SNPEFF_JAR"; exit 1; }
    [ -f "$SNPEFF_CONFIG" ] || { fail "Missing SNPEFF_CONFIG: $SNPEFF_CONFIG"; exit 1; }
}

check_illumina_inputs() {
    log "Checking Illumina inputs..."
    check_common_inputs
    [ -d "$ILLUMINA_INPUT_DIR" ] || { fail "Missing ILLUMINA_INPUT_DIR: $ILLUMINA_INPUT_DIR"; exit 1; }
    [ -f "$HUMAN_REF" ] || { fail "Missing HUMAN_REF: $HUMAN_REF"; exit 1; }
    [ -f "$ILLUMINA_SCRIPT" ] || { fail "Missing Illumina script: $ILLUMINA_SCRIPT"; exit 1; }
    mkdir -p "$ILLUMINA_OUTPUT_DIR"
    ok "Illumina inputs validated"
}

check_nanopore_inputs() {
    log "Checking Nanopore inputs..."
    check_common_inputs
    [ -d "$NANOPORE_DEMUX_DIR" ] || { fail "Missing NANOPORE_DEMUX_DIR: $NANOPORE_DEMUX_DIR"; exit 1; }
    [ -f "$NANOPORE_SCRIPT" ] || { fail "Missing Nanopore script: $NANOPORE_SCRIPT"; exit 1; }

    if [ -n "$NANOPORE_METADATA" ] && [ ! -f "$NANOPORE_METADATA" ]; then
        warn "NANOPORE_METADATA not found; barcode names will be used as sample names"
        NANOPORE_METADATA=""
    fi

    if [ ! -d "$CLAIR3_MODEL_PATH" ]; then
        fail "Missing CLAIR3_MODEL_PATH: $CLAIR3_MODEL_PATH"
        exit 1
    fi

    mkdir -p "$NANOPORE_OUTPUT_DIR"
    ok "Nanopore inputs validated"
}

############################################
# PIPELINE RUNNERS
############################################

run_illumina() {
    check_illumina_inputs

    log "Launching Illumina HRP2/HRP3 pipeline..."

    python3 "$ILLUMINA_SCRIPT" \
      --input_dir "$ILLUMINA_INPUT_DIR" \
      --output_dir "$ILLUMINA_OUTPUT_DIR" \
      --pf_ref "$PF_REF" \
      --human_ref "$HUMAN_REF" \
      --genes_bed "$GENES_BED" \
      --snpeff_jar "$SNPEFF_JAR" \
      --snpeff_db "$SNPEFF_DB" \
      --snpeff_config "$SNPEFF_CONFIG" \
      --threads "$THREADS" \
      --sample_workers "$SAMPLE_WORKERS"

    ok "Illumina pipeline completed"
}

run_nanopore() {
    check_nanopore_inputs

    log "Launching Nanopore HRP2/HRP3 pipeline..."

    CMD=(
      python3 "$NANOPORE_SCRIPT"
      --demux_dir "$NANOPORE_DEMUX_DIR"
      --output_dir "$NANOPORE_OUTPUT_DIR"
      --reference "$PF_REF"
      --genes_bed "$GENES_BED"
      --snpeff "$SNPEFF_JAR"
      --snpeff_db "$SNPEFF_DB"
      --snpeff_config "$SNPEFF_CONFIG"
      --threads "$THREADS"
      --sample_workers "$SAMPLE_WORKERS"
      --run_clair3 1
      --run_freebayes 1
      --run_bcftools 1
      --clair3_model_path "$CLAIR3_MODEL_PATH"
      --clair3_platform "$CLAIR3_PLATFORM"
    )

    if [ -n "$NANOPORE_METADATA" ]; then
        CMD+=(--metadata "$NANOPORE_METADATA")
    fi

    if [ "$RUN_DEPLOID" -eq 1 ]; then
        CMD+=(--run_deploid)

        if [ -n "$POPULATION_VCF" ]; then
            CMD+=(--population_vcf "$POPULATION_VCF")
        fi

        if [ -n "$GENETIC_MAP" ]; then
            CMD+=(--genetic_map "$GENETIC_MAP")
        fi

        if [ -n "$EXCLUDE_VCF" ]; then
            CMD+=(--exclude_vcf "$EXCLUDE_VCF")
        fi
    fi

    "${CMD[@]}"

    ok "Nanopore pipeline completed"
}

############################################
# MAIN
############################################

main() {
    echo -e "${BLUE}====================================================${RESET}"
    echo -e "${BLUE}  PF HRP2/HRP3 Deletion Pipeline${RESET}"
    echo -e "${BLUE}  Platform: ${PLATFORM}${RESET}"
    echo -e "${BLUE}====================================================${RESET}"

    case "$PLATFORM" in
        illumina|ILLUMINA)
            run_illumina
            ;;
        nanopore|NANOPORE|ont|ONT)
            run_nanopore
            ;;
        *)
            fail "Unknown platform: $PLATFORM"
            echo "Usage:"
            echo "  bash main.sh illumina"
            echo "  bash main.sh nanopore"
            exit 1
            ;;
    esac

    echo -e "${GREEN}====================================================${RESET}"
    echo -e "${GREEN}  ALL DONE${RESET}"
    echo -e "${GREEN}====================================================${RESET}"
}

main "$@"