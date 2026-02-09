#!/usr/bin/env bash
set -e

source config/sample_config.sh

echo "[INFO] Running FastQC on raw FASTQ files"
echo "[INFO] Sample ID: ${SAMPLE_ID}"

mkdir -p "${BASE_DIR}/fastqc_raw"

fastqc \
  "${R1_FASTQ}" \
  "${R2_FASTQ}" \
  -o "${BASE_DIR}/fastqc_raw"

echo "[INFO] FastQC completed successfully for ${SAMPLE_ID}"
