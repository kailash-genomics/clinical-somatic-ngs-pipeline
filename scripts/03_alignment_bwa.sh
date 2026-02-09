#!/usr/bin/env bash
set -e

source config/sample_config.sh

echo "[INFO] Starting BWA-MEM alignment"
echo "[INFO] Sample ID: ${SAMPLE_ID}"
echo "[INFO] Reference: ${REFERENCE_FA}"

mkdir -p "${ALIGNMENT_DIR}"

# Define outputs
SORTED_BAM="${ALIGNMENT_DIR}/${SAMPLE_ID}.sorted.bam"

# Align and sort directly (no intermediate SAM)
bwa mem -t "${THREADS}" \
  "${REFERENCE_FA}" \
  "${TRIMMED_R1}" \
  "${TRIMMED_R2}" \
| samtools sort -@ "${THREADS}" -o "${SORTED_BAM}"

# Index sorted BAM
samtools index "${SORTED_BAM}"

echo "[INFO] Alignment completed: ${SORTED_BAM}"
