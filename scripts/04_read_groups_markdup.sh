#!/usr/bin/env bash
set -e

source config/sample_config.sh

echo "[INFO] Adding Read Groups and marking duplicates"
echo "[INFO] Sample ID: ${SAMPLE_ID}"

mkdir -p "${ALIGNMENT_DIR}"
mkdir -p "${METRICS_DIR}"

# Input BAM from alignment step
INPUT_BAM="${ALIGNMENT_DIR}/${SAMPLE_ID}.sorted.bam"

# Output files
RG_BAM="${ALIGNMENT_DIR}/${SAMPLE_ID}.RG.sorted.bam"
DEDUP_BAM="${ALIGNMENT_DIR}/${SAMPLE_ID}.dedup.bam"
METRICS_FILE="${METRICS_DIR}/${SAMPLE_ID}.duplication_metrics.txt"

# Add or replace read groups (GATK best practice)
gatk AddOrReplaceReadGroups \
  -I "${INPUT_BAM}" \
  -O "${RG_BAM}" \
  -RGID "${SAMPLE_ID}" \
  -RGLB "TARGETED_PANEL" \
  -RGPL "ILLUMINA" \
  -RGPU "UNIT1" \
  -RGSM "${SAMPLE_ID}"

# Mark duplicates
gatk MarkDuplicates \
  -I "${RG_BAM}" \
  -O "${DEDUP_BAM}" \
  -M "${METRICS_FILE}" \
  --CREATE_INDEX true

echo "[INFO] Read groups added and duplicates marked"
echo "[INFO] Final BAM: ${DEDUP_BAM}"
