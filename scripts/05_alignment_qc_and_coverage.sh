#!/usr/bin/env bash
set -e

source config/sample_config.sh

echo "[INFO] Running alignment QC and targeted coverage analysis"
echo "[INFO] Sample ID: ${SAMPLE_ID}"

mkdir -p "${METRICS_DIR}"
mkdir -p "${LOG_DIR}"

DEDUP_BAM="${ALIGNMENT_DIR}/${SAMPLE_ID}.dedup.bam"

# Basic BAM integrity check
samtools quickcheck "${DEDUP_BAM}"

# Alignment QC
samtools flagstat "${DEDUP_BAM}" > "${METRICS_DIR}/${SAMPLE_ID}.flagstat.txt"
samtools idxstats "${DEDUP_BAM}" > "${METRICS_DIR}/${SAMPLE_ID}.idxstats.txt"

# Targeted coverage using BED
samtools depth \
  -b "${PANEL_BED}" \
  "${DEDUP_BAM}" \
  > "${METRICS_DIR}/${SAMPLE_ID}.on_target_coverage.txt"

# Mean on-target depth
awk '{sum+=$3; count++} END {if(count>0) print "Mean_on_target_depth\t"sum/count; else print "Mean_on_target_depth\t0"}' \
  "${METRICS_DIR}/${SAMPLE_ID}.on_target_coverage.txt" \
  > "${METRICS_DIR}/${SAMPLE_ID}.mean_on_target_depth.txt"

echo "[INFO] Alignment QC and coverage completed"
