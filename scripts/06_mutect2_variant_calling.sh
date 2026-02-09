#!/usr/bin/env bash
set -e

source config/sample_config.sh

echo "[INFO] Starting somatic variant calling with Mutect2"
echo "[INFO] Sample ID: ${SAMPLE_ID}"
echo "[INFO] Analysis type: ${ANALYSIS_TYPE}"

mkdir -p "${VARIANT_DIR}"

DEDUP_BAM="${ALIGNMENT_DIR}/${SAMPLE_ID}.dedup.bam"
UNFILTERED_VCF="${VARIANT_DIR}/${SAMPLE_ID}.unfiltered.vcf.gz"

gatk Mutect2 \
  -R "${REFERENCE_FA}" \
  -I "${DEDUP_BAM}" \
  -tumor "${SAMPLE_ID}" \
  -L "${PANEL_BED}" \
  --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
  -O "${UNFILTERED_VCF}"

echo "[INFO] Mutect2 completed"
echo "[INFO] Output VCF: ${UNFILTERED_VCF}"
