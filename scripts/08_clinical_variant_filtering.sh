#!/usr/bin/env bash
set -e

source config/sample_config.sh

echo "[INFO] Starting clinical variant filtering"
echo "[INFO] Sample ID: ${SAMPLE_ID}"
echo "[INFO] Minimum depth: ${MIN_DP}"
echo "[INFO] Minimum AF: ${MIN_AF}"

mkdir -p "${VARIANT_DIR}"

ANNOTATED_VCF="${VARIANT_DIR}/${SAMPLE_ID}.filtered.norm.snpeff.vcf.gz"
CLINICAL_VCF="${VARIANT_DIR}/${SAMPLE_ID}.clinical.pass.vcf.gz"

bcftools view \
  -i "FILTER=\"PASS\" && INFO/DP>=${MIN_DP} && FORMAT/AF>=${MIN_AF}" \
  -R "${PANEL_BED}" \
  "${ANNOTATED_VCF}" \
  -Oz -o "${CLINICAL_VCF}"

bcftools index "${CLINICAL_VCF}"

echo "[INFO] Clinical variant filtering completed"
echo "[INFO] Final clinical VCF: ${CLINICAL_VCF}"
