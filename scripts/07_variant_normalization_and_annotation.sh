#!/usr/bin/env bash
set -e

source config/sample_config.sh

echo "[INFO] Starting variant normalization and annotation"
echo "[INFO] Sample ID: ${SAMPLE_ID}"

mkdir -p "${VARIANT_DIR}"

UNFILTERED_VCF="${VARIANT_DIR}/${SAMPLE_ID}.unfiltered.vcf.gz"
NORMALIZED_VCF="${VARIANT_DIR}/${SAMPLE_ID}.filtered.norm.vcf.gz"
ANNOTATED_VCF="${VARIANT_DIR}/${SAMPLE_ID}.filtered.norm.snpeff.vcf.gz"

# Normalize and left-align variants
bcftools norm \
  -f "${REFERENCE_FA}" \
  -m -both \
  "${UNFILTERED_VCF}" \
  -Oz -o "${NORMALIZED_VCF}"

bcftools index "${NORMALIZED_VCF}"

# Annotate variants with SnpEff
snpEff \
  -canon \
  -hgvs \
  "${GENOME_BUILD}" \
  "${NORMALIZED_VCF}" \
  | bgzip -c > "${ANNOTATED_VCF}"

bcftools index "${ANNOTATED_VCF}"

echo "[INFO] Variant normalization and annotation completed"
echo "[INFO] Annotated VCF: ${ANNOTATED_VCF}"
