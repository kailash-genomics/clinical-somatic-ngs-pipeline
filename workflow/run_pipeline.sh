#!/usr/bin/env bash
set -e

echo "=========================================="
echo " Clinical Somatic NGS Pipeline"
echo "=========================================="

source config/sample_config.sh

echo "[INFO] Sample ID      : ${SAMPLE_ID}"
echo "[INFO] Assay          : ${ASSAY_TYPE}"
echo "[INFO] Genome build   : ${GENOME_BUILD}"
echo "[INFO] Analysis type : ${ANALYSIS_TYPE}"
echo "------------------------------------------"

bash scripts/00_setup_colab.sh
bash scripts/01_fastqc.sh
bash scripts/02_fastp_trimming.sh
bash scripts/03_alignment_bwa.sh
bash scripts/04_read_groups_markdup.sh
bash scripts/05_alignment_qc_and_coverage.sh
bash scripts/06_mutect2_variant_calling.sh
bash scripts/07_variant_normalization_and_annotation.sh
bash scripts/08_clinical_variant_filtering.sh

echo "=========================================="
echo "[INFO] Pipeline completed successfully"
echo "=========================================="
