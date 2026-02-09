#!/usr/bin/env bash
set -e

echo "[INFO] Running fastp for adapter trimming and quality filtering"

BASE_DIR="/content/drive/MyDrive/Somatic_Demo_Analysis"
SAMPLE_ID="DEMO_TUMOR_001"

# Input FASTQ files
R1_IN="${BASE_DIR}/raw_fastq/DEMO_TUMOR_001_R1.fastq.gz"
R2_IN="${BASE_DIR}/raw_fastq/DEMO_TUMOR_001_R2.fastq.gz"

# Output FASTQ files
R1_OUT="${BASE_DIR}/trimmed_fastq/DEMO_TUMOR_001_trimmed_R1.fastq.gz"
R2_OUT="${BASE_DIR}/trimmed_fastq/DEMO_TUMOR_001_trimmed_R2.fastq.gz"

# fastp QC outputs
HTML_OUT="${BASE_DIR}/fastqc_trimmed/DEMO_TUMOR_001_fastp.html"
JSON_OUT="${BASE_DIR}/fastqc_trimmed/DEMO_TUMOR_001_fastp.json"

echo "[INFO] Sample ID: ${SAMPLE_ID}"
echo "[INFO] Input R1: ${R1_IN}"
echo "[INFO] Input R2: ${R2_IN}"

# Sanity check
if [[ ! -f "${R1_IN}" || ! -f "${R2_IN}" ]]; then
  echo "[ERROR] Input FASTQ files not found"
  exit 1
fi

mkdir -p "${BASE_DIR}/trimmed_fastq"
mkdir -p "${BASE_DIR}/fastqc_trimmed"

# Run fastp
fastp \
  -i "${R1_IN}" \
  -I "${R2_IN}" \
  -o "${R1_OUT}" \
  -O "${R2_OUT}" \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --html "${HTML_OUT}" \
  --json "${JSON_OUT}"

echo "[INFO] fastp trimming completed successfully for ${SAMPLE_ID}"
