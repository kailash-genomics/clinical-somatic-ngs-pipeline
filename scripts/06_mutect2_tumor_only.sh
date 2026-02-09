#!/usr/bin/env bash
set -e

echo "[INFO] Starting Mutect2 tumor-only variant calling"

BASE_DIR="/content/drive/MyDrive/Somatic_Demo_Analysis"
SAMPLE_ID="DEMO_TUMOR_001"

# Inputs
BAM="${BASE_DIR}/alignment/DEMO_TUMOR_001.dedup.bam"
REF="${BASE_DIR}/reference/hg19.p13.plusMT.no_alt_analysis_set/hg19.p13.plusMT.no_alt_analysis_set.fa"
BED="${BASE_DIR}/bed/FINAL_92gene_panel.bed"

# Outputs
OUT_DIR="${BASE_DIR}/variants"
UNFILTERED_VCF="${OUT_DIR}/DEMO_TUMOR_001.unfiltered.vcf.gz"
FILTERED_VCF="${OUT_DIR}/DEMO_TUMOR_001.filtered.vcf.gz"

mkdir -p "$OUT_DIR"

# Sanity checks
for f in "$BAM" "$REF" "$BED"; do
  if [[ ! -f "$f" ]]; then
    echo "[ERROR] Required file not found: $f"
    exit 1
  fi
done

#####################################
# Mutect2 (tumor-only)
#####################################
gatk Mutect2 \
  -R "$REF" \
  -I "$BAM" \
  -tumor "$SAMPLE_ID" \
  -L "$BED" \
  --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
  -O "$UNFILTERED_VCF"

#####################################
# Filter Mutect2 calls
#####################################
gatk FilterMutectCalls \
  -R "$REF" \
  -V "$UNFILTERED_VCF" \
  -O "$FILTERED_VCF"

echo "[INFO] Mutect2 calling and filtering completed successfully"
