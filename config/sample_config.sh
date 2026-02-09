#!/usr/bin/env bash

# ===============================
# SAMPLE METADATA
# ===============================
SAMPLE_ID="DEMO_TUMOR_001"
ASSAY_TYPE="Targeted solid tumor panel (92 genes)"
GENOME_BUILD="hg19"
ANALYSIS_TYPE="Tumor-only somatic"

# ===============================
# BASE DIRECTORIES
# ===============================
BASE_DIR="/content/drive/MyDrive/Somatic_Demo_Analysis"

RAW_FASTQ_DIR="${BASE_DIR}/raw_fastq"
TRIMMED_FASTQ_DIR="${BASE_DIR}/trimmed_fastq"
ALIGNMENT_DIR="${BASE_DIR}/alignment"
VARIANT_DIR="${BASE_DIR}/variants"
METRICS_DIR="${BASE_DIR}/metrics"
LOG_DIR="${BASE_DIR}/logs"
REPORT_DIR="${BASE_DIR}/reports"

# ===============================
# INPUT FILES
# ===============================
R1_FASTQ="${RAW_FASTQ_DIR}/DEMO_TUMOR_001_R1.fastq.gz"
R2_FASTQ="${RAW_FASTQ_DIR}/DEMO_TUMOR_001_R2.fastq.gz"

# ===============================
# REFERENCE & PANEL
# ===============================
REFERENCE_FA="${BASE_DIR}/reference/hg19.p13.plusMT.no_alt_analysis_set/hg19.p13.plusMT.no_alt_analysis_set.fa"
PANEL_BED="${BASE_DIR}/bed/FINAL_92gene_panel.bed"

# ===============================
# THREADING
# ===============================
THREADS=4

# ===============================
# CLINICAL FILTERS
# ===============================
MIN_DP=100
MIN_AF=0.05

