# Clinical Somatic NGS Pipeline  
Tumor-Only Targeted Sequencing (92-Gene Panel)

## Overview

This repository contains a production-style, end-to-end clinical somatic NGS analysis pipeline designed for tumor-only targeted DNA sequencing panels used in solid tumor profiling.

The pipeline reflects clinical laboratory best practices and demonstrates:

- End-to-end ownership of NGS data analysis  
- SOP-driven, modular pipeline design  
- Robust quality control, filtering, and reporting logic  
- Reproducibility, traceability, and data governance awareness  

This project is intended as:

- A professional bioinformatics portfolio  
- An internal validation–style reference pipeline  
- A reference implementation for targeted somatic analysis workflows  

---

## Analysis Scope

Analysis type: Tumor-only somatic variant detection  
Assay: Targeted solid tumor panel (92 genes)  
Genome build: hg19 / GRCh37  
Variant caller: GATK Mutect2 (tumor-only mode)  
Reporting focus: High-confidence, clinically relevant variants  

---

## Pipeline Workflow (End-to-End)

1. Raw FASTQ quality control (FastQC)  
2. Adapter trimming and quality filtering (fastp)  
3. Alignment to reference genome (BWA-MEM)  
4. BAM sorting and indexing (samtools)  
5. Read group assignment and duplicate marking (GATK)  
6. Alignment QC and on-target coverage analysis  
7. Somatic variant calling (Mutect2, tumor-only mode)  
8. Variant normalization and left-alignment (bcftools)  
9. Functional annotation (SnpEff)  
10. Clinical variant filtering  
   - PASS variants  
   - Depth (DP) threshold  
   - Allele fraction (AF) threshold  
   - Target panel BED restriction  
11. Final Excel-ready variant report generation (TSV)  

---

## Toolchain

FastQC – Raw read quality assessment  
fastp – Adapter trimming and read filtering  
BWA-MEM – Alignment to reference genome  
samtools – BAM processing and alignment QC  
GATK 4.x – Read groups, duplicate marking, Mutect2  
bcftools – Variant normalization and filtering  
SnpEff – Functional annotation  

Exact tool versions used are documented within the pipeline scripts.

---

## Repository Structure

clinical-somatic-ngs-pipeline/  
├── workflow/           One-command pipeline runner  
├── scripts/            Modular pipeline steps (00–08)  
├── config/             Centralized configuration  
├── docs/               Supporting documentation  
├── example_results/    Safe example outputs  
├── validation_notes/   Interpretation and validation notes  
├── .gitignore          Data governance rules  
└── README.md  

---

## Configuration-Driven Design

All sample-specific and environment-specific parameters are centralized in:

config/sample_config.sh

This includes:

- Sample ID  
- Input FASTQ paths  
- Reference genome path  
- Panel BED file  
- Reporting thresholds (DP, AF)  

Only this file needs to be modified to re-run the pipeline on a new sample.

---

## How to Run

Google Colab or Linux environment:

bash workflow/run_pipeline.sh

The pipeline is fully modular and can also be executed step-by-step using individual scripts if required.

---

## Example Outputs

To demonstrate real execution without exposing patient data, the repository includes:

- A sample final variant report (TSV)  
- A summary of coverage metrics  

Located in:

example_results/

---

## Data Governance and Ethics

- No raw sequencing data is included  
- No BAM, VCF, or reference genome files are tracked  
- No patient-identifiable information is present  
- .gitignore explicitly prevents accidental data leakage  

This repository is safe for public sharing and professional review.

---

## Notes on Clinical Interpretation

Variant interpretation (pathogenicity classification, therapeutic relevance, and clinical actionability) is a manual expert process and is not automated in this pipeline.

Typical external resources include:

ClinVar  
COSMIC  
OncoKB  
NCCN / ESMO guidelines  

---

## Disclaimer

This pipeline is provided for demonstration and educational purposes only.  
Clinical deployment requires laboratory validation, quality management approval, and regulatory compliance.

---

## Author

Kailash  
Clinical Laboratory Scientist | Biotechnology  
NGS and Molecular Diagnostics
