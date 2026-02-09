# Clinical Somatic NGS Pipeline (Tumor-Only, Targeted Panel)

## Overview
This repository contains a production-style, end-to-end clinical somatic NGS analysis pipeline
designed for tumor-only targeted DNA sequencing panels (solid tumors).

The pipeline follows clinical laboratory best practices and demonstrates:
- End-to-end ownership of NGS data analysis
- SOP-driven pipeline design
- Proper QC, filtering, and reporting logic
- Reproducibility and data governance awareness

This project is intended as:
- A professional portfolio
- An internal validation-style demonstration
- A reference implementation for targeted somatic analysis

---

## Analysis Scope
- Analysis type: Tumor-only somatic variant detection  
- Assay: Targeted solid tumor panel (92 genes)  
- Genome build: hg19 / GRCh37  
- Variant caller: GATK Mutect2  
- Reporting focus: High-confidence, clinically relevant variants  

---

## Pipeline Steps (End-to-End)

1. Raw FASTQ quality control (FastQC)
2. Adapter trimming and quality filtering (fastp)
3. Alignment to reference genome (BWA-MEM)
4. BAM sorting and indexing
5. Read group assignment and duplicate marking (GATK)
6. Alignment QC and on-target coverage analysis
7. Somatic variant calling (Mutect2, tumor-only mode)
8. Variant normalization and left-alignment (bcftools)
9. Functional annotation (SnpEff)
10. Clinical variant filtering (PASS, depth, allele fraction, panel BED)
11. Final Excel-ready variant report generation

---

## Toolchain

| Tool | Purpose |
|------|--------|
| FastQC | Raw read quality assessment |
| fastp | Adapter trimming and read filtering |
| BWA-MEM | Alignment to reference genome |
| samtools | BAM processing and QC |
| GATK 4.x | Read groups, duplicate marking, Mutect2 |
| bcftools | Variant normalization and filtering |
| SnpEff | Functional annotation |

Exact versions used are documented in the pipeline scripts.

---

## Repository Structure

clinical-somatic-ngs-pipeline/
- workflow/           One-command pipeline runner  
- scripts/            Modular analysis steps  
- config/             Centralized configuration  
- docs/               Supporting documentation  
- example_results/    Safe example outputs  
- validation_notes/   Interpretation notes  
- .gitignore          Data governance rules  

---

## Configuration-Driven Design

All sample-specific and environment-specific parameters are defined in:

config/sample_config.sh

This includes:
- Sample ID
- Input FASTQ paths
- Reference genome
- Panel BED
- Reporting thresholds (DP, AF)

Changing one file is sufficient to re-run the pipeline on a new sample.

---

## How to Run

Google Colab / Linux environment:

bash workflow/run_pipeline.sh

The pipeline is modular and can also be executed step-by-step if required.

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
- All patient-identifiable data is excluded
- .gitignore explicitly prevents accidental data leakage

---

## Notes on Clinical Interpretation

Variant interpretation (pathogenicity classification, therapeutic relevance)
is a manual expert process and is not automated in this pipeline.

Typical external resources include:
- ClinVar
- COSMIC
- OncoKB
- NCCN / ESMO guidelines

---

## Disclaimer
This pipeline is provided for demonstration and educational purposes.
Clinical deployment requires laboratory validation and regulatory approval.

---

## Author
Kailash  
Clinical Laboratory Scientist | Biotechnology | NGS & Molecular Diagnostics
