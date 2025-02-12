# Metagenomic Data Analysis Pipeline

This repository contains the steps and scripts for analyzing metagenomic data, from raw reads to functional and taxonomic annotations. Below is a detailed breakdown of the pipeline.

## Table of Contents
1. Filtering Steps  
2. Metagenomic Assembly  
3. Gene Prediction and Quantification  
4. Non-Redundant Gene-Set Construction  
5. Functional Annotation  
   - KEGG Analysis  
   - CAZy Analysis  
   - Antibiotic Resistance Analysis  
   - Xenobiotics Degradation Analysis  
6. Taxonomic Analysis  
7. Genome Reconstruction  
8. Statistical Analysis  

## 1. Filtering Steps
### FastQC
Run FastQC to assess the quality of raw reads:
```bash
qsub 1_Fastqc_only.sh
```
### Ambiguity Filtering
Filter out reads with ambiguous bases:
```bash
qsub 2_Ambiguity.sh
perl ~/Softwares/NGSQCToolkit_v2.3.3/Trimming/AmbiguityFiltering.pl -i ../sample_R1.fastq -irev ../sample_R2.fastq -c 1 -t5 -t3
```
