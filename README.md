# Metagenomic Data Analysis Pipeline

This repository contains the steps and scripts for analyzing metagenomic data, from raw reads to functional and taxonomic annotations. Below is a detailed breakdown of the pipeline.

## Steps
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

### Homopolymer Removal
Remove homopolymers using Prinseq:
```bash
qsub 3_Prinseq.sh
prinseq-lite.pl -fastq ../sample_R1.fastq_trimmed -fastq2 ../sample_R2.fastq_trimmed -custom_params "AAT 10;T 70%;A 15;G 70%;C 15"
```

### Barcode Removal
Trim adapters and low-quality bases using Trimmomatic:
```bash
qsub 4_Trimmomatic.sh
java -jar $TRIMMOMATIC/trimmomatic.jar PE -phred33 input_R1.fastq input_R2.fastq output_R1_P.fq output_R1_U.fq output_R2_P.fq output_R2_U.fq ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq2-PE.fa:2:40:15 LEADING:10 TRAILING:10 SLIDINGWINDOW:10:20 MINLEN:80
```

### Host Contamination Removal
Remove host-derived reads using Bowtie2:
```bash
qsub 5_Host_contaRemoval.sh
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 input_R1.fastq -2 input_R2.fastq -S output.sam
samtools view -bS output.sam > output.bam
samtools view -b -f 12 -F 256 output.bam > output_unmapped.bam
samtools sort -n output_unmapped.bam > output_sorted.bam
bedtools bamtofastq -i output_sorted.bam -fq output_host_removed_r1.fastq -fq2 output_host_removed_r2.fastq
```

## 2. Metagenomic Assembly
### SPAdes Assembly
Assemble reads using SPAdes:
```bash
qsub 6_Spades_Assembly.sh
spades.py --meta -1 input_R1.fastq -2 input_R2.fastq -o assembly_output
```
### Contigs Filtering
Filter contigs based on length:
```bash
perl Filter_contigs_read_length.pl
```

## 3. Gene Prediction and Quantification
### Gene Prediction
Predict genes using Prodigal:
```bash
qsub 7_Prodigal.sh
prodigal -i contigs.fasta -o genes.gff -a proteins.faa -p meta
```
### Gene Quantification
Quantify genes using BWA and SAMtools:
```bash
qsub newBWA.sh
bwa index genes.fna
bwa mem genes.fna input_R1.fastq input_R2.fastq | samtools view -Sb - > output.bam
samtools sort output.bam -o sorted_output.bam
samtools index sorted_output.bam
```

## 4. Non-Redundant Gene-Set Construction
### CD-HIT
Create a non-redundant gene set:
```bash
qsub 9_cdhit.sh
cd-hit-est -i genes.fna -o genes_cdhit.fna -c 0.95 -n 10 -aS 0.9 -d 0 -T 48 -M 60000
```

## 5. Functional Annotation
### KEGG Analysis
Annotate genes using KEGG:
```bash
blastp -query genes.faa -db KEGG_DB (pre-built) -out kegg_out.txt -num_threads 24
```
### CAZy Analysis
Annotate genes using CAZyDB ((pre-built)):
```bash
blastp -query genes.faa -db CAZyDB -out cazy_out.txt -num_threads 24
```
### Antibiotic Resistance Analysis
Annotate genes using ARDB ((pre-built)):
```bash
blastp -query genes.faa -db ARDB -out ardb_out.txt -num_threads 24
```
### Xenobiotics Degradation Analysis
Annotate genes for xenobiotics degradation:
```bash
blastp -query genes.faa -db Xenobiotics_DB(pre-built) -out xenobiotics_out.txt -num_threads 24
```

## 6. Taxonomic Analysis
### Taxonomic Assignment
Assign taxonomy using BLAST against NCBI and HMP databases (pre-built):
```bash
blastn -query genes.fna -db NCBI_HMP_DB -out taxonomic_out.txt -num_threads 50
```

## 7. Genome Reconstruction
### Cross-Assembly
Assemble reads using Megahit:
```bash
qsub Megahit_human_assembly.sh
```
### Binning and Quality Assessment
Bin contigs using MetaBAT and assess quality using CheckM:
```bash
qsub binning_quality_check.sh
```

## 8. Statistical Analysis
### Perform statistical analyses on metagenomic data:
Assemble reads using Megahit:
```bash
Rscript statistical_analysis.R
```
This pipeline provides a comprehensive workflow for metagenomic data analysis, covering data preprocessing, assembly, annotation, taxonomic classification, and genome reconstruction. Please make sure all dependencies are correctly installed before running the scripts. 
- All intermediate scripts can be found in the **Final_Scripts** folder!!

### Contact
For any questions, please contact: 👉 Ashok K. Sharma; ashoks773@gmail.com
