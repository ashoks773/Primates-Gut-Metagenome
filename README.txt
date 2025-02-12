# Metagenomic Data Analysis Pipeline

## Overview
This repository contains a pipeline for metagenomic data analysis, including preprocessing, assembly, gene prediction, and quantification.

## Filtering Steps

### 1. Quality Control
Run FastQC to assess read quality:
```
qsub 1_Fastqc_only.sh
```

### 2. Ambiguity Filtering
```sh
qsub 2_Ambiguity.sh
perl ~/Softwares/NGSQCToolkit_v2.3.3/Trimming/AmbiguityFiltering.pl \
    -i ../HU100_S17_R1_001.fastq -irev ../HU100_S17_R2_001.fastq -c 1 -t5 -t3
```

### 3. Homopolymer Removal
```sh
qsub 3_Prinseq.sh
prinseq-lite.pl -fastq ../HU100_S17_R1_001.fastq_trimmed \
    -fastq2 ../HU100_S17_R2_001.fastq_trimmed -custom_params "AAT 10;T 70%;A 15;G 70%;C 15"
```

### 4. Barcode Removal
```sh
qsub 4_Trimmomatic.sh
java -jar $TRIMMOMATIC/trimmomatic.jar PE -phred33 \
    ../HU100_S17_R1_001_prinseq_good_OLyg.fastq \
    ../HU100_S17_R2_001_prinseq_good_Rtjq.fastq \
    HU100_S17_R1_001_prinseq_good_OLyg.fastq.PE \
    HU100_S17_R1_001_prinseq_good_OLyg.fastq.UN \
    HU100_S17_R2_001_prinseq_good_Rtjq.fastq.PE \
    HU100_S17_R2_001_prinseq_good_Rtjq.fastq.UN \
    ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq2-PE.fa:2:40:15 \
    LEADING:10 TRAILING:10 SLIDINGWINDOW:10:20 MINLEN:80
```

### 5. Host Contamination Removal
```sh
qsub 5_Host_contaRemoval.sh
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex \
    -1 ../HU100_S17_R1_001_prinseq_good_OLyg.fastq.PE \
    -2 ../HU100_S17_R2_001_prinseq_good_Rtjq.fastq.PE \
    -S HU100_S17_prinseq_good_Rtjq_mapped_and_unmapped.sam

samtools view -bS HU100_S17_prinseq_good_Rtjq_mapped_and_unmapped.sam > HU100_S17_prinseq_good_Rtjq_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU100_S17_prinseq_good_Rtjq_mapped_and_unmapped.bam > HU100_S17_prinseq_good_Rtjq_bothEnds_Unmapped.bam
samtools sort -n HU100_S17_prinseq_good_Rtjq_bothEnds_Unmapped.bam > HU100_S17_prinseq_good_Rtjq_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU100_S17_prinseq_good_Rtjq_bothEnds_Unmapped_sorted.bam \
    -fq HU100_S17_prinseq_good_Rtjq_host_removed_r1.fastq \
    -fq2 HU100_S17_prinseq_good_Rtjq_host_removed_r2.fastq
```

## Metagenomic Assembly
### 6. SPAdes Assembly
```sh
qsub 6_Spades_Assembly.sh
/home/gomeza/sharm646/Softwares/SPAdes-3.12.0-Linux/bin/spades.py --meta \
    -1 ../HU100_S17_prinseq_good_Rtjq_host_removed_r1.fastq \
    -2 ../HU100_S17_prinseq_good_Rtjq_host_removed_r2.fastq \
    -s ../HU100_S17_R2_001_prinseq_good_Rtjq_Nomates.fastq \
    -m 950 -t 30 -o HU100_S17_prinseq_good_Rtjq_metaspades_out_result
```

### 7. Contigs Filtering
```sh
perl Filter_contigs_read_length.pl
```

## Gene Prediction
```sh
qsub 7_Prodigal.sh
prodigal -i HU100_S17_scaffolds.fasta.Filtered.fna \
    -o HU100_S17.Genes -a HU100_S17.proteins.faa -p meta
```

## Gene Catalogue Construction
### 8. Combining Human and Gorilla Data
```sh
perl 8_index.pl *fna
cat *indexed_GENES.fna > Combined_Human_Indexed.fna
perl 8_index.pl *fna
cat *indexed_GENES.fna > Combined_Gorilla_Indexed.fna
cat Combined_Human_Indexed.fna Combined_Gorilla_Indexed.fna > Combined_Total_Genes.fna
```

### 9. Construction of Non-Redundant Gene Set
```sh
qsub 9_cdhit.sh
cd-hit-est -i ../Combined_Total_Genes.fna -o Combined_Total_Genes.cdhit.fna \
    -c 0.95 -n 10 -aS 0.9 -d 0 -T 48 -M 60000
```

## Gene Quantification
```sh
bwa index Combined_Total_Genes.cdhit.fna
bwa mem -t 144 -M Combined_Total_Genes.cdhit.fna \
    /home/gomeza/sharm646/Metagenomic/Human/N_Trimmed/Homopolymer/Trimmomatic/Host_removal/HU100_S17_prinseq_good_Rtjq_host_removed_r1.fastq \
    /home/gomeza/sharm646/Metagenomic/Human/N_Trimmed/Homopolymer/Trimmomatic/Host_removal/HU100_S17_prinseq_good_Rtjq_host_removed_r2.fastq | samtools view -Sb - > HU100_S17.bam
```

### 10. Generate Gene Count Matrix
```sh
samtools idxstats 2009_106.sort.bam | grep -v "*" | cut -f1 > gene_names
samtools idxstats 2009_106.sort.bam | grep -v "*" | cut -f3 > Counts1
samtools idxstats 2009_107.sort.bam | grep -v "*" | cut -f3 > Counts2
paste gene_names Counts1 Counts2 ... > Gene_count_matrix.tab
```

### 11. Gene Abundance Calculation
```sh
perl seq_single_line.pl Combined_Total_Genes.cdhit.fna
mv Combined_Total_Genes.cdhit.fna.single Combined_Total_Genes.cdhit.single.fna
perl Calculate_gene_abundance.pl Combined_Total_Genes.cdhit.single.fna Gene_count_matrix.tab
perl Filter_genes.pl GENE_ABUNDANCE_TABLE
```

## Normalization and Analysis
```r
Rscript Gene_proportions.R
Gene_abundance <- read.csv("GENE_ABUNDANCE_TABLE.filtered", row.names = 1, header = TRUE, sep = "\t")
Gene_abundance_normalized <- Gene_abundance / colSums(Gene_abundance)
```

## Contact
For any questions, please reach out to Ashok Kumar Sharma at ashoks773@gmail.com
