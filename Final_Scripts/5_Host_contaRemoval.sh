#!/bin/bash -l
#PBS -l walltime=48:00:00,nodes=6:ppn=15,mem=60gb
#PBS -m abe
#PBS -e Hostremoval.error
#PBS -o Hostremoval.out

module load bowtie2
module load samtools
module load bedtools

cd ~/Metagenomic/Human/N_Trimmed/Homopolymer/Trimmomatic/Host_removal

bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU100_S17_R1_001_prinseq_good_OLyg.fastq.PE -2 ../HU100_S17_R2_001_prinseq_good_Rtjq.fastq.PE -S HU100_S17_prinseq_good_Rtjq_mapped_and_unmapped.sam
samtools view -bS HU100_S17_prinseq_good_Rtjq_mapped_and_unmapped.sam > HU100_S17_prinseq_good_Rtjq_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU100_S17_prinseq_good_Rtjq_mapped_and_unmapped.bam > HU100_S17_prinseq_good_Rtjq_bothEnds_Unmapped.bam
samtools sort -n HU100_S17_prinseq_good_Rtjq_bothEnds_Unmapped.bam > HU100_S17_prinseq_good_Rtjq_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU100_S17_prinseq_good_Rtjq_bothEnds_Unmapped_sorted.bam -fq HU100_S17_prinseq_good_Rtjq_host_removed_r1.fastq -fq2 HU100_S17_prinseq_good_Rtjq_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU101_S15_R1_001_prinseq_good_JAKY.fastq.PE -2 ../HU101_S15_R2_001_prinseq_good_I9fo.fastq.PE -S HU101_S15_prinseq_good_I9fo_mapped_and_unmapped.sam
samtools view -bS HU101_S15_prinseq_good_I9fo_mapped_and_unmapped.sam > HU101_S15_prinseq_good_I9fo_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU101_S15_prinseq_good_I9fo_mapped_and_unmapped.bam > HU101_S15_prinseq_good_I9fo_bothEnds_Unmapped.bam
samtools sort -n HU101_S15_prinseq_good_I9fo_bothEnds_Unmapped.bam > HU101_S15_prinseq_good_I9fo_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU101_S15_prinseq_good_I9fo_bothEnds_Unmapped_sorted.bam -fq HU101_S15_prinseq_good_I9fo_host_removed_r1.fastq -fq2 HU101_S15_prinseq_good_I9fo_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU103_S22_R1_001_prinseq_good_wZ5w.fastq.PE -2 ../HU103_S22_R2_001_prinseq_good_UH26.fastq.PE -S HU103_S22_prinseq_good_UH26_mapped_and_unmapped.sam
samtools view -bS HU103_S22_prinseq_good_UH26_mapped_and_unmapped.sam > HU103_S22_prinseq_good_UH26_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU103_S22_prinseq_good_UH26_mapped_and_unmapped.bam > HU103_S22_prinseq_good_UH26_bothEnds_Unmapped.bam
samtools sort -n HU103_S22_prinseq_good_UH26_bothEnds_Unmapped.bam > HU103_S22_prinseq_good_UH26_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU103_S22_prinseq_good_UH26_bothEnds_Unmapped_sorted.bam -fq HU103_S22_prinseq_good_UH26_host_removed_r1.fastq -fq2 HU103_S22_prinseq_good_UH26_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU107_S1_R1_001_prinseq_good_fkQ1.fastq.PE -2 ../HU107_S1_R2_001_prinseq_good_7_eb.fastq.PE -S HU107_S1_prinseq_good_7_eb_mapped_and_unmapped.sam
samtools view -bS HU107_S1_prinseq_good_7_eb_mapped_and_unmapped.sam > HU107_S1_prinseq_good_7_eb_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU107_S1_prinseq_good_7_eb_mapped_and_unmapped.bam > HU107_S1_prinseq_good_7_eb_bothEnds_Unmapped.bam
samtools sort -n HU107_S1_prinseq_good_7_eb_bothEnds_Unmapped.bam > HU107_S1_prinseq_good_7_eb_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU107_S1_prinseq_good_7_eb_bothEnds_Unmapped_sorted.bam -fq HU107_S1_prinseq_good_7_eb_host_removed_r1.fastq -fq2 HU107_S1_prinseq_good_7_eb_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU108_S18_R1_001_prinseq_good_QlzH.fastq.PE -2 ../HU108_S18_R2_001_prinseq_good_BAHK.fastq.PE -S HU108_S18_prinseq_good_BAHK_mapped_and_unmapped.sam
samtools view -bS HU108_S18_prinseq_good_BAHK_mapped_and_unmapped.sam > HU108_S18_prinseq_good_BAHK_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU108_S18_prinseq_good_BAHK_mapped_and_unmapped.bam > HU108_S18_prinseq_good_BAHK_bothEnds_Unmapped.bam
samtools sort -n HU108_S18_prinseq_good_BAHK_bothEnds_Unmapped.bam > HU108_S18_prinseq_good_BAHK_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU108_S18_prinseq_good_BAHK_bothEnds_Unmapped_sorted.bam -fq HU108_S18_prinseq_good_BAHK_host_removed_r1.fastq -fq2 HU108_S18_prinseq_good_BAHK_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU109_S19_R1_001_prinseq_good_BxAO.fastq.PE -2 ../HU109_S19_R2_001_prinseq_good_sWzB.fastq.PE -S HU109_S19_prinseq_good_sWzB_mapped_and_unmapped.sam
samtools view -bS HU109_S19_prinseq_good_sWzB_mapped_and_unmapped.sam > HU109_S19_prinseq_good_sWzB_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU109_S19_prinseq_good_sWzB_mapped_and_unmapped.bam > HU109_S19_prinseq_good_sWzB_bothEnds_Unmapped.bam
samtools sort -n HU109_S19_prinseq_good_sWzB_bothEnds_Unmapped.bam > HU109_S19_prinseq_good_sWzB_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU109_S19_prinseq_good_sWzB_bothEnds_Unmapped_sorted.bam -fq HU109_S19_prinseq_good_sWzB_host_removed_r1.fastq -fq2 HU109_S19_prinseq_good_sWzB_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU112_S2_R1_001_prinseq_good_OBtH.fastq.PE -2 ../HU112_S2_R2_001_prinseq_good_YZ7R.fastq.PE -S HU112_S2_prinseq_good_YZ7R_mapped_and_unmapped.sam
samtools view -bS HU112_S2_prinseq_good_YZ7R_mapped_and_unmapped.sam > HU112_S2_prinseq_good_YZ7R_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU112_S2_prinseq_good_YZ7R_mapped_and_unmapped.bam > HU112_S2_prinseq_good_YZ7R_bothEnds_Unmapped.bam
samtools sort -n HU112_S2_prinseq_good_YZ7R_bothEnds_Unmapped.bam > HU112_S2_prinseq_good_YZ7R_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU112_S2_prinseq_good_YZ7R_bothEnds_Unmapped_sorted.bam -fq HU112_S2_prinseq_good_YZ7R_host_removed_r1.fastq -fq2 HU112_S2_prinseq_good_YZ7R_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU113_S3_R1_001_prinseq_good_SsyW.fastq.PE -2 ../HU113_S3_R2_001_prinseq_good_exsy.fastq.PE -S HU113_S3_prinseq_good_exsy_mapped_and_unmapped.sam
samtools view -bS HU113_S3_prinseq_good_exsy_mapped_and_unmapped.sam > HU113_S3_prinseq_good_exsy_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU113_S3_prinseq_good_exsy_mapped_and_unmapped.bam > HU113_S3_prinseq_good_exsy_bothEnds_Unmapped.bam
samtools sort -n HU113_S3_prinseq_good_exsy_bothEnds_Unmapped.bam > HU113_S3_prinseq_good_exsy_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU113_S3_prinseq_good_exsy_bothEnds_Unmapped_sorted.bam -fq HU113_S3_prinseq_good_exsy_host_removed_r1.fastq -fq2 HU113_S3_prinseq_good_exsy_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU114_S27_R1_001_prinseq_good_7ooY.fastq.PE -2 ../HU114_S27_R2_001_prinseq_good___JQ.fastq.PE -S HU114_S27_prinseq_good___JQ_mapped_and_unmapped.sam
samtools view -bS HU114_S27_prinseq_good___JQ_mapped_and_unmapped.sam > HU114_S27_prinseq_good___JQ_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU114_S27_prinseq_good___JQ_mapped_and_unmapped.bam > HU114_S27_prinseq_good___JQ_bothEnds_Unmapped.bam
samtools sort -n HU114_S27_prinseq_good___JQ_bothEnds_Unmapped.bam > HU114_S27_prinseq_good___JQ_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU114_S27_prinseq_good___JQ_bothEnds_Unmapped_sorted.bam -fq HU114_S27_prinseq_good___JQ_host_removed_r1.fastq -fq2 HU114_S27_prinseq_good___JQ_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU115_S28_R1_001_prinseq_good_p8Ck.fastq.PE -2 ../HU115_S28_R2_001_prinseq_good_c8rI.fastq.PE -S HU115_S28_prinseq_good_c8rI_mapped_and_unmapped.sam
samtools view -bS HU115_S28_prinseq_good_c8rI_mapped_and_unmapped.sam > HU115_S28_prinseq_good_c8rI_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU115_S28_prinseq_good_c8rI_mapped_and_unmapped.bam > HU115_S28_prinseq_good_c8rI_bothEnds_Unmapped.bam
samtools sort -n HU115_S28_prinseq_good_c8rI_bothEnds_Unmapped.bam > HU115_S28_prinseq_good_c8rI_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU115_S28_prinseq_good_c8rI_bothEnds_Unmapped_sorted.bam -fq HU115_S28_prinseq_good_c8rI_host_removed_r1.fastq -fq2 HU115_S28_prinseq_good_c8rI_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU117_S4_R1_001_prinseq_good_uxvZ.fastq.PE -2 ../HU117_S4_R2_001_prinseq_good_5grf.fastq.PE -S HU117_S4_prinseq_good_5grf_mapped_and_unmapped.sam
samtools view -bS HU117_S4_prinseq_good_5grf_mapped_and_unmapped.sam > HU117_S4_prinseq_good_5grf_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU117_S4_prinseq_good_5grf_mapped_and_unmapped.bam > HU117_S4_prinseq_good_5grf_bothEnds_Unmapped.bam
samtools sort -n HU117_S4_prinseq_good_5grf_bothEnds_Unmapped.bam > HU117_S4_prinseq_good_5grf_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU117_S4_prinseq_good_5grf_bothEnds_Unmapped_sorted.bam -fq HU117_S4_prinseq_good_5grf_host_removed_r1.fastq -fq2 HU117_S4_prinseq_good_5grf_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU118_S5_R1_001_prinseq_good_RXre.fastq.PE -2 ../HU118_S5_R2_001_prinseq_good_MBZg.fastq.PE -S HU118_S5_prinseq_good_MBZg_mapped_and_unmapped.sam
samtools view -bS HU118_S5_prinseq_good_MBZg_mapped_and_unmapped.sam > HU118_S5_prinseq_good_MBZg_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU118_S5_prinseq_good_MBZg_mapped_and_unmapped.bam > HU118_S5_prinseq_good_MBZg_bothEnds_Unmapped.bam
samtools sort -n HU118_S5_prinseq_good_MBZg_bothEnds_Unmapped.bam > HU118_S5_prinseq_good_MBZg_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU118_S5_prinseq_good_MBZg_bothEnds_Unmapped_sorted.bam -fq HU118_S5_prinseq_good_MBZg_host_removed_r1.fastq -fq2 HU118_S5_prinseq_good_MBZg_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU119_S20_R1_001_prinseq_good_Kh_O.fastq.PE -2 ../HU119_S20_R2_001_prinseq_good_NSH3.fastq.PE -S HU119_S20_prinseq_good_NSH3_mapped_and_unmapped.sam
samtools view -bS HU119_S20_prinseq_good_NSH3_mapped_and_unmapped.sam > HU119_S20_prinseq_good_NSH3_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU119_S20_prinseq_good_NSH3_mapped_and_unmapped.bam > HU119_S20_prinseq_good_NSH3_bothEnds_Unmapped.bam
samtools sort -n HU119_S20_prinseq_good_NSH3_bothEnds_Unmapped.bam > HU119_S20_prinseq_good_NSH3_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU119_S20_prinseq_good_NSH3_bothEnds_Unmapped_sorted.bam -fq HU119_S20_prinseq_good_NSH3_host_removed_r1.fastq -fq2 HU119_S20_prinseq_good_NSH3_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU122_S6_R1_001_prinseq_good_NRrq.fastq.PE -2 ../HU122_S6_R2_001_prinseq_good_jNH_.fastq.PE -S HU122_S6_prinseq_good_jNH__mapped_and_unmapped.sam
samtools view -bS HU122_S6_prinseq_good_jNH__mapped_and_unmapped.sam > HU122_S6_prinseq_good_jNH__mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU122_S6_prinseq_good_jNH__mapped_and_unmapped.bam > HU122_S6_prinseq_good_jNH__bothEnds_Unmapped.bam
samtools sort -n HU122_S6_prinseq_good_jNH__bothEnds_Unmapped.bam > HU122_S6_prinseq_good_jNH__bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU122_S6_prinseq_good_jNH__bothEnds_Unmapped_sorted.bam -fq HU122_S6_prinseq_good_jNH__host_removed_r1.fastq -fq2 HU122_S6_prinseq_good_jNH__host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU123_S7_R1_001_prinseq_good_aikJ.fastq.PE -2 ../HU123_S7_R2_001_prinseq_good_aY1_.fastq.PE -S HU123_S7_prinseq_good_aY1__mapped_and_unmapped.sam
samtools view -bS HU123_S7_prinseq_good_aY1__mapped_and_unmapped.sam > HU123_S7_prinseq_good_aY1__mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU123_S7_prinseq_good_aY1__mapped_and_unmapped.bam > HU123_S7_prinseq_good_aY1__bothEnds_Unmapped.bam
samtools sort -n HU123_S7_prinseq_good_aY1__bothEnds_Unmapped.bam > HU123_S7_prinseq_good_aY1__bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU123_S7_prinseq_good_aY1__bothEnds_Unmapped_sorted.bam -fq HU123_S7_prinseq_good_aY1__host_removed_r1.fastq -fq2 HU123_S7_prinseq_good_aY1__host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU124_S23_R1_001_prinseq_good_zHKO.fastq.PE -2 ../HU124_S23_R2_001_prinseq_good_k2jQ.fastq.PE -S HU124_S23_prinseq_good_k2jQ_mapped_and_unmapped.sam
samtools view -bS HU124_S23_prinseq_good_k2jQ_mapped_and_unmapped.sam > HU124_S23_prinseq_good_k2jQ_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU124_S23_prinseq_good_k2jQ_mapped_and_unmapped.bam > HU124_S23_prinseq_good_k2jQ_bothEnds_Unmapped.bam
samtools sort -n HU124_S23_prinseq_good_k2jQ_bothEnds_Unmapped.bam > HU124_S23_prinseq_good_k2jQ_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU124_S23_prinseq_good_k2jQ_bothEnds_Unmapped_sorted.bam -fq HU124_S23_prinseq_good_k2jQ_host_removed_r1.fastq -fq2 HU124_S23_prinseq_good_k2jQ_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU500_S8_R1_001_prinseq_good_HWRS.fastq.PE -2 ../HU500_S8_R2_001_prinseq_good_O0kA.fastq.PE -S HU500_S8_prinseq_good_O0kA_mapped_and_unmapped.sam
samtools view -bS HU500_S8_prinseq_good_O0kA_mapped_and_unmapped.sam > HU500_S8_prinseq_good_O0kA_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU500_S8_prinseq_good_O0kA_mapped_and_unmapped.bam > HU500_S8_prinseq_good_O0kA_bothEnds_Unmapped.bam
samtools sort -n HU500_S8_prinseq_good_O0kA_bothEnds_Unmapped.bam > HU500_S8_prinseq_good_O0kA_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU500_S8_prinseq_good_O0kA_bothEnds_Unmapped_sorted.bam -fq HU500_S8_prinseq_good_O0kA_host_removed_r1.fastq -fq2 HU500_S8_prinseq_good_O0kA_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU502_S24_R1_001_prinseq_good_5oX9.fastq.PE -2 ../HU502_S24_R2_001_prinseq_good_H9AM.fastq.PE -S HU502_S24_prinseq_good_H9AM_mapped_and_unmapped.sam
samtools view -bS HU502_S24_prinseq_good_H9AM_mapped_and_unmapped.sam > HU502_S24_prinseq_good_H9AM_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU502_S24_prinseq_good_H9AM_mapped_and_unmapped.bam > HU502_S24_prinseq_good_H9AM_bothEnds_Unmapped.bam
samtools sort -n HU502_S24_prinseq_good_H9AM_bothEnds_Unmapped.bam > HU502_S24_prinseq_good_H9AM_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU502_S24_prinseq_good_H9AM_bothEnds_Unmapped_sorted.bam -fq HU502_S24_prinseq_good_H9AM_host_removed_r1.fastq -fq2 HU502_S24_prinseq_good_H9AM_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU503_S9_R1_001_prinseq_good_FefH.fastq.PE -2 ../HU503_S9_R2_001_prinseq_good_Sw1C.fastq.PE -S HU503_S9_prinseq_good_Sw1C_mapped_and_unmapped.sam
samtools view -bS HU503_S9_prinseq_good_Sw1C_mapped_and_unmapped.sam > HU503_S9_prinseq_good_Sw1C_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU503_S9_prinseq_good_Sw1C_mapped_and_unmapped.bam > HU503_S9_prinseq_good_Sw1C_bothEnds_Unmapped.bam
samtools sort -n HU503_S9_prinseq_good_Sw1C_bothEnds_Unmapped.bam > HU503_S9_prinseq_good_Sw1C_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU503_S9_prinseq_good_Sw1C_bothEnds_Unmapped_sorted.bam -fq HU503_S9_prinseq_good_Sw1C_host_removed_r1.fastq -fq2 HU503_S9_prinseq_good_Sw1C_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU504_S25_R1_001_prinseq_good_KXbX.fastq.PE -2 ../HU504_S25_R2_001_prinseq_good_zYS3.fastq.PE -S HU504_S25_prinseq_good_zYS3_mapped_and_unmapped.sam
samtools view -bS HU504_S25_prinseq_good_zYS3_mapped_and_unmapped.sam > HU504_S25_prinseq_good_zYS3_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU504_S25_prinseq_good_zYS3_mapped_and_unmapped.bam > HU504_S25_prinseq_good_zYS3_bothEnds_Unmapped.bam
samtools sort -n HU504_S25_prinseq_good_zYS3_bothEnds_Unmapped.bam > HU504_S25_prinseq_good_zYS3_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU504_S25_prinseq_good_zYS3_bothEnds_Unmapped_sorted.bam -fq HU504_S25_prinseq_good_zYS3_host_removed_r1.fastq -fq2 HU504_S25_prinseq_good_zYS3_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU505_S21_R1_001_prinseq_good_2Vgk.fastq.PE -2 ../HU505_S21_R2_001_prinseq_good_rtT4.fastq.PE -S HU505_S21_prinseq_good_rtT4_mapped_and_unmapped.sam
samtools view -bS HU505_S21_prinseq_good_rtT4_mapped_and_unmapped.sam > HU505_S21_prinseq_good_rtT4_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU505_S21_prinseq_good_rtT4_mapped_and_unmapped.bam > HU505_S21_prinseq_good_rtT4_bothEnds_Unmapped.bam
samtools sort -n HU505_S21_prinseq_good_rtT4_bothEnds_Unmapped.bam > HU505_S21_prinseq_good_rtT4_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU505_S21_prinseq_good_rtT4_bothEnds_Unmapped_sorted.bam -fq HU505_S21_prinseq_good_rtT4_host_removed_r1.fastq -fq2 HU505_S21_prinseq_good_rtT4_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU506_S10_R1_001_prinseq_good_KcHI.fastq.PE -2 ../HU506_S10_R2_001_prinseq_good_pXwj.fastq.PE -S HU506_S10_prinseq_good_pXwj_mapped_and_unmapped.sam
samtools view -bS HU506_S10_prinseq_good_pXwj_mapped_and_unmapped.sam > HU506_S10_prinseq_good_pXwj_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU506_S10_prinseq_good_pXwj_mapped_and_unmapped.bam > HU506_S10_prinseq_good_pXwj_bothEnds_Unmapped.bam
samtools sort -n HU506_S10_prinseq_good_pXwj_bothEnds_Unmapped.bam > HU506_S10_prinseq_good_pXwj_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU506_S10_prinseq_good_pXwj_bothEnds_Unmapped_sorted.bam -fq HU506_S10_prinseq_good_pXwj_host_removed_r1.fastq -fq2 HU506_S10_prinseq_good_pXwj_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU508_S26_R1_001_prinseq_good_O5D0.fastq.PE -2 ../HU508_S26_R2_001_prinseq_good_0ext.fastq.PE -S HU508_S26_prinseq_good_0ext_mapped_and_unmapped.sam
samtools view -bS HU508_S26_prinseq_good_0ext_mapped_and_unmapped.sam > HU508_S26_prinseq_good_0ext_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU508_S26_prinseq_good_0ext_mapped_and_unmapped.bam > HU508_S26_prinseq_good_0ext_bothEnds_Unmapped.bam
samtools sort -n HU508_S26_prinseq_good_0ext_bothEnds_Unmapped.bam > HU508_S26_prinseq_good_0ext_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU508_S26_prinseq_good_0ext_bothEnds_Unmapped_sorted.bam -fq HU508_S26_prinseq_good_0ext_host_removed_r1.fastq -fq2 HU508_S26_prinseq_good_0ext_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU537_S11_R1_001_prinseq_good_7hPA.fastq.PE -2 ../HU537_S11_R2_001_prinseq_good_DLJi.fastq.PE -S HU537_S11_prinseq_good_DLJi_mapped_and_unmapped.sam
samtools view -bS HU537_S11_prinseq_good_DLJi_mapped_and_unmapped.sam > HU537_S11_prinseq_good_DLJi_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU537_S11_prinseq_good_DLJi_mapped_and_unmapped.bam > HU537_S11_prinseq_good_DLJi_bothEnds_Unmapped.bam
samtools sort -n HU537_S11_prinseq_good_DLJi_bothEnds_Unmapped.bam > HU537_S11_prinseq_good_DLJi_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU537_S11_prinseq_good_DLJi_bothEnds_Unmapped_sorted.bam -fq HU537_S11_prinseq_good_DLJi_host_removed_r1.fastq -fq2 HU537_S11_prinseq_good_DLJi_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU538_S12_R1_001_prinseq_good_hxY1.fastq.PE -2 ../HU538_S12_R2_001_prinseq_good_DEtR.fastq.PE -S HU538_S12_prinseq_good_DEtR_mapped_and_unmapped.sam
samtools view -bS HU538_S12_prinseq_good_DEtR_mapped_and_unmapped.sam > HU538_S12_prinseq_good_DEtR_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU538_S12_prinseq_good_DEtR_mapped_and_unmapped.bam > HU538_S12_prinseq_good_DEtR_bothEnds_Unmapped.bam
samtools sort -n HU538_S12_prinseq_good_DEtR_bothEnds_Unmapped.bam > HU538_S12_prinseq_good_DEtR_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU538_S12_prinseq_good_DEtR_bothEnds_Unmapped_sorted.bam -fq HU538_S12_prinseq_good_DEtR_host_removed_r1.fastq -fq2 HU538_S12_prinseq_good_DEtR_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU539_S13_R1_001_prinseq_good_FaKf.fastq.PE -2 ../HU539_S13_R2_001_prinseq_good_RsRp.fastq.PE -S HU539_S13_prinseq_good_RsRp_mapped_and_unmapped.sam
samtools view -bS HU539_S13_prinseq_good_RsRp_mapped_and_unmapped.sam > HU539_S13_prinseq_good_RsRp_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU539_S13_prinseq_good_RsRp_mapped_and_unmapped.bam > HU539_S13_prinseq_good_RsRp_bothEnds_Unmapped.bam
samtools sort -n HU539_S13_prinseq_good_RsRp_bothEnds_Unmapped.bam > HU539_S13_prinseq_good_RsRp_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU539_S13_prinseq_good_RsRp_bothEnds_Unmapped_sorted.bam -fq HU539_S13_prinseq_good_RsRp_host_removed_r1.fastq -fq2 HU539_S13_prinseq_good_RsRp_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU540_S14_R1_001_prinseq_good_uxdX.fastq.PE -2 ../HU540_S14_R2_001_prinseq_good_EiSO.fastq.PE -S HU540_S14_prinseq_good_EiSO_mapped_and_unmapped.sam
samtools view -bS HU540_S14_prinseq_good_EiSO_mapped_and_unmapped.sam > HU540_S14_prinseq_good_EiSO_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU540_S14_prinseq_good_EiSO_mapped_and_unmapped.bam > HU540_S14_prinseq_good_EiSO_bothEnds_Unmapped.bam
samtools sort -n HU540_S14_prinseq_good_EiSO_bothEnds_Unmapped.bam > HU540_S14_prinseq_good_EiSO_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU540_S14_prinseq_good_EiSO_bothEnds_Unmapped_sorted.bam -fq HU540_S14_prinseq_good_EiSO_host_removed_r1.fastq -fq2 HU540_S14_prinseq_good_EiSO_host_removed_r2.fastq
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU99_S16_R1_001_prinseq_good_rT3P.fastq.PE -2 ../HU99_S16_R2_001_prinseq_good_gv68.fastq.PE -S HU99_S16_prinseq_good_gv68_mapped_and_unmapped.sam
samtools view -bS HU99_S16_prinseq_good_gv68_mapped_and_unmapped.sam > HU99_S16_prinseq_good_gv68_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU99_S16_prinseq_good_gv68_mapped_and_unmapped.bam > HU99_S16_prinseq_good_gv68_bothEnds_Unmapped.bam
samtools sort -n HU99_S16_prinseq_good_gv68_bothEnds_Unmapped.bam > HU99_S16_prinseq_good_gv68_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU99_S16_prinseq_good_gv68_bothEnds_Unmapped_sorted.bam -fq HU99_S16_prinseq_good_gv68_host_removed_r1.fastq -fq2 HU99_S16_prinseq_good_gv68_host_removed_r2.fastq



cd ~/Metagenomic/Gorilla/N_Trimmed/Homopolymer/Trimmomatic/Host_removal

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2009_106_CGATGT_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2009_106_CGATGT_L001_R2_001.Combined.fastq_trimmed.PE -S 2009_106_CGATGT_mapped_and_unmapped.sam 
samtools view -bS 2009_106_CGATGT_mapped_and_unmapped.sam > 2009_106_CGATGT_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2009_106_CGATGT_mapped_and_unmapped.bam > 2009_106_CGATGTbothEndsUnmapped.bam 
samtools sort -n 2009_106_CGATGTbothEndsUnmapped.bam > 2009_106_CGATGTbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2009_106_CGATGTbothEndsUnmapped_sorted.bam -fq 2009_106_CGATGT_host_removed_r1.fastq -fq2 2009_106_CGATGT_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2009_107_TGACCA_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2009_107_TGACCA_L001_R2_001.Combined.fastq_trimmed.PE -S 2009_107_TGACCA_mapped_and_unmapped.sam 
samtools view -bS 2009_107_TGACCA_mapped_and_unmapped.sam > 2009_107_TGACCA_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2009_107_TGACCA_mapped_and_unmapped.bam > 2009_107_TGACCAbothEndsUnmapped.bam 
samtools bowtie2sort -n 2009_107_TGACCAbothEndsUnmapped.bam > 2009_107_TGACCAbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2009_107_TGACCAbothEndsUnmapped_sorted.bam -fq 2009_107_TGACCA_host_removed_r1.fastq -fq2 2009_107_TGACCA_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2009_109_ATCACG_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2009_109_ATCACG_L001_R2_001.Combined.fastq_trimmed.PE -S 2009_109_ATCACG_mapped_and_unmapped.sam 
samtools view -bS 2009_109_ATCACG_mapped_and_unmapped.sam > 2009_109_ATCACG_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2009_109_ATCACG_mapped_and_unmapped.bam > 2009_109_ATCACGbothEndsUnmapped.bam 
samtools sort -n 2009_109_ATCACGbothEndsUnmapped.bam > 2009_109_ATCACGbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2009_109_ATCACGbothEndsUnmapped_sorted.bam -fq 2009_109_ATCACG_host_removed_r1.fastq -fq2 2009_109_ATCACG_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2009_110_GGTAGC_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2009_110_GGTAGC_L001_R2_001.Combined.fastq_trimmed.PE -S 2009_110_GGTAGC_mapped_and_unmapped.sam 
samtools view -bS 2009_110_GGTAGC_mapped_and_unmapped.sam > 2009_110_GGTAGC_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2009_110_GGTAGC_mapped_and_unmapped.bam > 2009_110_GGTAGCbothEndsUnmapped.bam 
samtools sort -n 2009_110_GGTAGCbothEndsUnmapped.bam > 2009_110_GGTAGCbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2009_110_GGTAGCbothEndsUnmapped_sorted.bam -fq 2009_110_GGTAGC_host_removed_r1.fastq -fq2 2009_110_GGTAGC_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2009_111_ACTTGA_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2009_111_ACTTGA_L001_R2_001.Combined.fastq_trimmed.PE -S 2009_111_ACTTGA_mapped_and_unmapped.sam 
samtools view -bS 2009_111_ACTTGA_mapped_and_unmapped.sam > 2009_111_ACTTGA_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2009_111_ACTTGA_mapped_and_unmapped.bam > 2009_111_ACTTGAbothEndsUnmapped.bam 
samtools sort -n 2009_111_ACTTGAbothEndsUnmapped.bam > 2009_111_ACTTGAbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2009_111_ACTTGAbothEndsUnmapped_sorted.bam -fq 2009_111_ACTTGA_host_removed_r1.fastq -fq2 2009_111_ACTTGA_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2009_115_GATCAG_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2009_115_GATCAG_L001_R2_001.Combined.fastq_trimmed.PE -S 2009_115_GATCAG_mapped_and_unmapped.sam 
samtools view -bS 2009_115_GATCAG_mapped_and_unmapped.sam > 2009_115_GATCAG_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2009_115_GATCAG_mapped_and_unmapped.bam > 2009_115_GATCAGbothEndsUnmapped.bam 
samtools sort -n 2009_115_GATCAGbothEndsUnmapped.bam > 2009_115_GATCAGbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2009_115_GATCAGbothEndsUnmapped_sorted.bam -fq 2009_115_GATCAG_host_removed_r1.fastq -fq2 2009_115_GATCAG_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2009_118_TAGCTT_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2009_118_TAGCTT_L001_R2_001.Combined.fastq_trimmed.PE -S 2009_118_TAGCTT_mapped_and_unmapped.sam 
samtools view -bS 2009_118_TAGCTT_mapped_and_unmapped.sam > 2009_118_TAGCTT_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2009_118_TAGCTT_mapped_and_unmapped.bam > 2009_118_TAGCTTbothEndsUnmapped.bam 
samtools sort -n 2009_118_TAGCTTbothEndsUnmapped.bam > 2009_118_TAGCTTbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2009_118_TAGCTTbothEndsUnmapped_sorted.bam -fq 2009_118_TAGCTT_host_removed_r1.fastq -fq2 2009_118_TAGCTT_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2009_173_ACAGTG_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2009_173_ACAGTG_L001_R2_001.Combined.fastq_trimmed.PE -S 2009_173_ACAGTG_mapped_and_unmapped.sam 
samtools view -bS 2009_173_ACAGTG_mapped_and_unmapped.sam > 2009_173_ACAGTG_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2009_173_ACAGTG_mapped_and_unmapped.bam > 2009_173_ACAGTGbothEndsUnmapped.bam 
samtools sort -n 2009_173_ACAGTGbothEndsUnmapped.bam > 2009_173_ACAGTGbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2009_173_ACAGTGbothEndsUnmapped_sorted.bam -fq 2009_173_ACAGTG_host_removed_r1.fastq -fq2 2009_173_ACAGTG_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2009_193_GCCAAT_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2009_193_GCCAAT_L001_R2_001.Combined.fastq_trimmed.PE -S 2009_193_GCCAAT_mapped_and_unmapped.sam 
samtools view -bS 2009_193_GCCAAT_mapped_and_unmapped.sam > 2009_193_GCCAAT_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2009_193_GCCAAT_mapped_and_unmapped.bam > 2009_193_GCCAATbothEndsUnmapped.bam 
samtools sort -n 2009_193_GCCAATbothEndsUnmapped.bam > 2009_193_GCCAATbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2009_193_GCCAATbothEndsUnmapped_sorted.bam -fq 2009_193_GCCAAT_host_removed_r1.fastq -fq2 2009_193_GCCAAT_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2009_245_CAGATC_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2009_245_CAGATC_L001_R2_001.Combined.fastq_trimmed.PE -S 2009_245_CAGATC_mapped_and_unmapped.sam 
samtools view -bS 2009_245_CAGATC_mapped_and_unmapped.sam > 2009_245_CAGATC_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2009_245_CAGATC_mapped_and_unmapped.bam > 2009_245_CAGATCbothEndsUnmapped.bam 
samtools sort -n 2009_245_CAGATCbothEndsUnmapped.bam > 2009_245_CAGATCbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2009_245_CAGATCbothEndsUnmapped_sorted.bam -fq 2009_245_CAGATC_host_removed_r1.fastq -fq2 2009_245_CAGATC_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2009_265_CTTGTA_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2009_265_CTTGTA_L001_R2_001.Combined.fastq_trimmed.PE -S 2009_265_CTTGTA_mapped_and_unmapped.sam 
samtools view -bS 2009_265_CTTGTA_mapped_and_unmapped.sam > 2009_265_CTTGTA_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2009_265_CTTGTA_mapped_and_unmapped.bam > 2009_265_CTTGTAbothEndsUnmapped.bam 
samtools sort -n 2009_265_CTTGTAbothEndsUnmapped.bam > 2009_265_CTTGTAbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2009_265_CTTGTAbothEndsUnmapped_sorted.bam -fq 2009_265_CTTGTA_host_removed_r1.fastq -fq2 2009_265_CTTGTA_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2011_109_GTTTCG_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2011_109_GTTTCG_L001_R2_001.Combined.fastq_trimmed.PE -S 2011_109_GTTTCG_mapped_and_unmapped.sam 
samtools view -bS 2011_109_GTTTCG_mapped_and_unmapped.sam > 2011_109_GTTTCG_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2011_109_GTTTCG_mapped_and_unmapped.bam > 2011_109_GTTTCGbothEndsUnmapped.bam 
samtools sort -n 2011_109_GTTTCGbothEndsUnmapped.bam > 2011_109_GTTTCGbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2011_109_GTTTCGbothEndsUnmapped_sorted.bam -fq 2011_109_GTTTCG_host_removed_r1.fastq -fq2 2011_109_GTTTCG_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2011_118_CGTACG_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2011_118_CGTACG_L001_R2_001.Combined.fastq_trimmed.PE -S 2011_118_CGTACG_mapped_and_unmapped.sam 
samtools view -bS 2011_118_CGTACG_mapped_and_unmapped.sam > 2011_118_CGTACG_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2011_118_CGTACG_mapped_and_unmapped.bam > 2011_118_CGTACGbothEndsUnmapped.bam 
samtools sort -n 2011_118_CGTACGbothEndsUnmapped.bam > 2011_118_CGTACGbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2011_118_CGTACGbothEndsUnmapped_sorted.bam -fq 2011_118_CGTACG_host_removed_r1.fastq -fq2 2011_118_CGTACG_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2011_119_GAGTGG_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2011_119_GAGTGG_L001_R2_001.Combined.fastq_trimmed.PE -S 2011_119_GAGTGG_mapped_and_unmapped.sam 
samtools view -bS 2011_119_GAGTGG_mapped_and_unmapped.sam > 2011_119_GAGTGG_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2011_119_GAGTGG_mapped_and_unmapped.bam > 2011_119_GAGTGGbothEndsUnmapped.bam 
samtools sort -n 2011_119_GAGTGGbothEndsUnmapped.bam > 2011_119_GAGTGGbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2011_119_GAGTGGbothEndsUnmapped_sorted.bam -fq 2011_119_GAGTGG_host_removed_r1.fastq -fq2 2011_119_GAGTGG_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2011_149_GTCCGC_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2011_149_GTCCGC_L001_R2_001.Combined.fastq_trimmed.PE -S 2011_149_GTCCGC_mapped_and_unmapped.sam 
samtools view -bS 2011_149_GTCCGC_mapped_and_unmapped.sam > 2011_149_GTCCGC_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2011_149_GTCCGC_mapped_and_unmapped.bam > 2011_149_GTCCGCbothEndsUnmapped.bam 
samtools sort -n 2011_149_GTCCGCbothEndsUnmapped.bam > 2011_149_GTCCGCbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2011_149_GTCCGCbothEndsUnmapped_sorted.bam -fq 2011_149_GTCCGC_host_removed_r1.fastq -fq2 2011_149_GTCCGC_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2011_56_GGCTAC_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2011_56_GGCTAC_L001_R2_001.Combined.fastq_trimmed.PE -S 2011_56_GGCTAC_mapped_and_unmapped.sam 
samtools view -bS 2011_56_GGCTAC_mapped_and_unmapped.sam > 2011_56_GGCTAC_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2011_56_GGCTAC_mapped_and_unmapped.bam > 2011_56_GGCTACbothEndsUnmapped.bam 
samtools sort -n 2011_56_GGCTACbothEndsUnmapped.bam > 2011_56_GGCTACbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2011_56_GGCTACbothEndsUnmapped_sorted.bam -fq 2011_56_GGCTAC_host_removed_r1.fastq -fq2 2011_56_GGCTAC_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2011_60_AGTCAA_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2011_60_AGTCAA_L001_R2_001.Combined.fastq_trimmed.PE -S 2011_60_AGTCAA_mapped_and_unmapped.sam 
samtools view -bS 2011_60_AGTCAA_mapped_and_unmapped.sam > 2011_60_AGTCAA_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2011_60_AGTCAA_mapped_and_unmapped.bam > 2011_60_AGTCAAbothEndsUnmapped.bam 
samtools sort -n 2011_60_AGTCAAbothEndsUnmapped.bam > 2011_60_AGTCAAbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2011_60_AGTCAAbothEndsUnmapped_sorted.bam -fq 2011_60_AGTCAA_host_removed_r1.fastq -fq2 2011_60_AGTCAA_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2011_62_AGTTCC_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2011_62_AGTTCC_L001_R2_001.Combined.fastq_trimmed.PE -S 2011_62_AGTTCC_mapped_and_unmapped.sam 
samtools view -bS 2011_62_AGTTCC_mapped_and_unmapped.sam > 2011_62_AGTTCC_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2011_62_AGTTCC_mapped_and_unmapped.bam > 2011_62_AGTTCCbothEndsUnmapped.bam 
samtools sort -n 2011_62_AGTTCCbothEndsUnmapped.bam > 2011_62_AGTTCCbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2011_62_AGTTCCbothEndsUnmapped_sorted.bam -fq 2011_62_AGTTCC_host_removed_r1.fastq -fq2 2011_62_AGTTCC_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2011_66_ATGTCA_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2011_66_ATGTCA_L001_R2_001.Combined.fastq_trimmed.PE -S 2011_66_ATGTCA_mapped_and_unmapped.sam 
samtools view -bS 2011_66_ATGTCA_mapped_and_unmapped.sam > 2011_66_ATGTCA_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2011_66_ATGTCA_mapped_and_unmapped.bam > 2011_66_ATGTCAbothEndsUnmapped.bam 
samtools sort -n 2011_66_ATGTCAbothEndsUnmapped.bam > 2011_66_ATGTCAbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2011_66_ATGTCAbothEndsUnmapped_sorted.bam -fq 2011_66_ATGTCA_host_removed_r1.fastq -fq2 2011_66_ATGTCA_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2011_69_GTGAAA_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2011_69_GTGAAA_L001_R2_001.Combined.fastq_trimmed.PE -S 2011_69_GTGAAA_mapped_and_unmapped.sam 
samtools view -bS 2011_69_GTGAAA_mapped_and_unmapped.sam > 2011_69_GTGAAA_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2011_69_GTGAAA_mapped_and_unmapped.bam > 2011_69_GTGAAAbothEndsUnmapped.bam 
samtools sort -n 2011_69_GTGAAAbothEndsUnmapped.bam > 2011_69_GTGAAAbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2011_69_GTGAAAbothEndsUnmapped_sorted.bam -fq 2011_69_GTGAAA_host_removed_r1.fastq -fq2 2011_69_GTGAAA_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2011_70_GTGGCC_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2011_70_GTGGCC_L001_R2_001.Combined.fastq_trimmed.PE -S 2011_70_GTGGCC_mapped_and_unmapped.sam 
samtools view -bS 2011_70_GTGGCC_mapped_and_unmapped.sam > 2011_70_GTGGCC_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2011_70_GTGGCC_mapped_and_unmapped.bam > 2011_70_GTGGCCbothEndsUnmapped.bam 
samtools sort -n 2011_70_GTGGCCbothEndsUnmapped.bam > 2011_70_GTGGCCbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2011_70_GTGGCCbothEndsUnmapped_sorted.bam -fq 2011_70_GTGGCC_host_removed_r1.fastq -fq2 2011_70_GTGGCC_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2011_79_CCGTCC_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2011_79_CCGTCC_L001_R2_001.Combined.fastq_trimmed.PE -S 2011_79_CCGTCC_mapped_and_unmapped.sam 
samtools view -bS 2011_79_CCGTCC_mapped_and_unmapped.sam > 2011_79_CCGTCC_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2011_79_CCGTCC_mapped_and_unmapped.bam > 2011_79_CCGTCCbothEndsUnmapped.bam 
samtools sort -n 2011_79_CCGTCCbothEndsUnmapped.bam > 2011_79_CCGTCCbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2011_79_CCGTCCbothEndsUnmapped_sorted.bam -fq 2011_79_CCGTCC_host_removed_r1.fastq -fq2 2011_79_CCGTCC_host_removed_r2.fastq

bowtie2 -x ~/Databases/Gorilla/Gorilla_gorilla.gorGor4.cdna_BowtieIndex -1 ../2011_85_GTAGAG_L001_R1_001.Combined.fastq_trimmed.PE -2 ../2011_85_GTAGAG_L001_R2_001.Combined.fastq_trimmed.PE -S 2011_85_GTAGAG_mapped_and_unmapped.sam 
samtools view -bS 2011_85_GTAGAG_mapped_and_unmapped.sam > 2011_85_GTAGAG_mapped_and_unmapped.bam 
samtools view -b -f 12 -F 256 2011_85_GTAGAG_mapped_and_unmapped.bam > 2011_85_GTAGAGbothEndsUnmapped.bam 
samtools sort -n 2011_85_GTAGAGbothEndsUnmapped.bam > 2011_85_GTAGAGbothEndsUnmapped_sorted.bam 
bedtools bamtofastq -i 2011_85_GTAGAGbothEndsUnmapped_sorted.bam -fq 2011_85_GTAGAG_host_removed_r1.fastq -fq2 2011_85_GTAGAG_host_removed_r2.fastq

