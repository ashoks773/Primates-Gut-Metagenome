#!/bin/bash -l         
#PBS -l walltime=24:00:00,nodes=6:ppn=24,mem=62gb  
#PBS -m abe  
#PBS -e bwa2.error 
#PBS -o bwa2.out 

module load bwa
module load samtools

cd ~/Metagenomic/Total_Genes/Gene_quantification

bwa index Combined_Total_Genes.cdhit.fna

samtools view -h HU100_S17.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU100_S17.sort.bam -
samtools view -h HU101_S15.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU101_S15.sort.bam -
samtools view -h HU103_S22.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU103_S22.sort.bam -
samtools view -h HU107_S1.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU107_S1.sort.bam -
samtools view -h HU108_S18.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU108_S18.sort.bam -
samtools view -h HU109_S19.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU109_S19.sort.bam -
samtools view -h HU112_S2.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU112_S2.sort.bam -
samtools view -h HU113_S3.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU113_S3.sort.bam -
samtools view -h HU114_S27.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU114_S27.sort.bam -
samtools view -h HU115_S28.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU115_S28.sort.bam -
samtools view -h HU117_S4.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU117_S4.sort.bam -
samtools view -h HU118_S5.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU118_S5.sort.bam -
samtools view -h HU119_S20.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU119_S20.sort.bam -
samtools view -h HU122_S6.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU122_S6.sort.bam -
samtools view -h HU123_S7.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU123_S7.sort.bam -
samtools view -h HU124_S23.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU124_S23.sort.bam -
samtools view -h HU500_S8.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU500_S8.sort.bam -
samtools view -h HU502_S24.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU502_S24.sort.bam -
samtools view -h HU503_S9.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU503_S9.sort.bam -
samtools view -h HU504_S25.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU504_S25.sort.bam -
samtools view -h HU505_S21.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU505_S21.sort.bam -
samtools view -h HU506_S10.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU506_S10.sort.bam -
samtools view -h HU508_S26.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU508_S26.sort.bam -
samtools view -h HU537_S11.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU537_S11.sort.bam -
samtools view -h HU538_S12.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU538_S12.sort.bam -
samtools view -h HU539_S13.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU539_S13.sort.bam -
samtools view -h HU540_S14.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU540_S14.sort.bam -
samtools view -h HU99_S16.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU99_S16.sort.bam -
samtools view -h 2009_106.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2009_106.sort.bam -
samtools view -h 2009_107.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2009_107.sort.bam -
samtools view -h 2009_109.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2009_109.sort.bam -
samtools view -h 2009_110.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2009_110.sort.bam -
samtools view -h 2009_111.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2009_111.sort.bam -
samtools view -h 2009_115.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2009_115.sort.bam -
samtools view -h 2009_118.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2009_118.sort.bam -
samtools view -h 2009_173.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2009_173.sort.bam -
samtools view -h 2009_193.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2009_193.sort.bam -
samtools view -h 2009_245.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2009_245.sort.bam -
samtools view -h 2009_265.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2009_265.sort.bam -
samtools view -h 2011_109.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2011_109.sort.bam -
samtools view -h 2011_118.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2011_118.sort.bam -
samtools view -h 2011_119.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2011_119.sort.bam -
samtools view -h 2011_149.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2011_149.sort.bam -
samtools view -h 2011_56.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2011_56.sort.bam -
samtools view -h 2011_60.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2011_60.sort.bam -
samtools view -h 2011_62.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2011_62.sort.bam -
samtools view -h 2011_66.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2011_66.sort.bam -
samtools view -h 2011_69.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2011_69.sort.bam -
samtools view -h 2011_70.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2011_70.sort.bam -
samtools view -h 2011_79.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2011_79.sort.bam -
samtools view -h 2011_85.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o 2011_85.sort.bam -

samtools index HU100_S17.sort.bam
samtools index HU101_S15.sort.bam
samtools index HU103_S22.sort.bam
samtools index HU107_S1.sort.bam
samtools index HU108_S18.sort.bam
samtools index HU109_S19.sort.bam
samtools index HU112_S2.sort.bam
samtools index HU113_S3.sort.bam
samtools index HU114_S27.sort.bam
samtools index HU115_S28.sort.bam
samtools index HU117_S4.sort.bam
samtools index HU118_S5.sort.bam
samtools index HU119_S20.sort.bam
samtools index HU122_S6.sort.bam
samtools index HU123_S7.sort.bam
samtools index HU124_S23.sort.bam
samtools index HU500_S8.sort.bam
samtools index HU502_S24.sort.bam
samtools index HU503_S9.sort.bam
samtools index HU504_S25.sort.bam
samtools index HU505_S21.sort.bam
samtools index HU506_S10.sort.bam
samtools index HU508_S26.sort.bam
samtools index HU537_S11.sort.bam
samtools index HU538_S12.sort.bam
samtools index HU539_S13.sort.bam
samtools index HU540_S14.sort.bam
samtools index HU99_S16.sort.bam
samtools index 2009_106.sort.bam
samtools index 2009_107.sort.bam
samtools index 2009_109.sort.bam
samtools index 2009_110.sort.bam
samtools index 2009_111.sort.bam
samtools index 2009_115.sort.bam
samtools index 2009_118.sort.bam
samtools index 2009_173.sort.bam
samtools index 2009_193.sort.bam
samtools index 2009_245.sort.bam
samtools index 2009_265.sort.bam
samtools index 2011_109.sort.bam
samtools index 2011_118.sort.bam
samtools index 2011_119.sort.bam
samtools index 2011_149.sort.bam
samtools index 2011_56.sort.bam
samtools index 2011_60.sort.bam
samtools index 2011_62.sort.bam
samtools index 2011_66.sort.bam
samtools index 2011_69.sort.bam
samtools index 2011_70.sort.bam
samtools index 2011_79.sort.bam
samtools index 2011_85.sort.bam
