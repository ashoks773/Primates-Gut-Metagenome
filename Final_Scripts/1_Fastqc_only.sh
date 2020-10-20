#!/bin/bash -l        
#PBS -l walltime=12:00:00,nodes=3:ppn=8,mem=20gb 
#PBS -m abe 
#PBS -e fastqc.error
#PBS -o fastqc.out
#PBS -N FastQC

module load fastqc

cd ~/Metagenomic/Human
FASTQ="/home/gomeza/sharm646/Metagenomic/Human"

for F in 'ls *.fastq.gz'
do
    fastqc -o $FASTQ $F
done



