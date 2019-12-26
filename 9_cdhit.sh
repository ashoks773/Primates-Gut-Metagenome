#!/bin/bash -l         
#PBS -l walltime=72:00:00,nodes=2:ppn=24,mem=62gb  
#PBS -m abe  
#PBS -e chdhit.error 
#PBS -o cdhit.out 
 
module load cdhit
 
cd ~/Metagenomic/Total_Genes/Check

cd-hit-est -i ../Combined_Total_Genes.fna -o Combined_Total_Genes.cdhit.fna -c 0.95 -n 10 -aS 0.9 -d 0 -T 48 -M 60000
