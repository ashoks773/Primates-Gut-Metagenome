#!/bin/bash -l
#PBS -l walltime=24:00:00,nodes=1:ppn=24,mem=240gb 
#PBS -m abe 
#PBS -e CheckM.error
#PBS -o CheckM.out

cd ~/Metagenomic/Genome_Reconstruction/Megahit_Assembly/Metabat_Binning

#-- https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow
#Flow1---- Steps for Lineage-specific Workflow
#checkm tree bins_dir/bin bins_dir/bin_checkM -x .fa -t 24 [By default bins should be in .fna format if not then use -x to specify the format of bins]
checkm tree_qa bin_checkM -f marker_file
checkm lineage_set bin_checkM marker_file
checkm analyze marker_file bin bin_checkM_analyse -x .fa -t 24
checkm qa marker_file bin_checkM_analyse -o 1 -t 24
#Output
#Bin Id               Marker lineage             # genomes   # markers   # marker sets    0     1     2     3    4    5+   Completeness   Contamination   Strain heterogeneity  
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#bin.467     f__Bifidobacteriaceae (UID1462)         65         476           217         1    474    1     0    0    0       99.54            0.46              100.00   


#----Make a folder using selected bins (Completeness > 70 and Contamination < 10)
# Selected_bins = 130
#---- Bin merging
checkm taxon_set domain Bacteria bacterial.ms
checkm merge bacterial.ms selected_bins Merged_out -x .fa -t 24

#---Quality check
checkm unique selected_bins
checkm outliers bin_checkM_analyse selected_bins ../tetra_freq.txt outliers.tsv -x .fa

#---- To make plots
checkm bin_qa_plot --image_type png bin_checkM_analyse bin plot_folder -x .fa --dpi 75

#Flow2----- Steps for Taxonomic-specific Workflow 
checkm taxon_list --rank order  
checkm taxon_set order Bacteroidales Bacteroidales.ms 
checkm analyze Bacteroidales.ms bin bin_checkM_Bacteroidales -x .fa -t 24 
checkm qa marker_file bin_checkM_Bacteroidales
