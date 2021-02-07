#-- Metagenomic data analysis---- Filtering Steps
#-Fastqc
qsub 1_Fastqc_only.sh

#-Ambiguity filtering
qsub 2_Ambiguity.sh
perl ~/Softwares/NGSQCToolkit_v2.3.3/Trimming/AmbiguityFiltering.pl -i ../HU100_S17_R1_001.fastq -irev ../HU100_S17_R2_001.fastq -c 1 -t5 -t3

#-Homopolymer removal
qsub 3_Prinseq.sh
prinseq-lite.pl -fastq ../HU100_S17_R1_001.fastq_trimmed -fastq2 ../HU100_S17_R2_001.fastq_trimmed -custom_params "AAT 10;T 70%;A 15;G 70%;C 15"

#-Barcode removal
qsub 4_Trimmomatic.sh
java -jar $TRIMMOMATIC/trimmomatic.jar PE -phred33 ../HU100_S17_R1_001_prinseq_good_OLyg.fastq ../HU100_S17_R2_001_prinseq_good_Rtjq.fastq HU100_S17_R1_001_prinseq_good_OLyg.fastq.PE HU100_S17_R1_001_prinseq_good_OLyg.fastq.UN HU100_S17_R2_001_prinseq_good_Rtjq.fastq.PE HU100_S17_R2_001_prinseq_good_Rtjq.fastq.UN ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq2-PE.fa:2:40:15 LEADING:10 TRAILING:10 SLIDINGWINDOW:10:20 MINLEN:80

#-Make a singleton file from unpaired reads
cat HU100_S17_R1_001_prinseq_good_OLyg.fastq.UN HU100_S17_R2_001_prinseq_good_Rtjq.fastq.UN > HU100_S17_R2_001_prinseq_good_Rtjq_Nomates.fastq
Repeat same for all samples

#-Host contaminaion removal
qsub 5_Host_contaRemoval.sh
bowtie2 -x ~/Databases/Human/Homo_sapiens.GRCh38.cdna_BowtieIndex -1 ../HU100_S17_R1_001_prinseq_good_OLyg.fastq.PE -2 ../HU100_S17_R2_001_prinseq_good_Rtjq.fastq.PE -S HU100_S17_prinseq_good_Rtjq_mapped_and_unmapped.sam
samtools view -bS HU100_S17_prinseq_good_Rtjq_mapped_and_unmapped.sam > HU100_S17_prinseq_good_Rtjq_mapped_and_unmapped.bam
samtools view -b -f 12 -F 256 HU100_S17_prinseq_good_Rtjq_mapped_and_unmapped.bam > HU100_S17_prinseq_good_Rtjq_bothEnds_Unmapped.bam
samtools sort -n HU100_S17_prinseq_good_Rtjq_bothEnds_Unmapped.bam > HU100_S17_prinseq_good_Rtjq_bothEnds_Unmapped_sorted.bam
bedtools bamtofastq -i HU100_S17_prinseq_good_Rtjq_bothEnds_Unmapped_sorted.bam -fq HU100_S17_prinseq_good_Rtjq_host_removed_r1.fastq -fq2 HU100_S17_prinseq_good_Rtjq_host_removed_r2.fastq

#-Metagenomic Assembly
#-Spades
qsub 6_Spades_Assembly.sh
/home/gomeza/sharm646/Softwares/SPAdes-3.12.0-Linux/bin/spades.py --meta -1 ../HU100_S17_prinseq_good_Rtjq_host_removed_r1.fastq -2 ../HU100_S17_prinseq_good_Rtjq_host_removed_r2.fastq -s ../HU100_S17_R2_001_prinseq_good_Rtjq_Nomates.fastq -m 950 -t 30 -o HU100_S17_prinseq_good_Rtjq_metaspades_out_result

#-Contigs filtering
perl Filter_contigs_read_length.pl

#-Gene prediction 
qsub 7_Prodigal.sh
prodigal -i HU100_S17_scaffolds.fasta.Filtered.fna -o HU100_S17.Genes -a HU100_S17.proteins.faa -p meta

#-Combined Gene catalogue for Human and Gorilla data
perl 8_index.pl *fna
cat *indexed_GENES.fna > Combined_Human_Indexed.fna
perl 8_index.pl *fna
cat *indexed_GENES.fna > Combined_Gorilla_Indexed.fna
cat Combined_Human_Indexed.fna Combined_Gorilla_Indexed.fna > Combined_Total_Genes.fna

#--Construction of Non-Redundant Gene-set
qsub 9_cdhit.sh
cd-hit-est -i ../Combined_Total_Genes.fna -o Combined_Total_Genes.cdhit.fna -c 0.95 -n 10 -aS 0.9 -d 0 -T 48 -M 60000

#--Gene quantification
http://www.cbs.dtu.dk/courses/27626/Exercises/metagenomic.assembly.php
qsub newBWA.sh

1. bwa index Combined_Total_Genes.cdhit.fna
2. bwa mem -t 144 -M Combined_Total_Genes.cdhit.fna /home/gomeza/sharm646/Metagenomic/Human/N_Trimmed/Homopolymer/Trimmomatic/Host_removal/HU100_S17_prinseq_good_Rtjq_host_removed_r1.fastq /home/gomeza/sharm646/Metagenomic/Human/N_Trimmed/Homopolymer/Trimmomatic/Host_removal/HU100_S17_prinseq_good_Rtjq_host_removed_r2.fastq | samtools view -Sb - > HU100_S17.bam
3. samtools view -h HU100_S17.bam | perl read_count_bam.pl | samtools view -Su -q30 - | samtools sort -O BAM -o HU100_S17.sort.bam -
4. samtools index HU100_S17.sort.bam

sh files 1.sh 2.sh 3.sh

echo -e "Genes  2009_106        2009_107        2009_109        2009_110        2009_111        2009_115        2009_118        2009_173        2009_193        2009_245        2009_265        2011_109        2011_118        2011_119        2011_149        2011_56 2011_60 2011_62 2011_66 2011_69 2011_70 2011_79 2011_85 HU100_S17       HU101_S15       HU103_S22       HU107_S1        HU108_S18       HU109_S19       HU112_S2        HU113_S3        HU114_S27       HU115_S28       HU117_S4        HU118_S5        HU119_S20       HU122_S6        HU123_S7        HU124_S23       HU500_S8        HU502_S24       HU503_S9        HU504_S25       HU505_S21       HU506_S10       HU508_S26       HU537_S11       HU538_S12       HU539_S13       HU540_S14       HU99_S16" > header

samtools idxstats 2009_106.sort.bam | grep -v "*" | cut -f1 > gene_names

samtools idxstats 2009_106.sort.bam | grep -v "*" | cut -f3 > Counts1
samtools idxstats 2009_107.sort.bam | grep -v "*" | cut -f3 > Counts2
So on ...

paste gene_names Counts1 Counts2 Counts3 Counts4 Counts5 Counts6 Counts7 Counts8 Counts9 Counts10 Counts11 Counts12 Counts13 Counts14 Counts15 Counts16 Counts17 Counts18 Counts19 Counts20 Counts21 Counts22 Counts23 Counts24 Counts25 Counts26 Counts27 Counts28 Counts29 Counts30 Counts31 Counts32 Counts33 Counts34 Counts35 Counts36 Counts37 Counts38 Counts39 Counts40 Counts41 Counts42 Counts43 Counts44 Counts45 Counts46 Counts47 Counts48 Counts49 Counts50 Counts51 | cat header - > Gene_count_matrix.tab

#----Convert the total genes fasta file in single line sequences
perl seq_single_line.pl Combined_Total_Genes.cdhit.fna
mv Combined_Total_Genes.cdhit.fna.single Combined_Total_Genes.cdhit.single.fna

perl Calculate_gene_abundance.pl Combined_Total_Genes.cdhit.single.fna Gene_count_matrix.tab 
(# Total NR genes 4298551)

perl Filter_genes.pl GENE_ABUNDANCE_TABLE 
(# Total remained after filtering 3804658)

Rscript Gene_proportions.R
Gene_abundance <- read.csv (file = "GENE_ABUNDANCE_TABLE.filtered", row.names = 1, header = T, sep = "\t")
Gene_abundance_normalized <- Gene_abundance/colSums(Gene_abundance)[col(Gene_abundance)]		
write.table (Gene_abundance_normalized, file = "GENE_PROPORTION_TABLE", sep = "\t")

#---- Gene Richness analysis
1. Total number of genes/Total number of reads (For each sample)
2. Gene abundance table (filtered) was used for the rarefaction analysis (In steps.R)

#--- Pick NR-Genes (amino acid sequences) for functional annotations
perl Get_seqs.new.pl Total_Genes/Combined_Total_Proteins.faa Total_Genes/Combined_Total_Genes.cdhit.fna

#---- Kegg Analysis
/home/gomeza/sharm646/Databases/KEGG_Database
makeblastdb -in KEGG_GENOPEP_200314.single -dbtype prot
blastp -query xaa -db ~/Databases/KEGG_Database/KEGG_GENOPEP_200314.single -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -out xaa.Blast_out -num_thre
ads 24 (Query database was divided into multiple parts)

#-- Select tophit and then
cat *Blast_out > Total_NR_Genes_KeggOut
perl Extract_identity.pl Total_NR_Genes_KeggOut
(# 242349 Uniq kegg genes identified ((6.3% of total NR gene set))
perl genes_to_ko_new.pl Total_NR_Genes_KeggOut.besthit 
perl Script_add_EC_ID.pl ~/Databases/KEGG_Database/KEGG_2014_ko.list Total_NR_Genes_KeggOut.besthit.KO
perl Calculate_KO_Abundance.pl Total_NR_Genes_KeggOut.besthit.KO.EC_added ../Total_Genes/Gene_quantification/GENE_PROPORTION_TABLE
perl Calculate_Pathway_Abundance.pl KO_ABUNDANCE

Statistical analysis
Kegg.R

#-------------------CAZyDB_130918 Analysis
/home/gomeza/sharm646/Databases/CaZyDB_130918
cd-hit -i CAZyDB.07202017.fa -o CAZyDB.07202017.cdhit.fa -c 1 -n 5 -aS 1 -d 0 -T 48 -M 60000
makeblastdb -in CAZyDB.07202017.cdhit.fa -parse_seqids -dbtype prot

#Combined_Total_Proteins.NRset.faa divided in to multiple parts for annotation

blastp -query ~/Metagenomic/Final_NR_GeneSet/xae1 -db ~/Databases/CaZyDB_130918/CAZyDB.07202017.cdhit.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -out xae1.Blast_out -num_threads 24

#-- Select tophit and then
cat xaa xab ..... > Total_against_Cazy.BlastOut

perl Extract_identity.pl Total_against_Cazy.BlastOut
(# 19783 Uniq CaZy genes identified ((0.52% of total NR gene set))
perl genes_to_cazy.pl
perl Calculate_cazy_abundances.pl Total_against_Cazy.BlastOut.besthit.CazyAnot ../Total_Genes/Gene_quantification/GENE_PROPORTION_TABLE

Statistical analysis
Cazy.R

#------ Antibiotic resistance Analysis
#Database was downloaded from http://www.dantaslab.org/resfams/
#hmmpress ~/Databases/RESFAMS/Resfams-full.hmm
#hmmscan --tblout xaa.fa_hmmscan.out.tabular -o xaa.fa_hmmscan.out --cpu 120 /home/gomeza/sharm646/Databases/RESFAMS/Resfams-full.hmm xaa.fa

#-- Download ARDB database (7828 Genes -- Removed three genes having same annotations)
blastp -query ~/Metagenomic/Antibiotic_Resistance/xaa -db ~/Databases/ardbAnno1.0/resisGenes.pfasta.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs ppos" -out xaa.Blast_out -num_threads 24
for i in *Blast_out; do perl tophit.pl "$i"; done;
for i in *TopHit; do perl Extract_identity_ARDB.pl "$i"; done;
cat *besthit > Total_against_ARDB.BlastOut

perl genes_to_ardb_gene.pl
perl Calculate_ARDB_abundances.pl Total_against_ARDB.BlastOut.ARDB_anot ../Total_Genes/Gene_quantification/GENE_PROPORTION_TABLE

Statistical analysis
ARDB.R

#------ Xenobiotics Degradation Analysis
blastp -query ~/Metagenomic/Xenobiotics_analysis/xaa -db ~/Metagenomic/Xenobiotics_analysis/protein.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" -out xaa.Blast_out -num_threads 24
cat *Blast_out > Total_Xeno_BlastOut
perl tophit.pl Total_Xeno_BlastOut
perl Extract_identity_Xeno.pl Total_Xeno_BlastOut.Tophit 
cut -f1,2 Total_Xeno_BlastOut.TopHit.besthit > new
cut -f2 new | cut -d"|" -f2 > new1
paste new new1 > Gene_xeno_annot
perl Calculate_XenoEnzyme_abundances.pl Gene_xeno_annot ../Total_Genes/Gene_quantification/GENE_PROPORTION_TABLE

Statistical analysis
xeno.R

#-- Download ARDB database (7828 Genes -- Removed three genes having same annotations)
#--ppos in the command to get similarity (%) in the last column
blastp -query ~/Metagenomic/Antibiotic_Resistance/xaa -db ~/Databases/ardbAnno1.0/resisGenes.pfasta.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs ppos" -out xaa.Blast_out -num_threads 24
for i in *Blast_out; do perl tophit.pl "$i"; done;
for i in *TopHit; do perl Extract_identity_ARDB.pl "$i"; done;
cat *besthit > Total_against_ARDB.BlastOut

perl genes_to_ardb_gene.pl
perl Calculate_ARDB_abundances.pl Total_against_ARDB.BlastOut.ARDB_anot ../Total_Genes/Gene_quantification/GENE_PROPORTION_TABLE

Statistical analysis
ARDB.R

#--- Taxonomic analysis Using HMP+NCBI
blastn -query Final_H_MN_Genes_nt_seq.fasta -db REFERENCE_GENOMES_DATABASE_NCBI_HMP/Reference_Genome_Database -outfmt 6 -num_threads 50 -out Final_H_MN_Genes_nt_NCBI_HMP_BLASTOUT&
perl 1_select_besthit.pl  <Final_H_MN_Genes_nt_NCBI_HMP_BLASTOUT>    # Select best hits
perl Script_Annotate_Hits.pl <ALL_GENOMES_ID> <Combined_Gene_Catalogue.fna> <Final_H_MN_Genes_nt_NCBI_HMP_BLASTOUT.score>     # Annotate the Hits using th eScript
perl Script_for_LCA_Hits.pl <Taxon_Lineage_Final.txt> <Hits_Annotation.out>         # Annotate using LCA Method for Multiple hits
perl Script_FINAL_Assignment_Hits.pl <Taxon_Lineage_Final.txt> <Taxonomy_Assignment_LCA_Hits>       # Perform Final Assignment of Genes against Reference Genomes

Statistical analysis
Taxonomy.R

#--- Genome Reconstruction
1. Cross assembly
qsub Megahit_human_assembly.sh

2. Filter contigs
nvi-script-reformat-fasta contigs.fa -o contigs-fixed.fa -l 1000 --simplify-names

3. Calculate depth using metabat
jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam

4. Metabat binning
metabat2 -i contigs-fixed.fa -a depth.txt -o bins_dir/bin

5. CheckM for quality assessment of Bins
qsub checkM.sh

6. Selected bins and reference genomes for phylogenetic tree construction using marker genes
~/Softwares/ezTree-master/ezTree -list list_species -out genome.out -thread 120 -evalue 1e-5
# list_species contains name of the fasta files to be processed


#------ CaZy analysis on Reconstructed Prevotella and Treponema genomes
1. Gene prediction in the reconstructed Genomes
2. Genes were mapped against CAZyDB_130918 database
3. Counts of CAZymes in each bin were calculated based the best hits
Finally this scirpt was used to generate the final figures
CaZy_Ana_On_ReconstructedGenomes.R

#---- End
