Gene_abundance <- read.csv (file = "GENE_ABUNDANCE_TABLE.filtered", row.names = 1, header = T, sep = "\t")
Gene_abundance_normalized <- Gene_abundance/colSums(Gene_abundance)[col(Gene_abundance)]		
write.table (Gene_abundance_normalized, file = "GENE_PROPORTION_TABLE", sep = "\t")
