# Author: Ashok Kumar Sharma

#----Figure1: Gene Richness and diversity analsyis

#---Gene Count file used for the rarefaction and diversity analysis
#---Gene Count file used for the rarefaction and diversity analysis
Abundance <- read.csv(file = "../LargeFiles_Used/Gene_count_matrix.tab", row.names=1, header = T, sep = "\t")
colSums(Abundance)
new <- colSums(Abundance)
new <- data.frame(new)
View(new)
# The lowest value is 4498117 and hence, we will rarefy the data for 4.4 million (4400000) This is in one
# Sample named as 2011_60

#biom convert -i Gene_count_matrix.tab -o Gene_count_matrix.tab.biom --table-type="OTU table" --to-json
#multiple_rarefactions_even_depth.py -i Gene_count_matrix.tab.biom -o rarefied_files -d 4400000 -n 10
#alpha_diversity.py -i rarefied_files -m shannon,simpson,chao1,Observed -o alpha_diversity_out 

#perl calculate_mean.pl
#perl addLabels.pl ../../Metadata.txt Alpha_diversity_mean

setwd("~/Work/Metagenomic_Analysis/GeneRichness/Alpha_diversity_at_4.4mdepth_Using_Genecountsmatrix")

Color <- c("skyblue4", "wheat4")
Colors <- c("skyblue4", "skyblue", "wheat4", "wheat")

alpha_div_mean <- read.csv (file="Alpha_diversity_mean.Out", sep = "\t", header = TRUE, row.names = 1)

alpha_div_obsreved1 <- alpha_div_mean[, c(4,7)]
alpha_div_obsreved1 <- melt(alpha_div_obsreved1, id.vars = "Group1")
jpeg("observed1.jpg", height = 3.5, width = 1.8, units = 'in', res = 600)
ggplot(data = alpha_div_obsreved1, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("") + labs(x="",y="Observed genes") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
kruskalmc(alpha_div_obsreved1$value ~ alpha_div_obsreved1$Group1)

alpha_div_shannon1 <- alpha_div_mean[, c(1,7)]
alpha_div_shannon1_melted <- melt(alpha_div_shannon1, id.vars = "Group1")

jpeg("shannon1_new.jpg", height = 3.5, width = 1.8, units = 'in', res = 600)
ggplot(data = alpha_div_shannon1_melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("") + labs(x="",y="Shannon") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
kruskalmc(alpha_div_shannon1_melted$value ~ alpha_div_shannon1_melted$Group1)
