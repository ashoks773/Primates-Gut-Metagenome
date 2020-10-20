# Step 1: Principle coordinate analyis using Pathway_abundance_Table
Pathway <- read.csv(file = "PATHWAY_ABUNDANCE_TABLE", row.names=1, header = T, sep = "\t")
Pathway_proportions <- Pathway/colSums(Pathway)[col(Pathway)]
Pathway1 <- data.frame(t(Pathway_proportions))
write.table (Pathway1, file = "PATHWAY_ABUNDANCE_TABLE_relative_abundances", sep = "\t")
metadata <- read.csv (file = "../Metadata.txt", row.names=1, header = T, sep = "\t")

pathway_bray_pcoa <-pcoa(vegdist(Pathway1, "bray"))
Groups <- metadata[,1:3]
pathway_bray_pcoa$values[1:2,]
mds.var.per = round(pathway_bray_pcoa$values$Eigenvalues/sum(pathway_bray_pcoa$values$Eigenvalues)*100, 1)
Pathway_Bray_PCoA_MATRIX <- pathway_bray_pcoa$vectors[,1:2]
Pathway_Bray_PCoA_MATRIX <- data.frame(Pathway_Bray_PCoA_MATRIX)
Pathway_Bray_PCoA_MATRIX_New <- cbind(Pathway_Bray_PCoA_MATRIX, Groups)

Color <- c("darkgreen", "darkred")
Colors <- c("brown4", "orangered3", "green4", "limegreen")

jpeg("Path_Bray_PCoA_Group.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Pathway_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
jpeg("Path_Bray_PCoA_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Pathway_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCoA1 - ", mds.var.per[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

distance_bray <-vegdist(Pathway1, "bray")
adonis(distance_bray ~ Group*Group1, data=metadata, permutations=99)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#Group      1  0.050977 0.050977  37.087 0.37857   0.01 **
 # Group1     2  0.019079 0.009540   6.940 0.14169   0.01 **

pcoa1 <- Pathway_Bray_PCoA_MATRIX_New[,c(1,4)]
pcoa1_Melted <- melt(pcoa1, id.vars = "Group1")
jpeg("Bray_PCoA1_Distances.jpg", height = 5, width = 5, units = 'in', res = 600)
ggplot(data = pcoa1_Melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA1") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + coord_flip() + theme(legend.position='none')
dev.off ()
pcoa2 <- Pathway_Bray_PCoA_MATRIX_New[,c(2,4)]
pcoa2_Melted <- melt(pcoa2, id.vars = "Group1")
jpeg("Bray_PCoA2_Distances.jpg", height = 5, width = 5, units = 'in', res = 600)
ggplot(data = pcoa2_Melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA2") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()

#---Gorilla only
metadata_gorilla <- read.csv (file = "../Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
pathway_gorilla <- Pathway1[1:23,]
pathway_bray_pcoa_gorilla <-pcoa(vegdist(pathway_gorilla, "bray"))
pathway_bray_pcoa_gorilla$values[1:2,]
mds.var.per_g = round(pathway_bray_pcoa_gorilla$values$Eigenvalues/sum(pathway_bray_pcoa_gorilla$values$Eigenvalues)*100, 1)
Pathway_Bray_PCoA_MATRIX_G <- pathway_bray_pcoa_gorilla$vectors[,1:2]
Pathway_Bray_PCoA_MATRIX_G <- data.frame(Pathway_Bray_PCoA_MATRIX_G)
Pathway_Bray_PCoA_MATRIX_G_New <- cbind(Pathway_Bray_PCoA_MATRIX_G, Groups)

Color <- c("green4", "limegreen")
jpeg("Bray_PCoA_Gorilla_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Pathway_Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
jpeg("Bray_PCoA_Gorilla_Group2.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Pathway_Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1, shape=Group2)) + geom_point(size=3) + geom_point(aes(shape=Group2)) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted") + geom_point(shape=18)
dev.off ()

dist_bray_gorilla <-vegdist(pathway_gorilla, "bray")
adonis(dist_bray_gorilla ~ Group2*Group1, data=metadata_gorilla, permutations=99)
#Df SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)  
#Group2         1 0.0029329 0.00293295  3.8866 0.14474   0.03 *
#  Group1         1 0.0021860 0.00218600  2.8968 0.10788   0.05 *

#---Human only
metadata_human <- read.csv (file = "../Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_human[,2:3]
pathway_human <- Pathway1[24:51,]
pathway_bray_pcoa_human <-pcoa(vegdist(pathway_human, "bray"))
pathway_bray_pcoa_human$values[1:2,]
mds.var.per_g = round(pathway_bray_pcoa_human$values$Eigenvalues/sum(pathway_bray_pcoa_human$values$Eigenvalues)*100, 1)
Pathway_Bray_PCoA_MATRIX_H <- pathway_bray_pcoa_human$vectors[,1:2]
Pathway_Bray_PCoA_MATRIX_H <- data.frame(Pathway_Bray_PCoA_MATRIX_H)
Pathway_Bray_PCoA_MATRIX_H_New <- cbind(Pathway_Bray_PCoA_MATRIX_H, Groups)

Color <- c("brown4", "orangered3")
jpeg("Bray_PCoA_Human_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Pathway_Bray_PCoA_MATRIX_H_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

dist_bray_human <-vegdist(pathway_human, "bray")
adonis(dist_bray_human ~ Group1, data=metadata_human, permutations=99)
#Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)   
#Group1     1  0.016972 0.0169718  9.5008 0.26762   0.01 **
  
#Step 2: ----- 

#Identification of significantly discriminating pathways Human vs Gorilla
# grep -v "NA" KO_ABUNDANCE > KO_ABUNDANCE_removedNA
KO <- read.csv(file = "../KO_ABUNDANCE_removedNA", row.names=1, header = T, sep = "\t")
KO_proportions <- KO/colSums(KO)[col(KO)]
KO1 <- data.frame(t(KO_proportions))
Groups <- metadata$Group
KO_Groups <- cbind (KO1, Groups)

source('~/Work/Metagenomic_Analysis/Final_Analysis/KEGG_Analysis/Odds_Ratio_Calculation.R')
Calculate_Odds_Ratio(KO_Groups)
source('~/Work/Metagenomic_Analysis/Final_Analysis/KEGG_Analysis/WilcoxPairwise.R')
Wilcoxon_pairwise(KO_Groups)

# cut -f1,2 WILCOXON_TEST_TABLE_Human_Gorilla > P_VALUE_TABLE_Human_Gorilla
# cut -f1,5 ODDS_RATIO_TABLE_HumanvsGorilla > FOLD_TABLE_Human_Gorilla
# perl Script.pl FOLD_TABLE_Human_Gorilla P_VALUE_TABLE_Human_Gorilla
# perl Make_Gene_Table.pl Kegg_database/KEGG_2014_map_title.tab Kegg_database/KEGG_2014_ko.list P_VALUE_TABLE_Human_Gorilla

KO_FC_TABLE <- read.csv(file = "FOLD_TABLE_Human_Gorilla.arranged_out", sep = "\t", header = TRUE, row.names = 1)
KO_P_VALUE_TABLE <- read.csv(file = "P_VALUE_TABLE_Human_Gorilla.arranged_out", sep = "\t", header = TRUE, row.names = 1)

GENE_SET_TABLE <- read.csv(file = "GENE_SET_COLLECTION_Human_Gorilla.txt", sep = "\t", header = FALSE)
library(piano)
myGSC <- loadGSC(GENE_SET_TABLE)
GSA_RESULT <- runGSA(KO_P_VALUE_TABLE, KO_FC_TABLE, gsc = myGSC, geneSetStat="reporter", signifMethod="nullDist", adjMethod="fdr", gsSizeLim=c(10,500))
GSAsummaryTable(GSA_RESULT, save = TRUE, file = "GSA_RESULT")
nw <- networkPlot(GSA_RESULT,class = "distinct",direction = "both", label="numbers", significance=0.001)
nw$geneSets

#--- Fold change and p-value calculation for all pathways
metadata <- read.csv (file = "../../Metadata.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata$Group
Groups <- data.frame(Groups)
Pathway1_groups <- cbind(Pathway1, Groups)
source('~/Work/Metagenomic_Analysis/Final_Analysis/KEGG_Analysis/ReporterFeature_Human_Gorilla/Odds_Ratio_Calculation.R')
Calculate_Odds_Ratio(Pathway1_groups)
source('~/Work/Metagenomic_Analysis/Final_Analysis/KEGG_Analysis/ReporterFeature_Human_Gorilla/WilcoxPairwise.R')
Wilcoxon_pairwise(Pathway1_groups)

#----- Identification of significantly discriminating pathways Dry vs Wet
setwd("~/Work/Metagenomic_Analysis/Final_Analysis/KEGG_Analysis/ReporterFeature_Dry_Wet")
metadata_gorilla <- read.csv (file = "../../Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2]
Groups <- data.frame (Groups)
KO_gorilla <- KO1[1:23,]
KO_gorilla_groups <- cbind (KO_gorilla, Groups)

source('~/Work/Metagenomic_Analysis/Final_Analysis/KEGG_Analysis/ReporterFeature_Dry_Wet/Odds_Ratio_Calculation.R')
Calculate_Odds_Ratio(KO_gorilla_groups)
source('~/Work/Metagenomic_Analysis/Final_Analysis/KEGG_Analysis/ReporterFeature_Dry_Wet/WilcoxPairwise.R')
Wilcoxon_pairwise(KO_gorilla_groups)

#cut -f1,2 WILCOXON_TEST_TABLE_Dry_Wet > P_VALUE_TABLE_Dry_Wet
#cut -f1,5 ODDS_RATIO_TABLE_DryvsWet > FOLD_TABLE_Dry_Wet
#perl Script.pl FOLD_TABLE_Dry_Wet P_VALUE_TABLE_Dry_Wet
# perl Make_Gene_Table.pl ../Kegg_database/KEGG_2014_map_title.tab ../Kegg_database/KEGG_2014_ko.list P_VALUE_TABLE_Dry_Wet
# grep -v "NaN" P_VALUE_TABLE_Dry_Wet.arranged_out > P_VALUE_TABLE_Dry_Wet.arranged_NAremoved.out
# grep -v "NaN" FOLD_TABLE_Dry_Wet.arranged_out > FOLD_TABLE_Dry_Wet.arranged_NAremoved.out 

DW_KO_FC_TABLE <- read.csv(file = "FOLD_TABLE_Dry_Wet.arranged_NAremoved.out", sep = "\t", header = TRUE, row.names = 1)
DW_KO_P_VALUE_TABLE <- read.csv(file = "P_VALUE_TABLE_Dry_Wet.arranged_NAremoved.out", sep = "\t", header = TRUE, row.names = 1)

DW_GENE_SET_TABLE <- read.csv(file = "GENE_SET_COLLECTION_Dry_Wet.txt", sep = "\t", header = FALSE)
library(piano)
myGSC_dw <- loadGSC(DW_GENE_SET_TABLE)
DW_GSA_RESULT <- runGSA(DW_KO_P_VALUE_TABLE, DW_KO_FC_TABLE, gsc = myGSC_dw, geneSetStat="reporter", signifMethod="nullDist", adjMethod="fdr", gsSizeLim=c(10,500))
GSAsummaryTable(DW_GSA_RESULT, save = TRUE, file = "DW_GSA_RESULT")
#nw <- networkPlot(DW_GSA_RESULT,class = "distinct",direction = "both", label="numbers", significance=0.001)
#nw$geneSets
source('~/Work/Metagenomic_Analysis/Final_Analysis/KEGG_Analysis/ReporterFeature_Dry_Wet/lefseplotDW.R')

#----- Identification of significantly discriminating pathways BaAka vs Bantu
setwd("~/Work/Metagenomic_Analysis/Final_Analysis/KEGG_Analysis/ReporterFeature_BaAka_Bantu")
metadata_human <- read.csv (file = "../../Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_human[,2]
Groups <- data.frame (Groups)
KO_human <- KO1[24:51,]
KO_human_groups <- cbind (KO_human, Groups)

source('~/Work/Metagenomic_Analysis/Final_Analysis/KEGG_Analysis/ReporterFeature_BaAka_Bantu/Odds_Ratio_Calculation.R')
Calculate_Odds_Ratio(KO_human_groups)
source('~/Work/Metagenomic_Analysis/Final_Analysis/KEGG_Analysis/ReporterFeature_BaAka_Bantu/WilcoxPairwise.R')
Wilcoxon_pairwise(KO_human_groups)

#cut -f1,2 WILCOXON_TEST_TABLE_BaAka_Bantu > P_VALUE_TABLE_BaAka_Bantu
#cut -f1,5 ODDS_RATIO_TABLE_BaAkavsBantu > FOLD_TABLE_BaAka_Bantu
#perl Script.pl FOLD_TABLE_BaAka_Bantu P_VALUE_TABLE_BaAka_Bantu
#perl Make_Gene_Table.pl ../Kegg_database/KEGG_2014_map_title.tab ../Kegg_database/KEGG_2014_ko.list P_VALUE_TABLE_BaAka_Bantu
#grep -v "NaN" FOLD_TABLE_BaAka_Bantu.arranged_out > FOLD_TABLE_BaAka_Bantu.arranged_NAremoved.out
#grep -v "NaN" P_VALUE_TABLE_BaAka_Bantu.arranged_out > P_VALUE_TABLE_BaAka_Bantu.arranged_NAremoved.out

BB_KO_FC_TABLE <- read.csv(file = "FOLD_TABLE_BaAka_Bantu.arranged_NAremoved.out", sep = "\t", header = TRUE, row.names = 1)
BB_KO_P_VALUE_TABLE <- read.csv(file = "P_VALUE_TABLE_BaAka_Bantu.arranged_NAremoved.out", sep = "\t", header = TRUE, row.names = 1)

BB_GENE_SET_TABLE <- read.csv(file = "GENE_SET_COLLECTION_BaAka_Bantu.txt", sep = "\t", header = FALSE)
library(piano)
myGSC_bb <- loadGSC(BB_GENE_SET_TABLE)
BB_GSA_RESULT <- runGSA(BB_KO_P_VALUE_TABLE, BB_KO_FC_TABLE, gsc = myGSC_bb, geneSetStat="reporter", signifMethod="nullDist", adjMethod="fdr", gsSizeLim=c(10,500))
GSAsummaryTable(BB_GSA_RESULT, save = TRUE, file = "BB_GSA_RESULT")

jpeg("BaAkavsBantu.jpg", height = 2, width = 2, units = 'in', res = 600)
networkPlot(BB_GSA_RESULT,class = "distinct",direction = "both", significance=0.01)
dev.off ()
source('~/Work/Metagenomic_Analysis/Final_Analysis/KEGG_Analysis/ReporterFeature_BaAka_Bantu/lefseplotBB.R')
