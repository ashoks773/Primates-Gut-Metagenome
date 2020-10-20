# Author: Ashok Kumar Sharma

# Step 1: Principle coordinate analyis using Pathway_abundance_Table
Pathway <- read.csv(file = "../../../../KEGG_Analysis/PATHWAY_ABUNDANCE_TABLE", row.names=1, header = T, sep = "\t")
Pathway_proportions <- Pathway/colSums(Pathway)[col(Pathway)]
Pathway1 <- data.frame(t(Pathway_proportions))
#write.table (Pathway1, file = "PATHWAY_ABUNDANCE_TABLE_relative_abundances", sep = "\t")
metadata <- read.csv (file = "../../../../Metadata.txt", row.names=1, header = T, sep = "\t")

pathway_bray_pcoa <-pcoa(vegdist(Pathway1, "bray"))
Groups <- metadata[,1:4]
pathway_bray_pcoa$values[1:2,]
mds.var.per1 <- round(pathway_bray_pcoa$values$Rel_corr_eig[1]*100, 2)
mds.var.per2 <- round(pathway_bray_pcoa$values$Rel_corr_eig[2]*100, 2)
#mds.var.per = round(pathway_bray_pcoa$values$Eigenvalues/sum(pathway_bray_pcoa$values$Eigenvalues)*100, 1)
Pathway_Bray_PCoA_MATRIX <- pathway_bray_pcoa$vectors[,1:2]
Pathway_Bray_PCoA_MATRIX <- data.frame(Pathway_Bray_PCoA_MATRIX)
Pathway_Bray_PCoA_MATRIX_New <- cbind(Pathway_Bray_PCoA_MATRIX, Groups)

Color <- c("skyblue4", "wheat4")
Colors <- c("skyblue4", "skyblue", "wheat4", "wheat")

jpeg("Path_Bray_PCoA_Group1_try.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(Pathway_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCoA1 - ", mds.var.per1, "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per2, "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()


distance_bray <-vegdist(Pathway1, "bray")
adonis(distance_bray ~ Group*Group1, data=metadata, permutations=999)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#Group      1  0.050977 0.050977  37.087 0.37857   0.01 **
# Group1     2  0.019079 0.009540   6.940 0.14169   0.01 **

pcoa2 <- Pathway_Bray_PCoA_MATRIX_New[,c(2,4)]
pcoa2_Melted <- melt(pcoa2, id.vars = "Group1")
jpeg("Path_Bray_PCoA2_Distances.jpg", height = 4, width = 2, units = 'in', res = 600)
ggplot(data = pcoa2_Melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCo2") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
kruskalmc(pcoa2$Axis.2 ~ pcoa2$Group1)
#obs.dif critical.dif difference
#BaAka-Bantu 33.285714     14.82396       TRUE
#BaAka-Dry   14.876623     15.80240      FALSE
#BaAka-Wet   23.119048     15.42927       TRUE
#Bantu-Dry   18.409091     15.80240       TRUE
#Bantu-Wet   10.166667     15.42927      FALSE
#Dry-Wet      8.242424     16.37157      FALSE

#--- To check difference of BaAka Bantu with Gorilla - Irrespective of seasons
pcoa2_newgroup <- Pathway_Bray_PCoA_MATRIX_New[,c(2,6)]
pcoa2_newgroup_Melted <- melt(pcoa2_newgroup, id.vars = "GroupNew")
Color <- c("skyblue4", "skyblue", "wheat4")
jpeg("Path Bray_PCoA2_Distances_NewGroup.jpg", height = 4, width = 2, units = 'in', res = 600)
ggplot(data = pcoa2_newgroup_Melted, aes(x=GroupNew, y=value, fill=GroupNew)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCo2") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
kruskalmc(pcoa2_newgroup$Axis.2 ~ pcoa2_newgroup$GroupNew)
#obs.dif critical.dif difference
#BaAka-Bantu   33.28571     13.45140       TRUE
#BaAka-Gorilla 19.17702     12.06395       TRUE
#Bantu-Gorilla 14.10870     12.06395       TRUE

#----- Inter-individual distances
distance_bray <-vegdist(Pathway1, "bray")
distance_bray <- as.matrix(distance_bray)
Bray_dist_column <- melt(distance_bray)
#write.table (Bray_dist_column, file="Bray-distances.txt", sep="\t")
#perl Tag_withLabels.pl ../Metadata.txt Bray-distances.txt

Colors <- c("skyblue4", "skyblue", "wheat4", "wheat")
Bray_dist_groups <- read.csv("../../../../KEGG_Analysis/Bray-distances_formatted.txt.txt", sep = "\t", header=T)
Bray_dist_groups_Melted <- melt(Bray_dist_groups, id.vars = "Groups")
jpeg("Path_Bray_Distances_interindividual.jpg", height = 3.5, width = 2, units = 'in', res = 600)
ggplot(data = Bray_dist_groups_Melted, aes(x=Groups, y=value, fill=Groups)) + geom_boxplot() + ggtitle("") + labs(x="",y="Bray-Curtis distances") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 12, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()

kruskalmc(Bray_dist_groups$Bray_dist ~ Bray_dist_groups$Groups)
#Multiple comparison test after Kruskal-Wallis 
#p.value: 0.05 
#Comparisons
#obs.dif critical.dif difference
#BaAka_BaAka-Bantu_Bantu 154.12245     50.58353       TRUE
#BaAka_BaAka-Dry_Dry     110.94329     57.89367       TRUE
#BaAka_BaAka-Wet_Wet      24.16128     54.96071      FALSE
#Bantu_Bantu-Dry_Dry     265.06574     57.89367       TRUE
#Bantu_Bantu-Wet_Wet     129.96117     54.96071       TRUE
#Dry_Dry-Wet_Wet         135.10457     61.75486       TRUE


#---Gorilla only
metadata_gorilla <- read.csv (file = "../../../../Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
pathway_gorilla <- Pathway1[1:23,]
pathway_bray_pcoa_gorilla <-pcoa(vegdist(pathway_gorilla, "bray"))
pathway_bray_pcoa_gorilla$values[1:2,]
mds.var.per_g1 <- round(pathway_bray_pcoa_gorilla$values$Rel_corr_eig[1]*100, 2)
mds.var.per_g2 <- round(pathway_bray_pcoa_gorilla$values$Rel_corr_eig[2]*100, 2)
#mds.var.per_g = round(pathway_bray_pcoa_gorilla$values$Eigenvalues/sum(pathway_bray_pcoa_gorilla$values$Eigenvalues)*100, 1)
Pathway_Bray_PCoA_MATRIX_G <- pathway_bray_pcoa_gorilla$vectors[,1:2]
Pathway_Bray_PCoA_MATRIX_G <- data.frame(Pathway_Bray_PCoA_MATRIX_G)
Pathway_Bray_PCoA_MATRIX_G_New <- cbind(Pathway_Bray_PCoA_MATRIX_G, Groups)

Color <- c("wheat4", "wheat")
jpeg("Path_Bray_PCoA_Gorilla_Group1.jpg", height = 2.5, width = 3.5, units = 'in', res = 600)
ggplot(Pathway_Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g1, "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g2, "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

dist_bray_gorilla <-vegdist(pathway_gorilla, "bray")
adonis(dist_bray_gorilla ~ Group2*Group1, data=metadata_gorilla, permutations=99)
#Df SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)  
#Group2         1 0.0029329 0.00293295  3.8866 0.14474   0.03 *
# Group1         1 0.0021860 0.00218600  2.8968 0.10788   0.05 *

#Step 2: ----- 

#Identification of significantly discriminating pathways Human vs Gorilla
# grep -v "NA" KO_ABUNDANCE > KO_ABUNDANCE_removedNA
KO <- read.csv(file = "KO_ABUNDANCE_removedNA", row.names=1, header = T, sep = "\t")
KO_proportions <- KO/colSums(KO)[col(KO)]
KO1 <- data.frame(t(KO_proportions))
write.table (KO1, file = "Kegg_KO_relative_abundances", sep = "\t")
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

###################################################
#---- Selected Pathway plots-- For --- Gorilla
###################################################
Pathways_new <- read.csv (file = "../../../../KEGG_Analysis/Pathway_plots/Check_pathways_new2.txt", sep = "\t", row.names = 1, header =T)
Pathways_new <- data.frame(t(Pathways_new))
path2 <- sqrt(Pathways_new)
path3 <- path2*100

metadata_gorilla <- read.csv (file = "../../../../Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
pathway_gorilla <- path3[1:23,]
Color <- c( "wheat4", "wheat")

Pathways_sqrt_per_new <- cbind(pathway_gorilla, Groups)

ValineLeucineIsoLeucineDegradation <- Pathways_sqrt_per_new[,c(8,38)]
ValineLeucineIsoLeucineDegradation_melt <- melt(ValineLeucineIsoLeucineDegradation, id.vars = "Group1")
jpeg("G_ValineLeucineIsoLeucineDegradation.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = ValineLeucineIsoLeucineDegradation_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Valine leucine and isoleucine.degradation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Valine..leucine.and.isoleucine.degradation ~ Pathways_sqrt_per_new$Group1)

Styrene.degradation <- Pathways_sqrt_per_new[,c(14,38)]
Styrene.degradation_melt <- melt(Styrene.degradation, id.vars = "Group1")
jpeg("G_Styrene.degradation.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Styrene.degradation_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Styrene.degradation_melt") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Styrene.degradation ~ Pathways_sqrt_per_new$Group1)

Streptomycin.biosynthesis <- Pathways_sqrt_per_new[,c(34,38)]
Streptomycin.biosynthesis_melt <- melt(Streptomycin.biosynthesis, id.vars = "Group1")
jpeg("G_Streptomycin.biosynthesis_new.jpg", height = 2.5, width = 1, units = 'in', res = 600)
#ggplot(data = Streptomycin.biosynthesis_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Streptomycin.biosynthesis") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
ggplot(data = Streptomycin.biosynthesis, aes(x=Group1, y=Streptomycin.biosynthesis, fill=Group1)) + geom_boxplot() + ggtitle("Streptomycin.biosynthesis") + labs(x="",y="sqrt (Relative abundance)%") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Streptomycin.biosynthesis ~ Pathways_sqrt_per_new$Group1)

Two.component.system <- Pathways_sqrt_per_new[,c(35,38)]
Two.component.system_melt <- melt(Two.component.system, id.vars = "Group1")
jpeg("G_Two.component.system.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Two.component.system_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Two.component.system") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Two.component.system ~ Pathways_sqrt_per_new$Group1)
#W = 28, p-value = 0.01879

ABC.transporters <- Pathways_sqrt_per_new[,c(36,38)]
ABC.transporters_melt <- melt(ABC.transporters, id.vars = "Group1")
jpeg("G_ABC.transporters.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = ABC.transporters_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("ABC.transporters") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$ABC.transporters ~ Pathways_sqrt_per_new$Group1)
#W = 32, p-value = 0.03742

Phosphotransferase.system..PTS <- Pathways_sqrt_per_new[,c(37,38)]
Phosphotransferase.system..PTS_melt <- melt(Phosphotransferase.system..PTS, id.vars = "Group1")
jpeg("G_Phosphotransferase.system..PTS.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Phosphotransferase.system..PTS_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Phosphotransferase.system..PTS") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Phosphotransferase.system..PTS. ~ Pathways_sqrt_per_new$Group1)

###################################################
#---- Selected Pathway plots-- For --- Human
###################################################
Pathways_new <- read.csv (file = "../../../../KEGG_Analysis/Pathway_plots/Check_pathways_new2.txt", sep = "\t", row.names = 1, header =T)
Pathways_new <- data.frame(t(Pathways_new))
path2 <- sqrt(Pathways_new)
path3 <- path2*100

metadata_human <- read.csv (file = "../../../../Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_human[,2:3]
pathway_human <- path3[24:51,]
Color <- c("skyblue4", "skyblue")

Pathways_sqrt_per_new <- cbind(pathway_human, Groups)

ValineLeucineIsoLeucineDegradation <- Pathways_sqrt_per_new[,c(8,38)]
ValineLeucineIsoLeucineDegradation_melt <- melt(ValineLeucineIsoLeucineDegradation, id.vars = "Group1")
jpeg("H_ValineLeucineIsoLeucineDegradation.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = ValineLeucineIsoLeucineDegradation_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Valine leucine and isoleucine.degradation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Valine..leucine.and.isoleucine.degradation ~ Pathways_sqrt_per_new$Group1)
#W = 0, p-value = 4.985e-08

Styrene.degradation <- Pathways_sqrt_per_new[,c(14,38)]
Styrene.degradation_melt <- melt(Styrene.degradation, id.vars = "Group1")
jpeg("H_Styrene.degradation.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Styrene.degradation_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Styrene.degradation_melt") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Styrene.degradation ~ Pathways_sqrt_per_new$Group1)
#W = 34, p-value = 0.002486

Streptomycin.biosynthesis <- Pathways_sqrt_per_new[,c(34,38)]
Streptomycin.biosynthesis_melt <- melt(Streptomycin.biosynthesis, id.vars = "Group1")
jpeg("H_Streptomycin.biosynthesis_new.jpg", height = 2.5, width = 1, units = 'in', res = 600)
#ggplot(data = Streptomycin.biosynthesis_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Streptomycin.biosynthesis") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
ggplot(data = Streptomycin.biosynthesis, aes(x=Group1, y=Streptomycin.biosynthesis, fill=Group1)) + geom_boxplot() + ggtitle("Streptomycin.biosynthesis") + labs(x="",y="sqrt (Relative abundance)%") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Streptomycin.biosynthesis ~ Pathways_sqrt_per_new$Group1)
#W = 179, p-value = 5.973e-05

Two.component.system <- Pathways_sqrt_per_new[,c(35,38)]
Two.component.system_melt <- melt(Two.component.system, id.vars = "Group1")
jpeg("H_Two.component.system.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Two.component.system_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Two.component.system") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Two.component.system ~ Pathways_sqrt_per_new$Group1)
#W = 40, p-value = 0.006743

ABC.transporters <- Pathways_sqrt_per_new[,c(36,38)]
ABC.transporters_melt <- melt(ABC.transporters, id.vars = "Group1")
jpeg("H_ABC.transporters.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = ABC.transporters_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("ABC.transporters") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$ABC.transporters ~ Pathways_sqrt_per_new$Group1)
#W = 12, p-value = 1.356e-05

Phosphotransferase.system..PTS <- Pathways_sqrt_per_new[,c(37,38)]
Phosphotransferase.system..PTS_melt <- melt(Phosphotransferase.system..PTS, id.vars = "Group1")
jpeg("H_Phosphotransferase.system..PTS.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Phosphotransferase.system..PTS_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Phosphotransferase.system..PTS") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Phosphotransferase.system..PTS. ~ Pathways_sqrt_per_new$Group1)
#W = 45, p-value = 0.01411