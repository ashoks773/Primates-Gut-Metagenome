# Author: Ashok Kumar Sharma

# Step 1: Principle coordinate analyis using Cazy Family table
Fam <- read.csv(file = "Fam_ABUNDANCE", row.names=1, header = T, sep = "\t")
Fam_proportions <- Fam/colSums(Fam)[col(Fam)]
Fam1 <- data.frame(t(Fam_proportions))
#write.table (Fam1, file = "Fam_ABUNDANCE_relative_abundances", sep = "\t")
metadata <- read.csv (file = "Final_Data/Metadata.txt", row.names=1, header = T, sep = "\t")

#--- total
Fam_bray_pcoa <-pcoa(vegdist(Fam1, "bray"))
Groups <- metadata[,1:4]
Fam_bray_pcoa$values[1:2,]
mds.var.per1 <- round(Fam_bray_pcoa$values$Rel_corr_eig[1]*100, 2)
mds.var.per2 <- round(Fam_bray_pcoa$values$Rel_corr_eig[2]*100, 2)
#mds.var.per = round(Fam_bray_pcoa$values$Eigenvalues/sum(Fam_bray_pcoa$values$Eigenvalues)*100, 1)
Fam_Bray_PCoA_MATRIX <- Fam_bray_pcoa$vectors[,1:2]
Fam_Bray_PCoA_MATRIX <- data.frame(Fam_Bray_PCoA_MATRIX)
Fam_Bray_PCoA_MATRIX_New <- cbind(Fam_Bray_PCoA_MATRIX, Groups)

Color <- c("skyblue4", "wheat4")
Colors <- c("skyblue4", "skyblue", "wheat4", "wheat")


jpeg("Cazy_Bray_PCoA_Group1.jpg", height = 2.5, width = 3.5, units = 'in', res = 600)
ggplot(Fam_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCoA1 - ", mds.var.per1, "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per2, "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

distance_bray <-vegdist(Fam1, "bray")
adonis(distance_bray ~ Group*Group1, data=metadata, permutations=99)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#Group      1   0.21125 0.211253 22.1124 0.25330   0.01 **
#Group1     2   0.17373 0.086863  9.0922 0.20831   0.01 **

#----- Plot Axis2
pcoa2 <- Fam_Bray_PCoA_MATRIX_New[,c(2,4)]
pcoa2_Melted <- melt(pcoa2, id.vars = "Group1")
jpeg("Cazy_Bray_PCoA2_Distances.jpg", height = 3.5, width = 2, units = 'in', res = 600)
ggplot(data = pcoa2_Melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCo2") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
kruskalmc(pcoa2$Axis.2 ~ pcoa2$Group1)
#obs.dif critical.dif difference
#BaAka-Bantu 17.357143     14.82396       TRUE
#BaAka-Dry   15.220779     15.80240      FALSE
#BaAka-Wet    9.940476     15.42927      FALSE
#Bantu-Dry   32.577922     15.80240       TRUE
#Bantu-Wet   27.297619     15.42927       TRUE
#Dry-Wet      5.280303     16.37157      FALSE

#--- To check difference of BaAka Bantu with Gorilla - Irrespective of seasons
#pcoa2_newgroup <- Fam_Bray_PCoA_MATRIX_New[,c(2,6)]
#pcoa2_newgroup_Melted <- melt(pcoa2_newgroup, id.vars = "GroupNew")
#Color <- c("brown4", "orangered3", "darkgreen")
#jpeg("Bray_PCoA2_Distances_NewGroup.jpg", height = 4, width = 2, units = 'in', res = 600)
#ggplot(data = pcoa2_newgroup_Melted, aes(x=GroupNew, y=value, fill=GroupNew)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCo2") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
##dev.off ()
#kruskalmc(pcoa2_newgroup$Axis.2 ~ pcoa2_newgroup$GroupNew)

#----- Inter-individual distances
distance_bray <-vegdist(Fam1, "bray")
distance_bray <- as.matrix(distance_bray)
Bray_dist_column <- melt(distance_bray)
#write.table (Bray_dist_column, file="Bray-distances.txt", sep="\t")
#perl Tag_withLabels.pl ../Metadata.txt Bray-distances.txt

Colors <- c("skyblue4", "skyblue", "wheat4", "wheat")
Bray_dist_groups <- read.csv("Bray-distances_filtered.txt.txt", sep = "\t", header=T)
Bray_dist_groups_Melted <- melt(Bray_dist_groups, id.vars = "Groups")
jpeg("CAzy_Bray_Distances_interindividual.jpg", height = 3.5, width = 2, units = 'in', res = 600)
ggplot(data = Bray_dist_groups_Melted, aes(x=Groups, y=value, fill=Groups)) + geom_boxplot() + ggtitle("") + labs(x="",y="Bray-Curtis distances") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 10, colour = "black"), axis.text.y = element_text(size = 10, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()

kruskalmc(Bray_dist_groups$Bray.dist ~ Bray_dist_groups$Groups)
#Multiple comparison test after Kruskal-Wallis 
#p.value: 0.05 
#Comparisons
#obs.dif critical.dif difference
#BaAka_BaAka-Bantu_Bantu 205.20408     50.58353       TRUE
#BaAka_BaAka-Dry_Dry     136.95712     57.89367       TRUE
#BaAka_BaAka-Wet_Wet      24.16582     54.96071      FALSE
#Bantu_Bantu-Dry_Dry     342.16120     57.89367       TRUE
#Bantu_Bantu-Wet_Wet     181.03827     54.96071       TRUE
#Dry_Dry-Wet_Wet         161.12293     61.75486       TRUE

#------- Cazy Broad Category Relative abundances
cazy_broad_rel_abun <- read.csv("Cazy_BroadCat_Rel_Abun.txt", sep = "\t", header=T, row.names = 1)
cazy_broad_rel_abun <- data.frame(t(cazy_broad_rel_abun))
cazy_meta <- cbind(cazy_broad_rel_abun, metadata)
Cazy_Group <- cazy_meta[,c(1:6,8)]
#--- Plotted in Excel file name CazY_Fam_Rel_Abundance.xlsx

#---- Gorilla
Color <- c("wheat4", "wheat")
metadata_gorilla <- read.csv (file = "Final_Data/Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
fam_gorilla <- Fam1[1:23,]
fam_bray_pcoa_gorilla <-pcoa(vegdist(fam_gorilla, "bray"))
fam_bray_pcoa_gorilla$values[1:2,]
mds.var.per1 <- round(fam_bray_pcoa_gorilla$values$Rel_corr_eig[1]*100, 2)
mds.var.per2 <- round(fam_bray_pcoa_gorilla$values$Rel_corr_eig[2]*100, 2)
fam_Bray_PCoA_MATRIX_G <- fam_bray_pcoa_gorilla$vectors[,1:2]
fam_Bray_PCoA_MATRIX_G <- data.frame(fam_Bray_PCoA_MATRIX_G)
fam_Bray_PCoA_MATRIX_G_New <- cbind(fam_Bray_PCoA_MATRIX_G, Groups)
jpeg("Cazy_PCoA_Gorilla_Group1.jpg", height = 2.5, width = 3.2, units = 'in', res = 600)
ggplot(fam_Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per1, "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per2, "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

dist_bray_gorilla <-vegdist(fam_gorilla, "bray")
adonis(dist_bray_gorilla ~ Group2*Group1, data=metadata_gorilla, permutations=99)
#Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)   
#Group2         1  0.014959 0.0149586  3.2089 0.11038   0.01 **
# Group1         1  0.018447 0.0184471  3.9572 0.13612   0.01 **
# Group2:Group1  1  0.013546 0.0135456  2.9058 0.09995   0.01 **



# Step 2: Statistical Analysis 
#-Human vs Gorilla
setwd("~/Work/Metagenomic_Analysis/Final_Analysis/CaZY_analysis/Human_vs_Gorilla")
Groups <- metadata$Group
Groups <- data.frame(Groups)
Fam1_groups <- cbind (Fam1, Groups)
source('~/Work/Metagenomic_Analysis/Final_Analysis/CaZY_analysis/Human_vs_Gorilla/WilcoxPairwise.R')
Wilcoxon_pairwise(Fam1_groups)
source('~/Work/Metagenomic_Analysis/Final_Analysis/CaZY_analysis/Human_vs_Gorilla/Odds_Ratio_Calculation.R')
Calculate_Odds_Ratio(Fam1_groups)

#---- Enrichment Analysis
# grep -v "Inf" Qval_foldchange_Table > Qval_foldchange_Table_filtered
Fam_foldchange_pval <- read.csv (file = "Qval_foldchange_Table_filtered", sep = "\t", header =TRUE)

colnames(Fam_foldchange_pval) <- c("CazyFam", "Qval", "Fold")

Fam_foldchange_pval["group"] <- "NotSignificant"
Fam_foldchange_pval[which(Fam_foldchange_pval['Qval'] < 0.05 & abs(Fam_foldchange_pval['Fold']) < 2 ),"group"] <- "Significant"
Fam_foldchange_pval[which(Fam_foldchange_pval['Qval'] > 0.05 & abs(Fam_foldchange_pval['Fold']) > 2 ),"group"] <- "FoldChange"
Fam_foldchange_pval[which(Fam_foldchange_pval['Qval'] < 0.05 & abs(Fam_foldchange_pval['Fold']) > 2 ),"group"] <- "Significant&FoldChange"

png(filename="Odds_ratio.png", height = 10, width = 12, units = 'in', res = 600)
ggplot(Fam_foldchange_pval, aes(x=Fold, y= -log(Qval, 10))) + geom_point(data = subset (Fam_foldchange_pval, (-log(Qval, 10) < 50) & (-log(Qval, 10) > 0.001)), aes(fill = group), size = 3.0, pch = 21, colour = "black") + theme(text = element_text(family = "Times New Roman")) + theme(axis.line.x = element_line(colour = "black", size= 0.8), axis.line.y = element_line(colour = "black", size = 0.8), axis.text = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 15, face = "bold"), axis.title = element_text(size = 15, face = "bold", family = "Times New Roman")) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted") + theme_classic() + geom_text_repel(data = subset(Fam_foldchange_pval, ((-log(Qval,10) < 50)& (Fold > 2))| ((-log(Qval,10) < 50) & (Fold < -2))), aes(label = CazyFam))
dev.off ()

#-Dry vs Wet
setwd("~/Work/Metagenomic_Analysis/Final_Analysis/CaZY_analysis/Dry_vs_Wet")
Groups <- metadata_gorilla$Group1
fam_gorilla <- Fam1[1:23,]
Groups <- data.frame(Groups)
fam_gorilla_groups <- cbind (fam_gorilla, Groups)
source('~/Work/Metagenomic_Analysis/Final_Analysis/CaZY_analysis/Dry_vs_Wet/Odds_Ratio_Calculation.R')
Calculate_Odds_Ratio(fam_gorilla_groups)
source('~/Work/Metagenomic_Analysis/Final_Analysis/CaZY_analysis/Dry_vs_Wet/WilcoxPairwise.R')
Wilcoxon_pairwise(fam_gorilla_groups)

#---- Enrichment Analysis
#grep -v "NaN" Qval_foldchange_Table > Qval_foldchange_Table_filtered
#grep -v "Inf" Qval_foldchange_Table_filtered > Qval_foldchange_Table_filtered1
#mv Qval_foldchange_Table_filtered1 Qval_foldchange_Table_filtered
Fam_foldchange_pval <- read.csv (file = "Qval_foldchange_Table_filtered", sep = "\t", header =TRUE)

colnames(Fam_foldchange_pval) <- c("CazyFam", "Qval", "Fold")

Fam_foldchange_pval["group"] <- "NotSignificant"
Fam_foldchange_pval[which(Fam_foldchange_pval['Qval'] < 0.05 & abs(Fam_foldchange_pval['Fold']) < 2 ),"group"] <- "Significant"
Fam_foldchange_pval[which(Fam_foldchange_pval['Qval'] > 0.05 & abs(Fam_foldchange_pval['Fold']) > 2 ),"group"] <- "FoldChange"
Fam_foldchange_pval[which(Fam_foldchange_pval['Qval'] < 0.05 & abs(Fam_foldchange_pval['Fold']) > 2 ),"group"] <- "Significant&FoldChange"
# -log(0.05, 10) this is 1.301 (If this is greater than 1.3 Means it is significant)
# -log(0.01, 10) this is 2
png(filename="Dry_vs_Wet_Odds_ratio.png", height = 10, width = 12, units = 'in', res = 600)
ggplot(Fam_foldchange_pval, aes(x=Fold, y= -log(Qval, 10))) + geom_point(data = subset (Fam_foldchange_pval, (-log(Qval, 10) < 50) & (-log(Qval, 10) > 0.001)), aes(fill = group), size = 3.0, pch = 21, colour = "black") + theme(text = element_text(family = "Times New Roman")) + theme(axis.line.x = element_line(colour = "black", size= 0.8), axis.line.y = element_line(colour = "black", size = 0.8), axis.text = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 15, face = "bold"), axis.title = element_text(size = 15, face = "bold", family = "Times New Roman")) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted") + theme_classic() + geom_text_repel(data = subset(Fam_foldchange_pval, ((-log(Qval,10) > 1.30)& (Fold > 1.5))| ((-log(Qval,10) > 1.30) & (Fold < -2))), aes(label = CazyFam))
dev.off ()

#-Baka vs Bantu
setwd("~/Work/Metagenomic_Analysis/Final_Analysis/CaZY_analysis/BaAka_vs_Bantu")
Groups <- metadata_human$Group1
fam_human <- Fam1[24:51,]
Groups <- data.frame(Groups)
fam_human_groups <- cbind (fam_human, Groups)
source('~/Work/Metagenomic_Analysis/Final_Analysis/CaZY_analysis/BaAka_vs_Bantu/Odds_Ratio_Calculation.R')
Calculate_Odds_Ratio(fam_human_groups)
source('~/Work/Metagenomic_Analysis/Final_Analysis/CaZY_analysis/BaAka_vs_Bantu/WilcoxPairwise.R')
Wilcoxon_pairwise(fam_human_groups)

#---- Enrichment Analysis
#grep -v "NaN" Qval_foldchange_Table > Qval_foldchange_Table_filtered
#grep -v "Inf" Qval_foldchange_Table_filtered > Qval_foldchange_Table_filtered1
#mv Qval_foldchange_Table_filtered1 Qval_foldchange_Table_filtered
Fam_foldchange_pval <- read.csv (file = "Qval_foldchange_Table_filtered", sep = "\t", header =TRUE)

colnames(Fam_foldchange_pval) <- c("CazyFam", "Qval", "Fold")

Fam_foldchange_pval["group"] <- "NotSignificant"
Fam_foldchange_pval[which(Fam_foldchange_pval['Qval'] < 0.05 & abs(Fam_foldchange_pval['Fold']) < 0.5 ),"group"] <- "Significant"
Fam_foldchange_pval[which(Fam_foldchange_pval['Qval'] > 0.05 & abs(Fam_foldchange_pval['Fold']) > 0.5 ),"group"] <- "FoldChange"
Fam_foldchange_pval[which(Fam_foldchange_pval['Qval'] < 0.05 & abs(Fam_foldchange_pval['Fold']) > 0.5 ),"group"] <- "Significant&FoldChange"
# -log(0.05, 10) this is 1.301 (If this is greater than 1.3 Means it is significant)
# -log(0.01, 10) this is 2
png(filename="BaAka_vs_Bantu_Odds_ratio.png", height = 10, width = 12, units = 'in', res = 600)
ggplot(Fam_foldchange_pval, aes(x=Fold, y= -log(Qval, 10))) + geom_point(data = subset (Fam_foldchange_pval, (-log(Qval, 10) <50) & (-log(Qval, 10) > 0.001)), aes(fill = group), size = 3.0, pch = 21, colour = "black") + theme(text = element_text(family = "Times New Roman")) + theme(axis.line.x = element_line(colour = "black", size= 0.8), axis.line.y = element_line(colour = "black", size = 0.8), axis.text = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12, face = "bold"), legend.title = element_text(size = 15, face = "bold"), axis.title = element_text(size = 15, face = "bold", family = "Times New Roman")) + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted") + theme_classic() + geom_text_repel(data = subset(Fam_foldchange_pval, ((-log(Qval,10) > 1.30)& (Fold > 0.5))| ((-log(Qval,10) > 1.30) & (Fold < -0.5))), aes(label = CazyFam))
dev.off ()


#---HeatMap
write.table(Selected_cazy_Groups, file="Selected_cazy_Groups.txt", sep = "\t")
#--- From this file GH27 removed and samples formated by group in excel
Selected_cazy_new <- read.csv(file="Selected_cazy_Groups_format.txt", sep = "\t", row.names = 1, header =T)
categories <- Selected_cazy_new[,20:21]
Selected_cazy_pred <- Selected_cazy_new[,1:19]
rownames(categories) <- rownames(Selected_cazy_pred)
annot.color.col <- list('Group'=c('wheat4', 'skyblue4'),'Group1'=c("skyblue4", "skyblue", "wheat4", "wheat"))
aheatmap(sqrt(Selected_cazy_pred)*100, color = "-RdYlBu2:100", breaks = 0,  main = "Association", distfun = "spearman", hclustfun = "complete", fontsize=10,  filename="CaZy_Finally_Used.png", scale = "row", Rowv = NA, Colv = TRUE, border_color = NA, legend = TRUE, cexRow = .7, cexCol = .7, annRow = categories, annCol = NA, annColors = annot.color.col, annLegend = TRUE)


#--- Individual BoxPlots For Human
metadata_human <- read.csv (file = "Final_Data/Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Colors <- c("skyblue4", "skyblue")

#Digestion of animal Polysaccharides (GH13_32, GH13_23, GH13_41, GH13_16, PL13, PL21_1)
# Digestion of Starch, Glycogen and Mucopolysaccharide
Animal_poly <- read.csv(file="Animal_poly", sep = "\t", row.names = 1, header =T)
Animal <- colSums(Animal_poly)
Animal <- data.frame(Animal)
Animal <- Animal[24:51,]
Animal_Group <- cbind(Animal, metadata_human)
jpeg("Animal_h.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Animal_Group, aes(x=Group1, y=sqrt(Animal)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of animal polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size =10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Animal ~ Group1, data = Animal_Group)
#W = 106, p-value = 0.7345

#Digestion of plant Polysaccharides (GH27, GH43_8, PL4_2)
Plant_poly <- read.csv(file="plant_poly", sep = "\t", row.names = 1, header =T)
Plant <- colSums(Plant_poly)
Plant <- data.frame(Plant)
Plant <- Plant[24:51,]
Plant_Group <- cbind(Plant, metadata_human)
jpeg("Plant_h.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Plant_Group, aes(x=Group1, y=sqrt(Plant)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of plant polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Plant ~ Group1, data = Plant_Group)
#W = 183, p-value = 1.86e-05

#Digestion of fungal Polysaccharides (GH113)
#Digestion of Mannans, Galactomannans and Glucomannans
Fungal_poly <- read.csv(file="Fungal_poly", sep = "\t", row.names = 1, header =T)
Fungal <- colSums(Fungal_poly)
Fungal <- data.frame(Fungal)
Fungal <- Fungal[24:51,]
Fungal_Group <- cbind(Fungal, metadata_human)
jpeg("Fungal_h.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Fungal_Group, aes(x=Group1, y=sqrt(Fungal)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of fungal polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Fungal ~ Group1, data = Fungal_Group)
#W = 99, p-value = 0.982

#Digestion of Algal Polysaccharides (GH96)
Algal_poly <- read.csv(file="Algal_poly", sep = "\t", row.names = 1, header =T)
Algal <- colSums(Algal_poly)
Algal <- data.frame(Algal)
Algal <- Algal[24:51,]
Algal_Group <- cbind(Algal, metadata_human)
jpeg("Algal_h.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Algal_Group, aes(x=Group1, y=sqrt(Algal)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of algal polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Algal ~ Group1, data = Algal_Group)
#W = 112, p-value = 0.1649

#-- Lignin degradation (AA1, AA3_2, AA7)
Lignin_degradation <- read.csv(file="Lignin_degradation", sep = "\t", row.names = 1, header =T)
Lignin <- colSums(Lignin_degradation)
Lignin <- data.frame(Lignin)
Lignin <- Lignin[24:51,]
Lignin_Group <- cbind(Lignin, metadata_human)
jpeg("Lignin_h.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Lignin_Group, aes(x=Group1, y=sqrt(Lignin)*100, fill=Group1)) + geom_boxplot() + ggtitle("Lignin degradation") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size =10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Lignin ~ Group1, data = Lignin_Group)
#W = 132, p-value = 0.1227

#--#Biosynthesis of Polysaccharides (GT24, GT88, GT39, GT64, GT10, GT21)
Biosynthesis_poly <- read.csv(file="Biosynthesis_poly", sep = "\t", row.names = 1, header =T)
Biosynthesis <- colSums(Biosynthesis_poly)
Biosynthesis <- data.frame(Biosynthesis)
Biosynthesis <- Biosynthesis[24:51,]
Biosynthesis_Group <- cbind(Biosynthesis, metadata_human)
jpeg("Biosynthesis_h.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Biosynthesis_Group, aes(x=Group1, y=sqrt(Biosynthesis)*100, fill=Group1)) + geom_boxplot() + ggtitle("Biosynthesis of polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Biosynthesis ~ Group1, data = Biosynthesis_Group) 
#W = 65, p-value = 0.1371


####--- For Gorilla
metadata_gorilla <- read.csv (file = "Final_Data/Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Colors <- c("wheat4", "wheat")

#Digestion of animal Polysaccharides (GH13_32, GH13_23, GH13_41, GH13_16, PL13, PL21_1)
Animal_poly <- read.csv(file="Animal_poly", sep = "\t", row.names = 1, header =T)
Animal <- colSums(Animal_poly)
Animal <- data.frame(Animal)
Animal <- Animal[1:23,]
Animal_Group <- cbind(Animal, metadata_gorilla)
jpeg("Animal_g.jpg", height = 2.5, width = 1.2, units = 'in', res = 600)
ggplot(data = Animal_Group, aes(x=Group1, y=sqrt(Animal)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of animal polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold")) 
dev.off ()
wilcox.test(Animal ~ Group1, data = Animal_Group)
#W = 131, p-value = 2.958e-06

#Digestion of plant Polysaccharides (GH27, GH43_8, PL4_2)
Plant_poly <- read.csv(file="plant_poly", sep = "\t", row.names = 1, header =T)
Plant <- colSums(Plant_poly)
Plant <- data.frame(Plant)
Plant <- Plant[1:23,]
Plant_Group <- cbind(Plant, metadata_gorilla)
jpeg("Plant_g.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Plant_Group, aes(x=Group1, y=sqrt(Plant)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of plant polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Plant ~ Group1, data = Plant_Group)
#W = 99, p-value = 0.04388

#Digestion of fungal Polysaccharides (GH113)
Fungal_poly <- read.csv(file="Fungal_poly", sep = "\t", row.names = 1, header =T)
Fungal <- colSums(Fungal_poly)
Fungal <- data.frame(Fungal)
Fungal <- Fungal[1:23,]
Fungal_Group <- cbind(Fungal, metadata_gorilla)
jpeg("Fungal_g.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Fungal_Group, aes(x=Group1, y=sqrt(Fungal)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of fungal polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Fungal ~ Group1, data = Fungal_Group)
#W = 16, p-value = 0.001294

#Digestion of Algal Polysaccharides (GH96)
Algal_poly <- read.csv(file="Algal_poly", sep = "\t", row.names = 1, header =T)
Algal <- colSums(Algal_poly)
Algal <- data.frame(Algal)
Algal <- Algal[1:23,]
Algal_Group <- cbind(Algal, metadata_gorilla)
jpeg("Algal_g.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Algal_Group, aes(x=Group1, y=sqrt(Algal)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of algal polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Algal ~ Group1, data = Algal_Group)
#W = 122, p-value = 0.0002056

#-- Lignin degradation (AA1, AA3_2, AA7)
Lignin_degradation <- read.csv(file="Lignin_degradation", sep = "\t", row.names = 1, header =T)
Lignin <- colSums(Lignin_degradation)
Lignin <- data.frame(Lignin)
Lignin <- Lignin[1:23,]
Lignin_Group <- cbind(Lignin, metadata_gorilla)
jpeg("Lignin_g.jpg", height = 2.5, width = 1.2, units = 'in', res = 600)
ggplot(data = Lignin_Group, aes(x=Group1, y=sqrt(Lignin)*100, fill=Group1)) + geom_boxplot() + ggtitle("Lignin degradation") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Lignin ~ Group1, data = Lignin_Group)
#W = 103.5, p-value = 0.02161

#--#Biosynthesis of Polysaccharides (GT24, GT88, GT39, GT64, GT10, GT21)
Biosynthesis_poly <- read.csv(file="Biosynthesis_poly", sep = "\t", row.names = 1, header =T)
Biosynthesis <- colSums(Biosynthesis_poly)
Biosynthesis <- data.frame(Biosynthesis)
Biosynthesis <- Biosynthesis[1:23,]
Biosynthesis_Group <- cbind(Biosynthesis, metadata_gorilla)
jpeg("Biosynthesis_g.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Biosynthesis_Group, aes(x=Group1, y=sqrt(Biosynthesis)*100, fill=Group1)) + geom_boxplot() + ggtitle("Biosynthesis of polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold")) 
dev.off ()
wilcox.test(Biosynthesis ~ Group1, data = Biosynthesis_Group) 
#W = 0, p-value = 1.479e-06