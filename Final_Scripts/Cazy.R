# Step 1: Principle coordinate analyis using Cazy Family table
Fam <- read.csv(file = "Fam_ABUNDANCE", row.names=1, header = T, sep = "\t")
Fam_proportions <- Fam/colSums(Fam)[col(Fam)]
Fam1 <- data.frame(t(Fam_proportions))
write.table (Fam1, file = "Fam_ABUNDANCE_relative_abundances", sep = "\t")
metadata <- read.csv (file = "../Metadata.txt", row.names=1, header = T, sep = "\t")

#--- total
Fam_bray_pcoa <-pcoa(vegdist(Fam1, "bray"))
Groups <- metadata[,1:3]
Fam_bray_pcoa$values[1:2,]
mds.var.per = round(Fam_bray_pcoa$values$Eigenvalues/sum(Fam_bray_pcoa$values$Eigenvalues)*100, 1)
Fam_Bray_PCoA_MATRIX <- Fam_bray_pcoa$vectors[,1:2]
Fam_Bray_PCoA_MATRIX <- data.frame(Fam_Bray_PCoA_MATRIX)
Fam_Bray_PCoA_MATRIX_New <- cbind(Fam_Bray_PCoA_MATRIX, Groups)
Color <- c("darkgreen", "darkred")
Colors <- c("brown4", "orangered3", "green4", "limegreen")

jpeg("Cazy_Bray_PCoA_Group.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Fam_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
jpeg("Cazy_Bray_PCoA_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Fam_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCoA1 - ", mds.var.per[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

distance_bray <-vegdist(Fam1, "bray")
adonis(distance_bray ~ Group*Group1, data=metadata, permutations=99)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#Group      1   0.21125 0.211253 22.1124 0.25330   0.01 **
#Group1     2   0.17373 0.086863  9.0922 0.20831   0.01 **

#----- Inter-individual distances
distance_bray <-vegdist(Fam1, "bray")
distance_bray <- as.matrix(distance_bray)
Bray_dist_column <- melt(distance_bray)
write.table (Bray_dist_column, file="Bray-distances.txt", sep="\t")
#perl Tag_withLabels.pl ../Metadata.txt Bray-distances.txt

Colors <- c("brown4", "orangered3", "green4", "limegreen")
Bray_dist_groups <- read.csv("Bray-distances_filtered.txt.txt", sep = "\t", header=T)
Bray_dist_groups_Melted <- melt(Bray_dist_groups, id.vars = "Groups")
jpeg("Bray_Distances_interindividual.jpg", height = 4, width = 2, units = 'in', res = 600)
ggplot(data = Bray_dist_groups_Melted, aes(x=Groups, y=value, fill=Groups)) + geom_boxplot() + ggtitle("") + labs(x="",y="Bray-Curtis distances") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()

kruskalmc(Bray_dist_groups$Bray.dist ~ Bray_dist_groups$Groups)
#Multiple comparison test after Kruskal-Wallis 
#p.value: 0.05 
#Comparisons
#obs.dif critical.dif difference
#BaAka_BaAka-Bantu_Bantu 205.20408     50.58353       TRUE
#BaAka_BaAka-Dry_Dry     136.95712     57.89367       TRUE
B#aAka_BaAka-Wet_Wet      24.16582     54.96071      FALSE
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
Color <- c("green4", "limegreen")
metadata_gorilla <- read.csv (file = "../Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
fam_gorilla <- Fam1[1:23,]
fam_bray_pcoa_gorilla <-pcoa(vegdist(fam_gorilla, "bray"))
fam_bray_pcoa_gorilla$values[1:2,]
mds.var.per_g = round(fam_bray_pcoa_gorilla$values$Eigenvalues/sum(fam_bray_pcoa_gorilla$values$Eigenvalues)*100, 1)
fam_Bray_PCoA_MATRIX_G <- fam_bray_pcoa_gorilla$vectors[,1:2]
fam_Bray_PCoA_MATRIX_G <- data.frame(fam_Bray_PCoA_MATRIX_G)
fam_Bray_PCoA_MATRIX_G_New <- cbind(fam_Bray_PCoA_MATRIX_G, Groups)
jpeg("Cazy_PCoA_Gorilla_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(fam_Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
jpeg("Cazy_PCoA_Gorilla_Group2.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(fam_Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1, shape=Group2)) + geom_point(size=3) + geom_point(aes(shape=Group2)) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted") + geom_point(shape=18)
dev.off ()

dist_bray_gorilla <-vegdist(fam_gorilla, "bray")
adonis(dist_bray_gorilla ~ Group2*Group1, data=metadata_gorilla, permutations=99)
#Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)   
#Group2         1  0.014959 0.0149586  3.2089 0.11038   0.01 **
# Group1         1  0.018447 0.0184471  3.9572 0.13612   0.01 **
# Group2:Group1  1  0.013546 0.0135456  2.9058 0.09995   0.01 **

#---- Human
metadata_human <- read.csv (file = "../Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_human[,2:3]
fam_human <- Fam1[24:51,]
fam_bray_pcoa_human <-pcoa(vegdist(fam_human, "bray"))
fam_bray_pcoa_human$values[1:2,]
mds.var.per_g = round(fam_bray_pcoa_human$values$Eigenvalues/sum(fam_bray_pcoa_human$values$Eigenvalues)*100, 1)
fam_Bray_PCoA_MATRIX_H <- fam_bray_pcoa_human$vectors[,1:2]
fam_Bray_PCoA_MATRIX_H <- data.frame(fam_Bray_PCoA_MATRIX_H)
fam_Bray_PCoA_MATRIX_H_New <- cbind(fam_Bray_PCoA_MATRIX_H, Groups)

Color <- c("brown4", "orangered3")
jpeg("Cazy_PCoA_Human_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(fam_Bray_PCoA_MATRIX_H_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

dist_bray_human <-vegdist(fam_human, "bray")
adonis(dist_bray_human ~ Group1, data=metadata_human, permutations=99)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#Group1     1   0.15486 0.154857  12.114 0.31784   0.01 **


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


#---- Broad CaZy class
metadata <- read.csv (file = "../Metadata.txt", row.names=1, header = T, sep = "\t")
Group1 <- metadata$Group1
Group1 <- data.frame(Group1)

Color <- c("darkgreen", "darkred")
Colors <- c("brown4", "orangered3", "green4", "limegreen")

#------------------
#Digesion of starch and glycogen?
#Digestion of animal Polysaccharides (GH13_32, GH13_23, GH13_41, GH13_16, PL13, PL21_1)
Animal_poly <- read.csv(file="Animal_poly", sep = "\t", row.names = 1, header =T)
Animal <- colSums(Animal_poly)
Animal <- data.frame(Animal)
Animal_Group <- cbind(Animal, metadata)
jpeg("Animal.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = Animal_Group, aes(x=Group1, y=sqrt(Animal)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of animal polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
kruskal.test(Animal ~ Group1, data = Animal_Group)
#Kruskal-Wallis chi-squared = 27.92, df = 3, p-value = 3.774e-06

#Digestion of plant Polysaccharides (GH27, GH43_8, PL4_2)
Plant_poly <- read.csv(file="plant_poly", sep = "\t", row.names = 1, header =T)
Plant <- colSums(Plant_poly)
Plant <- data.frame(Plant)
Plant_Group <- cbind(Plant, metadata)
jpeg("Plant.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = Plant_Group, aes(x=Group1, y=sqrt(Plant)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of plant polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
kruskal.test(Plant ~ Group1, data = Plant_Group)
#Kruskal-Wallis chi-squared = 32.638, df = 3, p-value = 3.839e-07

#Digestion of fungal Polysaccharides (GH113)
Fungal_poly <- read.csv(file="Fungal_poly", sep = "\t", row.names = 1, header =T)
Fungal <- colSums(Fungal_poly)
Fungal <- data.frame(Fungal)
Fungal_Group <- cbind(Fungal, metadata)
jpeg("Fungal.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = Fungal_Group, aes(x=Group1, y=sqrt(Fungal)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of fungal polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
kruskal.test(Fungal ~ Group1, data = Fungal_Group)
#Kruskal-Wallis chi-squared = 28.306, df = 3, p-value = 3.133e-06

#Digestion of Algal Polysaccharides (GH96)
Algal_poly <- read.csv(file="Algal_poly", sep = "\t", row.names = 1, header =T)
Algal <- colSums(Algal_poly)
Algal <- data.frame(Algal)
Algal_Group <- cbind(Algal, metadata)
jpeg("Algal.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = Algal_Group, aes(x=Group1, y=sqrt(Algal)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of algal polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
kruskal.test(Algal ~ Group1, data = Algal_Group)
#Kruskal-Wallis chi-squared = 43.125, df = 3, p-value = 2.315e-09

#-- Lignin degradation (AA1, AA3_2, AA7)
Lignin_degradation <- read.csv(file="Lignin_degradation", sep = "\t", row.names = 1, header =T)
Lignin <- colSums(Lignin_degradation)
Lignin <- data.frame(Lignin)
Lignin_Group <- cbind(Lignin, metadata)
jpeg("Lignin.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = Lignin_Group, aes(x=Group1, y=sqrt(Lignin)*100, fill=Group1)) + geom_boxplot() + ggtitle("Lignin degradation") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
kruskal.test(Lignin ~ Group1, data = Lignin_Group)
#Kruskal-Wallis chi-squared = 10.862, df = 3, p-value = 0.0125

#--#Biosynthesis of Polysaccharides (GT24, GT88, GT39, GT64, GT10, GT21)
Biosynthesis_poly <- read.csv(file="Biosynthesis_poly", sep = "\t", row.names = 1, header =T)
Biosynthesis <- colSums(Biosynthesis_poly)
Biosynthesis <- data.frame(Biosynthesis)
Biosynthesis_Group <- cbind(Biosynthesis, metadata)
jpeg("Biosynthesis.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = Biosynthesis_Group, aes(x=Group1, y=sqrt(Biosynthesis)*100, fill=Group1)) + geom_boxplot() + ggtitle("Biosynthesis of polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) 
dev.off ()
kruskal.test(Biosynthesis ~ Group1, data = Biosynthesis_Group) 
#Kruskal-Wallis chi-squared = 20.971, df = 3, p-value = 0.0001067


#---- New HeatMap
Selected_cazy <- rbind(Animal_poly, Plant_poly, Algal_poly, Fungal_poly, Lignin_degradation, Biosynthesis_poly)
Selected_cazy <- data.frame(t(Selected_cazy))
Selected_cazy_Groups <- cbind (Selected_cazy, categories)
write.table(Selected_cazy_Groups, file="Selected_cazy_Groups.txt", sep = "\t")
#--- From this file GH27 removed and samples formated by group in excel
Selected_cazy_new <- read.csv(file="Selected_cazy_Groups_format.txt", sep = "\t", row.names = 1, header =T)
categories <- Selected_cazy_new[,20:21]
Selected_cazy_pred <- Selected_cazy_new[,1:19]
rownames(categories) <- rownames(Selected_cazy_pred)
annot.color.col <- list('Group'=c('darkgreen', 'darkred'),'Group1'=c("brown4", "orangered3", "green4", "limegreen"))
aheatmap(sqrt(Selected_cazy_pred)*100, color = "-RdYlBu2:100", breaks = 0,  main = "Association", distfun = "spearman", hclustfun = "complete", fontsize=10,  filename="CaZy_Finally_Used.png", scale = "row", Rowv = NA, Colv = TRUE, border_color = NA, legend = TRUE, cexRow = .7, cexCol = .7, annRow = categories, annCol = NA, annColors = annot.color.col, annLegend = TRUE)

# annRow = categories, annCol = NA, annColors = annot.color.col, annLegend = TRUE)

