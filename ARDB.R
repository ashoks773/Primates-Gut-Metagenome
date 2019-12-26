# Step 1: Principle coordinate analyis using ARDB Gene
ARG <- read.csv(file = "ARDB_GENE_ABUNDANCE", row.names=1, header = T, sep = "\t")
ARG_proportions <- ARG/colSums(ARG)[col(ARG)]
ARG1 <- data.frame(t(ARG_proportions))
write.table (ARG1, file = "ARG_relative_abundances", sep = "\t")
write.table (ARG_proportions, file = "ARG_relative_abundances_t", sep = "\t")
metadata <- read.csv (file = "../Metadata.txt", row.names=1, header = T, sep = "\t")

Color <- c("darkgreen", "darkred")
Colors <- c("brown4", "orangered3", "green4", "limegreen")

# Abundance plot
ARG2 <- data.frame(t(ARG))
ARG1_sum <- rowSums(ARG2)
ARG1_sum_Group <- cbind (ARG1_sum, metadata)
jpeg("ARG_cumulative_abundance_new.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = sqrt(ARG1_sum)*100, x = Group, fill = Group1), data = ARG1_sum_Group) + geom_boxplot() + ggtitle("Antibiotic resistance") + labs(x="",y="sqrt(Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
write.table (ARG1_sum_Group, file = "ARG_cumulative_abundances", sep = "\t")

ARG1_sum <- data.frame(ARG1_sum)
ARG_gorilla <- ARG1_sum[1:23,]
ARG_human <- ARG1_sum[24:51,]
#mean and standard deviation

#----- human and gorilla individual plots

#---- Gorilla
Color <- c("green4", "limegreen")
metadata_gorilla <- read.csv (file = "../Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
ARG2 <- data.frame(t(ARG))
ARG2_gorilla <- ARG2[1:23,]
ARG2_gorilla_sum <- rowSums(ARG2_gorilla)
ARG2_gorilla_sum_Group <- cbind (ARG2_gorilla_sum, metadata_gorilla)

jpeg("ARG_cumulative_abundance_new_Gorilla.jpg", height = 3, width = 2, units = 'in', res = 600)
ggplot(aes(y = sqrt(ARG2_gorilla_sum)*100, x = Group1, fill = Group1), data = ARG2_gorilla_sum_Group) + geom_boxplot() + ggtitle("Antibiotic Resistance") + labs(x="",y="sqrt(Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
wilcox.test(ARG2_gorilla_sum_Group$ARG2_gorilla_sum~ARG2_gorilla_sum_Group$Group1)
#W = 38, p-value = 0.09084
#---- Human
Color <- c("brown4", "orangered3")
metadata_human <- read.csv (file = "../Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_human[,2:3]
ARG2 <- data.frame(t(ARG))
ARG2_human <- ARG2[24:51,]
ARG2_human_sum <- rowSums(ARG2_human)
ARG2_human_sum_Group <- cbind (ARG2_human_sum, metadata_human)

jpeg("ARG_cumulative_abundance_new_Human.jpg", height = 3, width = 2, units = 'in', res = 600)
ggplot(aes(y = sqrt(ARG2_human_sum)*100, x = Group1, fill = Group1), data = ARG2_human_sum_Group) + geom_boxplot() + ggtitle("Antibiotic Resistance") + labs(x="",y="sqrt(Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
wilcox.test(ARG2_human_sum_Group$ARG2_human_sum~ARG2_human_sum_Group$Group1)
#W = 60, p-value = 0.08493

#------------------------------ total
ARG_bray_pcoa <-pcoa(vegdist(ARG1, "bray"))
Groups <- metadata[,1:3]
ARG_bray_pcoa$values[1:2,]
mds.var.per = round(ARG_bray_pcoa$values$Eigenvalues/sum(ARG_bray_pcoa$values$Eigenvalues)*100, 1)
ARG_Bray_PCoA_MATRIX <- ARG_bray_pcoa$vectors[,1:2]
ARG_Bray_PCoA_MATRIX <- data.frame(ARG_Bray_PCoA_MATRIX)
ARG_Bray_PCoA_MATRIX_New <- cbind(ARG_Bray_PCoA_MATRIX, Groups)


jpeg("ARG_Bray_PCoA_Group.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(ARG_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
jpeg("ARG_Bray_PCoA_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(ARG_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCoA1 - ", mds.var.per[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

distance_bray <-vegdist(ARG1, "bray")
adonis(distance_bray ~metadata$Group1)
adonis(distance_bray ~ Group1, data = metadata)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Group1     3    7.8239 2.60798  37.287 0.70414  0.001 ***

#---- Gorilla
Color <- c("green4", "limegreen")
metadata_gorilla <- read.csv (file = "../Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
ARG_gorilla <- ARG1[1:23,]
ARG_bray_pcoa_gorilla <-pcoa(vegdist(ARG_gorilla, "bray"))
ARG_bray_pcoa_gorilla$values[1:2,]
mds.var.per_g = round(ARG_bray_pcoa_gorilla$values$Eigenvalues/sum(ARG_bray_pcoa_gorilla$values$Eigenvalues)*100, 1)
ARG_Bray_PCoA_MATRIX_G <- ARG_bray_pcoa_gorilla$vectors[,1:2]
ARG_Bray_PCoA_MATRIX_G <- data.frame(ARG_Bray_PCoA_MATRIX_G)
ARG_Bray_PCoA_MATRIX_G_New <- cbind(ARG_Bray_PCoA_MATRIX_G, Groups)
jpeg("ARG_PCoA_Gorilla_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(ARG_Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
jpeg("ARG_PCoA_Gorilla_Group2.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(ARG_Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1, shape=Group2)) + geom_point(size=3) + geom_point(aes(shape=Group2)) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted") + geom_point(shape=18)
dev.off ()

dist_bray_gorilla <-vegdist(ARG_gorilla, "bray")
adonis(dist_bray_gorilla ~ Group2*Group1, data=metadata_gorilla, permutations=99)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#Group2         1   0.08859 0.088591  2.6800 0.10843   0.01 **
 # Group1         1   0.05592 0.055917  1.6915 0.06844   0.19  

#---- Human
metadata_human <- read.csv (file = "../Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_human[,2:3]
ARG_human <- ARG1[24:51,]
ARG_bray_pcoa_human <-pcoa(vegdist(ARG_human, "bray"))
ARG_bray_pcoa_human$values[1:2,]
mds.var.per_g = round(ARG_bray_pcoa_human$values$Eigenvalues/sum(ARG_bray_pcoa_human$values$Eigenvalues)*100, 1)
ARG_Bray_PCoA_MATRIX_H <- ARG_bray_pcoa_human$vectors[,1:2]
ARG_Bray_PCoA_MATRIX_H <- data.frame(ARG_Bray_PCoA_MATRIX_H)
ARG_Bray_PCoA_MATRIX_H_New <- cbind(ARG_Bray_PCoA_MATRIX_H, Groups)

Color <- c("brown4", "orangered3")
jpeg("Cazy_PCoA_Human_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(ARG_Bray_PCoA_MATRIX_H_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

dist_bray_human <-vegdist(ARG_human, "bray")
adonis(dist_bray_human ~ Group1, data=metadata_human, permutations=99)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Group1     1   0.60167 0.60167  6.1973 0.19248   0.01 **

#Step2--- Statistical Analysis
ARG_human <- ARG1[24:51,]
Group1 <- metadata_human$Group1
Group1 <- data.frame(Group1)
ARG1_group1 <- cbind (ARG_human, Group1)
source('~/Work/Metagenomic_Analysis/Final_Analysis/ARDB_Analysis/BaAka_Vs_Bantu/WilcoxPairwise.R')
Wilcoxon_pairwise(ARG1_group1)
source('~/Work/Metagenomic_Analysis/Final_Analysis/ARDB_Analysis/BaAka_Vs_Bantu/Odds_Ratio_Calculation.R')
Calculate_Odds_Ratio(ARG1_group1)

ARG_gorilla <- ARG1[1:23,]
Group1 <- metadata_gorilla$Group1
Group1 <- data.frame(Group1)
ARG1_group1 <- cbind (ARG_gorilla, Group1)
source('~/Work/Metagenomic_Analysis/Final_Analysis/ARDB_Analysis/Dry_Vs_Wet/WilcoxPairwise.R')
Wilcoxon_pairwise(ARG1_group1)
source('~/Work/Metagenomic_Analysis/Final_Analysis/ARDB_Analysis/Dry_Vs_Wet/Odds_Ratio_Calculation.R')
Calculate_Odds_Ratio(ARG1_group1)

#-- Selected Top12 from both comparisions
Selected12 <- read.csv(file="SelectedARG_relative_abundance.txt", sep = "\t", row.names = 1, header =TRUE)
ARG_selected <- data.frame(Selected12)
categories <- metadata[,1:2]
rownames(categories) <- colnames(ARG_selected)
annot.color.col <- list('Group'=c('darkgreen', 'darkred'),'Group1'=c("brown4", "orangered3", "green4", "limegreen"))
jpeg("DesHeatMap_Spearman_Using_Selected.jpg", height = 4, width = 8, units = 'in', res = 600)
aheatmap(sqrt(sqrt(ARG_selected)), color = "-RdYlBu2:100", breaks = 0,  main = "Association", distfun = "minkowski", hclustfun = "complete", fontsize=10, scale = "col", Rowv = TRUE, Colv = TRUE, border_color = NA, legend = TRUE, cexRow = .7, cexCol = .7, annCol = categories, annRow = NA, annColors = annot.color.col, annLegend = TRUE)
dev.off ()
#--- Two sample removed
Selected12 <- read.csv(file="SelectedARG_relative_abundance_twoRemoved.txt", sep = "\t", row.names = 1, header =TRUE)
ARG_selected <- data.frame(Selected12)
metadata <- read.csv (file = "Metadata_tworemoved.txt", row.names=1, header = T, sep = "\t")
categories <- metadata[,1:2]
rownames(categories) <- colnames(ARG_selected)
annot.color.col <- list('Group'=c('darkgreen', 'darkred'),'Group1'=c("brown4", "orangered3", "green4", "limegreen"))
jpeg("DesHeatMap_Manhattan_Using_Selected_twoRemoved.jpg", height = 4, width = 8, units = 'in', res = 600)
aheatmap(sqrt(sqrt(ARG_selected)), color = "-RdYlBu2:100", breaks = 0,  main = "Association", distfun = "manhattan", hclustfun = "complete", fontsize=10, scale = "col", Rowv = TRUE, Colv = TRUE, border_color = NA, legend = TRUE, cexRow = .7, cexCol = .7, annCol = categories, annRow = NA, annColors = annot.color.col, annLegend = TRUE)
dev.off ()

#--- New HeatMap odering seperately in Excel for Groups
Selected12 <- read.csv(file="SelectedARG_relative_abundance_twoRemoved.txt", sep = "\t", row.names = 1, header =TRUE)
ARG_selected <- data.frame(t(Selected12))
metadata <- read.csv (file = "Metadata_tworemoved.txt", row.names=1, header = T, sep = "\t")
categories <- metadata[,1:2]
ARG_selected_Groups <- cbind(ARG_selected, categories)
write.table (ARG_selected_Groups, file = "ARG_selected_Groups.txt", sep = "\t")

Selected12_formatted <- read.csv(file="ARG_selected_groups_formatted.txt", sep = "\t", row.names = 1, header =TRUE)
categories <- Selected12_formatted[,13:14]
Selected12_pred <- Selected12_formatted[,1:12]
rownames(categories) <- rownames(Selected12_pred)
annot.color.col <- list('Group'=c('darkgreen', 'darkred'),'Group1'=c("brown4", "orangered3", "green4", "limegreen"))
jpeg("Selected_HeatmapFinallyUsed.jpg", height = 9, width = 6, units = 'in', res = 600)
aheatmap(sqrt(sqrt(Selected12_pred)), color = "-RdYlBu2:100", breaks = 0,  main = "Association", distfun = "manhattan", hclustfun = "complete", fontsize=10, scale = "row", Rowv = NA, Colv = TRUE, border_color = NA, legend = TRUE, cexRow = .7, cexCol = .9, annRow = categories, annCol = NA, annColors = annot.color.col, annLegend = TRUE)
dev.off ()
