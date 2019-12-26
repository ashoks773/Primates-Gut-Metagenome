setwd("~/Work/Metagenomic_Analysis/Final_Analysis/XenobioticsDegradation")
Xeno <- read.csv(file = "XenoEnzyme_ABUNDANCE.txt", row.names=1, header = T, sep = "\t")
Xeno_proportions <- Xeno/colSums(Xeno)[col(Xeno)]
Xeno1 <- data.frame(t(Xeno_proportions))
write.table (Xeno_proportions, file = "Xeno_relative_abundances", sep = "\t")

metadata <- read.csv (file = "../Metadata.txt", row.names=1, header = T, sep = "\t")

# Step 1: Principle coordinate analyis using ARDB Gene
Xeno_bray_pcoa <-pcoa(vegdist(Xeno1, "bray"))
Groups <- metadata[,1:3]
Xeno_bray_pcoa$values[1:2,]
mds.var.per = round(Xeno_bray_pcoa$values$Eigenvalues/sum(Xeno_bray_pcoa$values$Eigenvalues)*100, 1)
Xeno_Bray_PCoA_MATRIX <- Xeno_bray_pcoa$vectors[,1:2]
Xeno_Bray_PCoA_MATRIX <- data.frame(Xeno_Bray_PCoA_MATRIX)
Xeno_Bray_PCoA_MATRIX_New <- cbind(Xeno_Bray_PCoA_MATRIX, Groups)

Color <- c("darkgreen", "darkred")
Colors <- c("brown4", "orangered3", "green4", "limegreen")

jpeg("Xeno_Bray_PCoA_Group.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Xeno_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
jpeg("Xeno_Bray_PCoA_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Xeno_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCoA1 - ", mds.var.per[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

distance_bray <-vegdist(Xeno1, "bray")
adonis(distance_bray ~ Group1, data=metadata, permutations=99)


#------- PCoA for Gorilla
#---- Gorilla
Color <- c("green4", "limegreen")
metadata_gorilla <- read.csv (file = "../Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
xeno1_gorilla <- Xeno1[1:23,]

Xeno_bray_pcoa_gorilla <-pcoa(vegdist(xeno1_gorilla, "bray"))
Xeno_bray_pcoa_gorilla$values[1:2,]
mds.var.per_g = round(Xeno_bray_pcoa_gorilla$values$Eigenvalues/sum(Xeno_bray_pcoa_gorilla$values$Eigenvalues)*100, 1)
Xeno_Bray_PCoA_MATRIX_G <- Xeno_bray_pcoa_gorilla$vectors[,1:2]
Xeno_Bray_PCoA_MATRIX_G <- data.frame(Xeno_Bray_PCoA_MATRIX_G)
Xeno_Bray_PCoA_MATRIX_G_New <- cbind(Xeno_Bray_PCoA_MATRIX_G, Groups)

jpeg("Xeno_PCoA_Gorilla_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Xeno_Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

gorilla_dist_bray <- vegdist(xeno1_gorilla, "bray")
adonis(gorilla_dist_bray ~ Group1, data=metadata_gorilla, permutations=99)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#Group1     1  0.033454 0.033454  4.1457 0.16487   0.01 **
  
#----- For total xenobiotics degradation ability plot
# Abundance plot
Xeno2 <- data.frame(t(Xeno))
Xeno2_sum <- rowSums(Xeno2)
Xeno2_sum_Group <- cbind (Xeno2_sum, metadata)
jpeg("Xeno_cumulative_abundance_new.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(aes(y = sqrt(Xeno2_sum)*100, x = Group, fill = Group1), data = Xeno2_sum_Group) + geom_boxplot() + ggtitle("Xenobiotics Degradation") + labs(x="",y="sqrt(Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
write.table (Xeno2_sum_Group, file = "Xeno_Total_cumulative_abundance", sep = "\t")
#-- Krukshal Wallis
Stats_xeno <- kruskalmc(Xeno2_sum_Group$Xeno2_sum~Xeno2_sum_Group$Group1, probs=0.05)
Stats_xeno_p0.05 <- Stats_xeno$dif.com
write.table(Stats_xeno_p0.05, file = "Stats_xeno_p0.05.txt", sep= "\t")

#---- Gorilla
Color <- c("green4", "limegreen")
metadata_gorilla <- read.csv (file = "../Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
Xeno2 <- data.frame(t(Xeno))
xeno2_gorilla <- Xeno2[1:23,]
xeno2_gorilla_sum <- rowSums(xeno2_gorilla)
xeno2_gorilla_sum_Group <- cbind (xeno2_gorilla_sum, metadata_gorilla)

jpeg("Xeno_cumulative_abundance_new_Gorilla.jpg", height = 3, width = 2, units = 'in', res = 600)
ggplot(aes(y = sqrt(xeno2_gorilla_sum)*100, x = Group1, fill = Group1), data = xeno2_gorilla_sum_Group) + geom_boxplot() + ggtitle("Xenobiotics Degradation") + labs(x="",y="sqrt(Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
wilcox.test(xeno2_gorilla_sum_Group$xeno2_gorilla_sum~xeno2_gorilla_sum_Group$Group1)
#W = 34, p-value = 0.05122

#---- Human
Color <- c("brown4", "orangered3")
metadata_human <- read.csv (file = "../Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_human[,2:3]
Xeno2 <- data.frame(t(Xeno))
xeno2_human <- Xeno2[24:51,]
xeno2_human_sum <- rowSums(xeno2_human)
xeno2_human_sum_Group <- cbind (xeno2_human_sum, metadata_human)

jpeg("Xeno_cumulative_abundance_new_Human.jpg", height = 3, width = 2, units = 'in', res = 600)
ggplot(aes(y = sqrt(xeno2_human_sum)*100, x = Group1, fill = Group1), data = xeno2_human_sum_Group) + geom_boxplot() + ggtitle("Xenobiotics Degradation") + labs(x="",y="sqrt(Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
wilcox.test(xeno2_human_sum_Group$xeno2_human_sum~xeno2_human_sum_Group$Group1, p=0.05)
#V = 18, p-value = 0.02954

#------- Significantly discrimintaing XDEx between Groilla group
Group1 <- metadata_gorilla$Group1
xeno1_gorilla <- Xeno1[1:23,]
Xeno_gorilla_Group1 <- cbind (xeno1_gorilla, Group1)

source('~/Work/Metagenomic_Analysis/Final_Analysis/XenobioticsDegradation/Dry_Wet/WilcoxPairwise.R')
Wilcoxon_pairwise(Xeno_gorilla_Group1)
source('~/Work/Metagenomic_Analysis/Final_Analysis/XenobioticsDegradation/Dry_Wet/Odds_Ratio_Calculation.R')
Calculate_Odds_Ratio(Xeno_gorilla_Group1)

Dry_wet <- read.csv(file = "Dry_wet_stats.txt", sep="\t", row.names = 1, header = T)
c1 <- row.names(Dry_wet)
c2 <- Dry_wet[,1]
c3 <- Dry_wet[,8]
c1 <- data.frame(c1)
c3 <- data.frame(c3)
c2 <- data.frame(c2)
diff_df <- cbind(c1, c2, c3)
colnames(diff_df) <- c("Xeno", "pvalue", "log2FoldChange")

library(calibrate)
png(filename="Odds_ratio.png", height = 4, width = 4, units = 'in', res = 600)
with(diff_df, plot(log2FoldChange, -log10(Qval), pch=20, main="Volcano plot", xlim=c(-3.5,4), ylim=c(0,5)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(diff_df, Qval<.05 ), points(log2FoldChange, -log10(Qval), pch=20, col="red"))
with(subset(diff_df, abs(log2FoldChange)>1), points(log2FoldChange, -log10(Qval), pch=20, col="orange"))
with(subset(diff_df, Qval<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(Qval), pch=20, col="green"))
# Label points with the textxy function from the calibrate plot
with(subset(diff_df, Qval<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(Qval), labs=Xeno, cex=.6))
dev.off ()

#--- Odds ratio New plot
library(dplyr)
diff_df_fil <- diff_df %>%
  mutate(threshold = factor(case_when(log2FoldChange > 1 & pvalue < 0.05 ~ "Dry",
                                      log2FoldChange < -1 & pvalue < 0.05 ~ "Wet",
                                      TRUE ~ "NS")))
png(filename="Dry_Wet_Odds_Ration_New.png", height = 3, width = 4, units = 'in', res = 600)
ggplot(data=diff_df_fil, aes(x=log2FoldChange, y=-log10(pvalue))) + 
  geom_point(aes(color=diff_df_fil$threshold), alpha=0.4, size=abs(diff_df_fil$log2FoldChange))+ 
  geom_vline(xintercept=c(-1.05,1.05), color="red", alpha=1)+ 
  geom_hline(yintercept=1.301, color="blue", alpha=1.0)+ 
  xlab("log2FoldChange")+ 
  ylab("-log10 (pvalue)")+ 
  theme_bw()+
  xlim(c(-6, 6)) +
  ylim(c(0, 5.5)) +
  scale_color_manual(name = "Threshold",
                     values = c("Dry" = "green4", "Wet" = "limegreen", "NS" = "darkgrey"))
dev.off ()


#--- Human
Xeno1_human <- Xeno1[24:51,]
Group1 <- metadata_human$Group1
Group1 <- data.frame(Group1)
Xeno1_human_group1 <- cbind (Xeno1_human, Group1)
source('~/Work/Metagenomic_Analysis/Final_Analysis/XenobioticsDegradation/BaAka_Bantu/WilcoxPairwise.R')
Wilcoxon_pairwise(Xeno1_human_group1)
source('~/Work/Metagenomic_Analysis/Final_Analysis/XenobioticsDegradation/BaAka_Bantu/Odds_Ratio_Calculation.R')
Calculate_Odds_Ratio(Xeno1_human_group1)


BaAka_Bantu <- read.csv(file = "BaAka_vs_Bantu_SelNames.txt", sep="\t", row.names = 1, header = T)
c1 <- row.names(BaAka_Bantu)
c2 <- BaAka_Bantu[,1]
c3 <- BaAka_Bantu[,7]
c1 <- data.frame(c1)
c3 <- data.frame(c3)
c2 <- data.frame(c2)
diff_df <- cbind(c1, c2, c3)
colnames(diff_df) <- c("Xeno", "Qval", "log2FoldChange")

library(calibrate)
png(filename="Odds_ratio_Names.png", height = 3.5, width = 3.5, units = 'in', res = 600)
with(diff_df, plot(log2FoldChange, -log10(Qval), pch=20, main="Volcano plot", xlim=c(-7,2), ylim=c(0,3.5)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(diff_df, Qval<.05 ), points(log2FoldChange, -log10(Qval), pch=20, col="red"))
with(subset(diff_df, abs(log2FoldChange)>1), points(log2FoldChange, -log10(Qval), pch=20, col="orange"))
with(subset(diff_df, Qval<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(Qval), pch=20, col="green"))
# Label points with the textxy function from the calibrate plot
with(subset(diff_df, Qval<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(Qval), labs=Xeno, cex=.3))
dev.off ()

#----- Odds ratio New plot
library(dplyr)
diff_df_fil <- diff_df %>%
  mutate(threshold = factor(case_when(log2FoldChange > 1 & pvalue < 0.05 ~ "BaAka",
                                      log2FoldChange < -1 & pvalue < 0.05 ~ "Bantu",
                                      TRUE ~ "NS")))
png(filename="BaAka_Bantu_Odds_Ration_New.png", height = 3, width = 4, units = 'in', res = 600)
ggplot(data=diff_df_fil, aes(x=log2FoldChange, y=-log10(pvalue))) + 
  geom_point(aes(color=diff_df_fil$threshold), alpha=0.4, size=abs(diff_df_fil$log2FoldChange))+ 
  geom_vline(xintercept=c(-1.5,1.5), color="red", alpha=1.0)+ 
  geom_hline(yintercept=1.301, color="blue", alpha=1.0)+ 
  xlab("log2FoldChange")+ 
  ylab("-log10 (pvalue)")+ 
  theme_bw()+
  xlim(c(-7, 7)) +
  ylim(c(0, 4)) +
  scale_color_manual(name = "Threshold",
                     values = c("BaAka" = "brown4", "Bantu" = "orangered3", "NS" = "darkgrey"))
dev.off ()


#-------- -Selected xenobiotic degradation genes
Sel_xeno_rel_abun <- read.csv(file = "Sel_xeno_Rel_abun.txt", row.names=1, header = T, sep = "\t")
Sel_xeno_rel_abun <- data.frame (t(Sel_xeno_rel_abun))
metadata <- read.csv (file = "../Metadata.txt", row.names=1, header = T, sep = "\t")

jpeg("P32320.jpg", height = 4, width = 2, units = 'in', res = 600)
boxplot(sqrt(Sel_xeno_rel_abun$P32320) ~ metadata$Group1, col=c("brown4", "orangered3", "green4", "limegreen"), cex.axis=0.6)
dev.off ()
kruskalmc(Sel_xeno_rel_abun$P32320 ~ metadata$Group1)
#BaAka-Bantu  7.857143     14.82396      FALSE
#Dry-Wet     13.431818     16.37157      FALSE

jpeg("O94808.jpg", height = 4, width = 2, units = 'in', res = 600)
boxplot(sqrt(Sel_xeno_rel_abun$O94808) ~ metadata$Group1, col=c("brown4", "orangered3", "green4", "limegreen"), cex.axis=0.6)
dev.off ()

jpeg("Q99707.jpg", height = 4, width = 2, units = 'in', res = 600)
boxplot(sqrt(Sel_xeno_rel_abun$Q99707) ~ metadata$Group1, col=c("brown4", "orangered3", "green4", "limegreen"), cex.axis=0.6)
dev.off ()
kruskalmc(Sel_xeno_rel_abun$Q99707 ~ metadata$Group1)
#BaAka-Bantu 24.428571     14.82396       TRUE
#Dry-Wet     18.712121     16.37157       TRUE

jpeg("P03697.jpg", height = 4, width = 2, units = 'in', res = 600)
boxplot(sqrt(Sel_xeno_rel_abun$P03697) ~ metadata$Group1, col=c("brown4", "orangered3", "green4", "limegreen"), cex.axis=0.6)
dev.off ()

jpeg("P49915.jpg", height = 4, width = 2, units = 'in', res = 600)
boxplot(sqrt(Sel_xeno_rel_abun$P49915) ~ metadata$Group1, col=c("brown4", "orangered3", "green4", "limegreen"), cex.axis=0.6)
dev.off ()
kruskalmc(Sel_xeno_rel_abun$P49915 ~ metadata$Group1)
#BaAka-Bantu  6.714286     14.82396      FALSE
#Dry-Wet     10.287879     16.37157      FALSE

jpeg("Q06830.jpg", height = 4, width = 2, units = 'in', res = 600)
boxplot(sqrt(Sel_xeno_rel_abun$Q06830) ~ metadata$Group1, col=c("brown4", "orangered3", "green4", "limegreen"), cex.axis=0.6)
dev.off ()

jpeg("P07203.jpg", height = 4, width = 2, units = 'in', res = 600)
boxplot(sqrt(Sel_xeno_rel_abun$P07203) ~ metadata$Group1, col=c("brown4", "orangered3", "green4", "limegreen"), cex.axis=0.6)
dev.off ()

jpeg("Q99798.jpg", height = 4, width = 2, units = 'in', res = 600)
boxplot(sqrt(Sel_xeno_rel_abun$Q99798) ~ metadata$Group1, col=c("brown4", "orangered3", "green4", "limegreen"), cex.axis=0.6)
dev.off ()
kruskalmc(Sel_xeno_rel_abun$Q99798 ~ metadata$Group1)
#BaAka-Bantu 20.0714286     14.82396       TRUE
#Dry-Wet     14.3939394     16.37157      FALSE



