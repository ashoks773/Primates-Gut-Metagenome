# Author: Ashok Kumar Sharma

setwd("~/Work/Metagenomic_Analysis/Final_Analysis/XenobioticsDegradation")
Xeno <- read.csv(file = "Final_Data/XenoEnzyme_ABUNDANCE.txt", row.names=1, header = T, sep = "\t")
Xeno_proportions <- Xeno/colSums(Xeno)[col(Xeno)]
Xeno1 <- data.frame(t(Xeno_proportions))
#write.table (Xeno_proportions, file = "Xeno_relative_abundances", sep = "\t")

metadata <- read.csv (file = "Final_Data/Metadata.txt", row.names=1, header = T, sep = "\t")

# Step 1: Principle coordinate analyis using ARDB Gene
Xeno_bray_pcoa <-pcoa(vegdist(Xeno1, "bray"))
Groups <- metadata[,1:4]
Xeno_bray_pcoa$values[1:2,]
mds.var.per1 <- round(Xeno_bray_pcoa$values$Rel_corr_eig[1]*100, 2)
mds.var.per2 <- round(Xeno_bray_pcoa$values$Rel_corr_eig[2]*100, 2)
#mds.var.per = round(Xeno_bray_pcoa$values$Eigenvalues/sum(Xeno_bray_pcoa$values$Eigenvalues)*100, 1)
Xeno_Bray_PCoA_MATRIX <- Xeno_bray_pcoa$vectors[,1:2]
Xeno_Bray_PCoA_MATRIX <- data.frame(Xeno_Bray_PCoA_MATRIX)
Xeno_Bray_PCoA_MATRIX_New <- cbind(Xeno_Bray_PCoA_MATRIX, Groups)

Color <- c("skyblue4", "wheat4")
Colors <- c("skyblue4", "skyblue", "wheat4", "wheat")

jpeg("Xeno_Bray_PCoA_Group1.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(Xeno_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCoA1 - ", mds.var.per1, "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per2, "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

distance_bray <-vegdist(Xeno1, "bray")
adonis(distance_bray ~ Group*Group1, data=metadata, permutations=99)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#Group      1   0.24726 0.247265 20.1542 0.24109   0.01 **
 # Group1     2   0.20174 0.100868  8.2216 0.19670   0.01 **
  
#----- Plot Axis2
pcoa2 <- Xeno_Bray_PCoA_MATRIX_New[,c(2,4)]
pcoa2_Melted <- melt(pcoa2, id.vars = "Group1")
jpeg("Bray_PCoA2_Distances.jpg", height = 4, width = 2, units = 'in', res = 600)
ggplot(data = pcoa2_Melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCo2") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
kruskalmc(pcoa2$Axis.2 ~ pcoa2$Group1)
#               obs.dif critical.dif difference
#BaAka-Bantu 20.5714286     14.82396       TRUE
#BaAka-Dry    0.3766234     15.80240      FALSE
#BaAka-Wet    0.6309524     15.42927      FALSE
#Bantu-Dry   20.9480519     15.80240       TRUE
#Bantu-Wet   19.9404762     15.42927       TRUE
#Dry-Wet      1.0075758     16.37157      FALSE

#--- To check difference of BaAka Bantu with Gorilla - Irrespective of seasons
#pcoa2_newgroup <- Xeno_Bray_PCoA_MATRIX_New[,c(2,6)]
#pcoa2_newgroup_Melted <- melt(pcoa2_newgroup, id.vars = "GroupNew")
#Color <- c("skyblue4", "skyblue", "wheat4")
#jpeg("Bray_PCoA2_Distances_NewGroup.jpg", height = 4, width = 2, units = 'in', res = 600)
#ggplot(data = pcoa2_newgroup_Melted, aes(x=GroupNew, y=value, fill=GroupNew)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCo2") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#dev.off ()
#kruskalmc(pcoa2_newgroup$Axis.2 ~ pcoa2_newgroup$GroupNew)
#obs.dif critical.dif difference
#BaAka-Bantu   20.5714286     13.45140       TRUE
#BaAka-Gorilla  0.1490683     12.06395      FALSE
#Bantu-Gorilla 20.4223602     12.06395       TRUE


#------- PCoA for Gorilla
#---- Gorilla
Color <- c("wheat4", "wheat")
metadata_gorilla <- read.csv (file = "Final_Data/Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
xeno1_gorilla <- Xeno1[1:23,]

Xeno_bray_pcoa_gorilla <-pcoa(vegdist(xeno1_gorilla, "bray"))
Xeno_bray_pcoa_gorilla$values[1:2,]
mds.var.per1 <- round(Xeno_bray_pcoa_gorilla$values$Rel_corr_eig[1]*100, 2)
mds.var.per2 <- round(Xeno_bray_pcoa_gorilla$values$Rel_corr_eig[2]*100, 2)
#mds.var.per_g = round(Xeno_bray_pcoa_gorilla$values$Eigenvalues/sum(Xeno_bray_pcoa_gorilla$values$Eigenvalues)*100, 1)
Xeno_Bray_PCoA_MATRIX_G <- Xeno_bray_pcoa_gorilla$vectors[,1:2]
Xeno_Bray_PCoA_MATRIX_G <- data.frame(Xeno_Bray_PCoA_MATRIX_G)
Xeno_Bray_PCoA_MATRIX_G_New <- cbind(Xeno_Bray_PCoA_MATRIX_G, Groups)

jpeg("Xeno_PCoA_Gorilla_Group1.jpg", height = 2.5, width = 3.5, units = 'in', res = 600)
ggplot(Xeno_Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per1, "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per2, "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

gorilla_dist_bray <- vegdist(xeno1_gorilla, "bray")
adonis(gorilla_dist_bray ~ Group*Group1, data=metadata_gorilla, permutations=99)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#Group1     1  0.033454 0.033454  4.1457 0.16487   0.01 **
  
#----- For total xenobiotics degradation ability plot
#---- Gorilla
Color <- c("wheat4", "wheat")
metadata_gorilla <- read.csv (file = "Final_Data/Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
Xeno2 <- data.frame(t(Xeno))
xeno2_gorilla <- Xeno2[1:23,]
xeno2_gorilla_sum <- rowSums(xeno2_gorilla)
xeno2_gorilla_sum_Group <- cbind (xeno2_gorilla_sum, metadata_gorilla)

jpeg("Xeno_cumulative_abundance_new_Gorilla.jpg", height = 3, width = 2, units = 'in', res = 600)
ggplot(aes(y = sqrt(xeno2_gorilla_sum)*100, x = Group1, fill = Group1), data = xeno2_gorilla_sum_Group) + geom_boxplot() + ggtitle("Xenobiotics Degradation") + labs(x="",y="sqrt(Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
wilcox.test(xeno2_gorilla_sum_Group$xeno2_gorilla_sum~xeno2_gorilla_sum_Group$Group1, alternative = "less", p.adjust.method="fdr")
#W = 34, p-value = 0.02561

#---- Human
Color <- c("skyblue4", "skyblue")
metadata_human <- read.csv (file = "Final_Data/Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_human[,2:3]
Xeno2 <- data.frame(t(Xeno))
xeno2_human <- Xeno2[24:51,]
xeno2_human_sum <- rowSums(xeno2_human)
xeno2_human_sum_Group <- cbind (xeno2_human_sum, metadata_human)

jpeg("Xeno_cumulative_abundance_new_Human.jpg", height = 3, width = 2, units = 'in', res = 600)
ggplot(aes(y = sqrt(xeno2_human_sum)*100, x = Group1, fill = Group1), data = xeno2_human_sum_Group) + geom_boxplot() + ggtitle("Xenobiotics Degradation") + labs(x="",y="sqrt(Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
wilcox.test(xeno2_human_sum_Group$xeno2_human_sum ~ xeno2_human_sum_Group$Group1, p.adjust.method="fdr", alternative = "less")
#W = 61, p-value = 0.04693

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


#BaAka_Bantu <- read.csv(file = "BaAka_vs_Bantu_SelNames.txt", sep="\t", row.names = 1, header = T)
#c1 <- row.names(BaAka_Bantu)
#c2 <- BaAka_Bantu[,1]
#c3 <- BaAka_Bantu[,7]
#c1 <- data.frame(c1)
#c3 <- data.frame(c3)
#c2 <- data.frame(c2)
#diff_df <- cbind(c1, c2, c3)
#colnames(diff_df) <- c("Xeno", "Qval", "log2FoldChange")

#library(calibrate)
#png(filename="Odds_ratio_Names.png", height = 3.5, width = 3.5, units = 'in', res = 600)
#with(diff_df, plot(log2FoldChange, -log10(Qval), pch=20, main="Volcano plot", xlim=c(-7,2), ylim=c(0,3.5)))
# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
#with(subset(diff_df, Qval<.05 ), points(log2FoldChange, -log10(Qval), pch=20, col="red"))
#with(subset(diff_df, abs(log2FoldChange)>1), points(log2FoldChange, -log10(Qval), pch=20, col="orange"))
#with(subset(diff_df, Qval<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(Qval), pch=20, col="green"))
# Label points with the textxy function from the calibrate plot
#with(subset(diff_df, Qval<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(Qval), labs=Xeno, cex=.3))
#dev.off ()

#----- Odds ratio New plot
#library(dplyr)
#diff_df_fil <- diff_df %>%
#  mutate(threshold = factor(case_when(log2FoldChange > 1 & pvalue < 0.05 ~ "BaAka",
#                                      log2FoldChange < -1 & pvalue < 0.05 ~ "Bantu",
#                                      TRUE ~ "NS")))
#png(filename="BaAka_Bantu_Odds_Ration_New.png", height = 3, width = 4, units = 'in', res = 600)
#ggplot(data=diff_df_fil, aes(x=log2FoldChange, y=-log10(pvalue))) + 
#  geom_point(aes(color=diff_df_fil$threshold), alpha=0.4, size=abs(diff_df_fil$log2FoldChange))+ 
#  geom_vline(xintercept=c(-1.5,1.5), color="red", alpha=1.0)+ 
#  geom_hline(yintercept=1.301, color="blue", alpha=1.0)+ 
#  xlab("log2FoldChange")+ 
#  ylab("-log10 (pvalue)")+ 
#  theme_bw()+
#  xlim(c(-7, 7)) +
#  ylim(c(0, 4)) +
#  scale_color_manual(name = "Threshold",
#                     values = c("BaAka" = "brown4", "Bantu" = "orangered3", "NS" = "darkgrey"))
#dev.off ()

#-----Selected xenobiotic degradation genes BoxPLots - Plotted Again
Sel_xeno_rel_abun <- read.csv(file = "Sel_xeno_Rel_abun.txt", row.names=1, header = T, sep = "\t")
Sel_xeno_rel_abun <- data.frame (t(Sel_xeno_rel_abun))
metadata <- read.csv (file = "Final_Data/Metadata.txt", row.names=1, header = T, sep = "\t")
Sel_xeno_rel_abun_Metadata <- cbind(Sel_xeno_rel_abun, metadata)

#-- P32320
p32320 <- Sel_xeno_rel_abun_Metadata[,c(1,10)]
p32320_Melted <- melt(p32320, id.vars = "Group1")
jpeg("p32320_new.jpg", height = 3.5, width = 1.8, units = 'in', res = 600)
ggplot(data = p32320_Melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("") + labs(x="",y="sqrt (relative abundance)") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
kruskalmc(p32320$P32320 ~ p32320$Group1)

##-- P49915
p49915 <- Sel_xeno_rel_abun_Metadata[,c(5,10)]
p49915_Melted <- melt(p49915, id.vars = "Group1")
jpeg("p49915_new.jpg", height = 3.5, width = 1.8, units = 'in', res = 600)
ggplot(data = p49915_Melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("") + labs(x="",y="sqrt (relative abundance)") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
kruskalmc(p49915$P49915 ~ p49915$Group1)

##-- Q99707
q99707 <- Sel_xeno_rel_abun_Metadata[,c(3,10)]
q99707_Melted <- melt(q99707, id.vars = "Group1")
jpeg("q99707_new.jpg", height = 3.5, width = 1.8, units = 'in', res = 600)
ggplot(data = q99707_Melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("") + labs(x="",y="sqrt (relative abundance)") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
kruskalmc(q99707$Q99707 ~ q99707$Group1)

##-- Q99798
q99798 <- Sel_xeno_rel_abun_Metadata[,c(8,10)]
q99798_Melted <- melt(q99798, id.vars = "Group1")
jpeg("q99798_new.jpg", height = 3.5, width = 1.8, units = 'in', res = 600)
ggplot(data = q99798_Melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("") + labs(x="",y="sqrt (relative abundance)") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
kruskalmc(q99798$Q99798 ~ q99798$Group1)

#--- To check significane using wilcox test
Sel_xeno_rel_abun_Metadata_Gorilla <- subset(Sel_xeno_rel_abun_Metadata, Group =="Gorilla")
wilcox.test(Sel_xeno_rel_abun_Metadata_Gorilla$P32320 ~Sel_xeno_rel_abun_Metadata_Gorilla$Group1, p.adjust.method="fdr")
Sel_xeno_rel_abun_Metadata_Human <- subset(Sel_xeno_rel_abun_Metadata, Group =="Human")
wilcox.test(Sel_xeno_rel_abun_Metadata_Human$P32320 ~Sel_xeno_rel_abun_Metadata_Human$Group1,p.adjust.method="fdr")


#--------- Supplementary Plots
