# Author: Ashok Kumar Sharma


##########################
#------ Figure 5a and 5b
#########################
HumanTaxa <- read.csv (file = "SelectedTaxa_human.txt", sep = "\t", row.names = 1, header =T)
Group1 <- HumanTaxa$Group1
Group1 <- data.frame(Group1)
Color <- c("skyblue4", "skyblue")

Prevotella_combined <- HumanTaxa[,c(16,34,46:54)]
Prevotella_combined <- rowSums(Prevotella_combined)
Prevotella_combined <- data.frame(Prevotella_combined)
Prevotella_combined_sqrtper <- sqrt(Prevotella_combined)*100
Prevotella_combined_Group <- cbind (Prevotella_combined_sqrtper, Group1)

jpeg("H_Prevotella_new.jpg", height = 3, width = 1.7, units = 'in', res = 600)
boxplot(Prevotella_combined_Group$Prevotella_combined ~ Prevotella_combined_Group$Group1, col=Color, ylab="sqrt(Relative abundance) %", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8),outline=FALSE)
stripchart(Prevotella_combined ~ Group1, vertical = TRUE, data = Prevotella_combined_Group, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off()
wilcox.test (Prevotella_combined_Group$Prevotella_combined ~ Prevotella_combined_Group$Group1)
#W = 195, p-value = 9.971e-08

Treponema_combined <- HumanTaxa[,c(2,20,55)]
Treponema_combined <- rowSums(Treponema_combined)
Treponema_combined <- data.frame(Treponema_combined)
Treponema_combined_sqrtper <- sqrt(Treponema_combined)*100
Treponema_combined_Group <- cbind (Treponema_combined_sqrtper, Group1)

jpeg("H_Spirochaetaceae.jpg", height = 3, width = 1.7, units = 'in', res = 600)
boxplot(Treponema_combined_Group$Treponema_combined ~ Treponema_combined_Group$Group1, col=Color, ylab="sqrt(Relative abundance) %", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8),outline=FALSE)
stripchart(Treponema_combined ~ Group1, vertical = TRUE, data = Treponema_combined_Group, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off()
wilcox.test (Treponema_combined_Group$Treponema_combined ~ Treponema_combined_Group$Group1)


GorillaTaxa <- read.csv (file = "SelectedTaxa_gorilla.txt", sep = "\t", row.names = 1, header =T)
Group1 <- GorillaTaxa$Group1
Group1 <- data.frame(Group1)
Color <- c("wheat4", "wheat")

Prevotella_combined <- GorillaTaxa[,c(16,34,46:54)]
Prevotella_combined <- rowSums(Prevotella_combined)
Prevotella_combined <- data.frame(Prevotella_combined)
Prevotella_combined_sqrtper <- sqrt(Prevotella_combined)*100
Prevotella_combined_Group <- cbind (Prevotella_combined_sqrtper, Group1)

jpeg("G_Prevotella.jpg", height = 3, width = 1.7, units = 'in', res = 600)
boxplot(Prevotella_combined_Group$Prevotella_combined ~ Prevotella_combined_Group$Group1, col=Color, ylab="sqrt(Relative abundance) %", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8),outline=FALSE)
stripchart(Prevotella_combined ~ Group1, vertical = TRUE, data = Prevotella_combined_Group, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off()
wilcox.test (Prevotella_combined_Group$Prevotella_combined ~ Prevotella_combined_Group$Group1)
#W = 108, p-value = 0.008625

Treponema_combined <- GorillaTaxa[,c(2,20,55)]
Treponema_combined <- rowSums(Treponema_combined)
Treponema_combined <- data.frame(Treponema_combined)
Treponema_combined_sqrtper <- sqrt(Treponema_combined)*100
Treponema_combined_Group <- cbind (Treponema_combined_sqrtper, Group1)

jpeg("G_Spirochaetaceae_new.jpg", height = 3, width = 1.7, units = 'in', res = 600)
boxplot(Treponema_combined_Group$Treponema_combined ~ Treponema_combined_Group$Group1, col=Color, ylab="sqrt(Relative abundance) %", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8),outline=FALSE)
stripchart(Treponema_combined ~ Group1, vertical = TRUE, data = Treponema_combined_Group, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off()
wilcox.test (Treponema_combined_Group$Treponema_combined ~ Treponema_combined_Group$Group1)
#W = 59, p-value = 0.6947


################################################################
#------- Prevo and Trepo Correlations with Plant PolySaccharides
################################################################
Plant_poly <- read.csv(file="Final_Data/plant_poly", sep = "\t", row.names = 1, header =T)
Plant <- colSums(Plant_poly)
Plant <- data.frame(Plant)
Fungal_poly <- read.csv(file="Final_Data/Fungal_poly", sep = "\t", row.names = 1, header =T)
Fungal <- colSums(Fungal_poly)
Fungal <- data.frame(Fungal)
Algal_poly <- read.csv(file="Final_Data/Algal_poly", sep = "\t", row.names = 1, header =T)
Algal <- colSums(Algal_poly)
Algal <- data.frame(Algal)
Lignin_degradation <- read.csv(file="Final_Data/Lignin_degradation", sep = "\t", row.names = 1, header =T)
Lignin <- colSums(Lignin_degradation)
Lignin <- data.frame(Lignin)
Biosynthesis_poly <- read.csv(file="Final_Data/Biosynthesis_poly", sep = "\t", row.names = 1, header =T)
Biosynthesis <- colSums(Biosynthesis_poly)
Biosynthesis <- data.frame(Biosynthesis)
Animal_poly <- read.csv(file="Final_Data/Animal_poly", sep = "\t", row.names = 1, header =T)
Animal <- colSums(Animal_poly)
Animal <- data.frame(Animal)

#--- Taxa
taxa_Sel <- read.csv (file = "Final_Data/Selected_common55Taxa_proportions.txt", sep = "\t", header = T, row.names = 1)
taxa_Sel <- data.frame(t(taxa_Sel))
Prevotella_combined <- taxa_Sel[,c(16,34,46:54)]
Prevotella_combined <- rowSums(Prevotella_combined)
Prevotella_combined <- sqrt(Prevotella_combined)*100

Spirochate_combined <- taxa_Sel[,c(2,20,55)]
Spirochate_combined <- rowSums(Spirochate_combined)
Spirochate_combined <- sqrt(Spirochate_combined)*100

#Actinobacteria <- sqrt(taxa_Sel[,32])*100

#- Metadata
metadata <- read.csv(file="Final_Data/Metadata.txt", sep = "\t", row.names=1, header = T)

#--- Combined Table
Cazy_Taxa_combined_group <- cbind (Plant, Fungal, Animal, Lignin, Biosynthesis, Algal, Prevotella_combined, Spirochate_combined, metadata)

#---- Correlations
library("ggpubr")
Colors <- c("skyblue4", "skyblue", "wheat4", "wheat")

jpeg("Plant_prevo_new.jpg", height = 2, width = 3, units = 'in', res = 600)
ggplot(data = Cazy_Taxa_combined_group, mapping = aes(x = Plant, y = Prevotella_combined)) + 
  geom_point(mapping = aes(color = Group1)) + scale_color_manual(values=Colors) + 
  scale_shape_manual(values =1:20) + theme_bw() +
  xlab("Digestion of plant polysaccharides") + ylab("sqrt(Prevetella)%") +
  geom_smooth()
dev.off ()
Correlation <- corr.test(Cazy_Taxa_combined_group$Plant, Cazy_Taxa_combined_group$Prevotella_combined, method="spearman", adjust="fdr", alpha =0.5, use = "pairwise")
#r=0.7, p = 7.753896e-09

jpeg("Plant_Spiro_new.jpg", height = 2, width = 3, units = 'in', res = 600)
ggplot(data = Cazy_Taxa_combined_group, mapping = aes(x = Plant, y = Spirochate_combined)) + 
  geom_point(mapping = aes(color = Group1)) + scale_color_manual(values=Colors) + 
  scale_shape_manual(values =1:20) + theme_bw() +
  xlab("Digestion of plant polysaccharides") + ylab("sqrt(Spirochaetaceae)%") +
  geom_smooth()
dev.off ()
Correlation <- corr.test(Cazy_Taxa_combined_group$Plant, Cazy_Taxa_combined_group$Spirochate_combined, method="spearman", adjust="fdr", alpha =0.5, use = "pairwise")
#r=-0.6210956, p = 1.156344e-06


########################################
#------ Analysis on Reconstructed Bins
########################################
setwd("~/Work/Metagenomic_Analysis/Final_Analysis/GenomeReconstuction/Total_Prevo_SpiroBins")
bin38 <- read.csv(file = "bin.38.txt", sep = "\t", row.names = 1, header = T)
bin47 <- read.csv(file = "bin.47.txt", sep = "\t", row.names = 1, header = T)
bin89 <- read.csv(file = "bin.89.txt", sep = "\t", row.names = 1, header = T)
bin126 <- read.csv(file = "bin.126.txt", sep = "\t", row.names = 1, header = T)
bin148 <- read.csv(file = "bin.148.txt", sep = "\t", row.names = 1, header = T)
bin460 <- read.csv(file = "bin.460.txt", sep = "\t", row.names = 1, header = T)
bin571 <- read.csv(file = "bin.571.txt", sep = "\t", row.names = 1, header = T)
bin487 <- read.csv(file = "bin.487.txt", sep = "\t", row.names = 1, header = T)
bin52 <- read.csv(file = "bin.52.txt", sep = "\t", row.names = 1, header = T)
bin174 <- read.csv(file = "bin.174.txt", sep = "\t", row.names = 1, header = T)
bin191 <- read.csv(file = "bin.191.txt", sep = "\t", row.names = 1, header = T)
bin217 <- read.csv(file = "bin.217.txt", sep = "\t", row.names = 1, header = T)
bin245 <- read.csv(file = "bin.245.txt", sep = "\t", row.names = 1, header = T)
bin372 <- read.csv(file = "bin.372.txt", sep = "\t", row.names = 1, header = T)
bin478 <- read.csv(file = "bin.478.txt", sep = "\t", row.names = 1, header = T)
bin633 <- read.csv(file = "bin.633.txt", sep = "\t", row.names = 1, header = T)
bin678 <- read.csv(file = "bin.678.txt", sep = "\t", row.names = 1, header = T)
bin727 <- read.csv(file = "bin.727.txt", sep = "\t", row.names = 1, header = T)

library(plyr)
mylist <- list(bin38,bin47,bin89,bin126,bin148,bin460,bin571,bin487,bin52,bin174,bin191,bin217,bin245,bin372,bin478,bin633,bin678,bin727)

for(i in 1:length(mylist)){
  colnames(mylist[[i]]) <- paste0( names(mylist)[i], colnames(mylist[[i]]) )
  mylist[[i]]$ROWNAMES  <- rownames(mylist[[i]])
}

Total_bins_cazy <- join_all( mylist, by="ROWNAMES", type="full" )
rownames(Total_bins_cazy) <- Total_bins_cazy$ROWNAMES; Total_bins_cazy$ROWNAMES <- NULL
Total_bins_cazy[is.na(Total_bins_cazy)] <- 0

#------
Total_bins_cazy_t <- data.frame(t(Total_bins_cazy))
cazy.sums <- colSums(Total_bins_cazy_t)
Total_bins_cazy_filtered <- Total_bins_cazy_t[ ,which(cazy.sums > 5)] #Remove CaZy if total sum is not more than 5
Total_bins_cazy_final <- dropspc (Total_bins_cazy_filtered, 3)

meta <- read.csv(file = "Bin_details.txt", sep = "\t", row.names = 1, header = T)

N.obs <- rowSums(Total_bins_cazy_final > 0)
Shannon <- diversity (Total_bins_cazy_final)
InvSimp <- diversity (Total_bins_cazy_final, "invsimpson")

total_div <- data.frame(N.obs, Shannon, InvSimp, meta)


###################
#------ Figure 5d
###################

jpeg("Observed.jpg", height = 4, width = 2, units = 'in', res = 600)
boxplot(total_div$N.obs ~ total_div$Taxa, col=c("#CD853F", "#BDB76B"), ylab="Observed CaZy", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(N.obs ~ Taxa, vertical = TRUE, data = total_div, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()
wilcox.test(total_div$N.obs ~ meta$Taxa)

jpeg("Shnnon.jpg", height = 4, width = 2, units = 'in', res = 600)
boxplot(total_div$Shannon ~ total_div$Taxa, col=c("#CD853F", "#BDB76B"), ylab="Shannon's H",  cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(Shannon ~ Taxa, vertical = TRUE, data = total_div, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()
wilcox.test(total_div$Shannon ~ meta$Taxa)

###################
#------ Figure 5c
###################

#----- Beta-diversity analysis
library("FactoMineR")
library("factoextra")

#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/

Total_bins_cazy_final_relab <- decostand(Total_bins_cazy_final, method = "total")*100
adonis(Total_bins_cazy_final_relab ~ meta$Taxa)
#meta$Taxa  1    1.2721 1.27208  12.518 0.43896  0.001 ***

res.pca <- PCA(Total_bins_cazy_final_relab, graph = FALSE)
var <- get_pca_var(res.pca)

#---- Finally Used
#var$coord: coordinates of variables to create a scatter plot
#var$cos2: represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.
#var$contrib: contains the contributions (in percentage) of the variables to the principal components. The contribution of a variable (var) to a given principal component is (in percentage) : (var.cos2 * 100) / (total cos2 of the component).

groups <- meta$Taxa
jpeg("PCA.jpg", height = 4, width = 6, units = 'in', res = 600)
fviz_pca_ind(res.pca,col.ind=groups,palette=c("#CD853F", "#BDB76B"),addEllipses=TRUE,ellipse.type="confidence",legend.title="Groups",repel =TRUE)
dev.off ()
jpeg("PCA_variables.jpg", height = 6, width = 8, units = 'in', res = 600)
fviz_pca_var(res.pca, col.var = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
dev.off ()

###################
#------ Figure 5e
###################
######################################
#----- Cumultaive CaZy Family Graphs
######################################
write.table(Total_bins_cazy_final_relab, file = "Cumulative_Relative_abundace.txt", sep = "\t")
Cumul_CaZy <- read.csv(file="Cumulative_Total_rel.txt", sep = "\t", row.names = 1, header = T)

Colors <- c("darkgreen", "darkred")

jpeg("CMB_new.jpg", height = 4, width = 2, units = 'in', res = 300)
boxplot(Cumul_CaZy$CBM ~ Cumul_CaZy$Taxa, col=c("#CD853F", "#BDB76B"), ylab="Average Relative Abundance%", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8), outlier=FALSE)
stripchart(CBM ~ Taxa, vertical = TRUE, data = Cumul_CaZy, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off()
wilcox.test (Cumul_CaZy$CBM ~ Cumul_CaZy$Taxa)
jpeg("CE.jpg", height = 4, width = 2, units = 'in', res = 300)
boxplot(Cumul_CaZy$CE ~ Cumul_CaZy$Taxa, col=c("#CD853F", "#BDB76B"), ylab="Average Relative Abundance%", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(CE ~ Taxa, vertical = TRUE, data = Cumul_CaZy, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off()
wilcox.test (Cumul_CaZy$CE ~ Cumul_CaZy$Taxa)
jpeg("GH.jpg", height = 4, width = 2, units = 'in', res = 300)
boxplot(Cumul_CaZy$GH ~ Cumul_CaZy$Taxa, col=c("#CD853F", "#BDB76B"), ylab="Average Relative Abundance%", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(GH ~ Taxa, vertical = TRUE, data = Cumul_CaZy, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off()
wilcox.test (Cumul_CaZy$GH ~ Cumul_CaZy$Taxa)
jpeg("GT.jpg", height = 4, width = 2, units = 'in', res = 300)
boxplot(Cumul_CaZy$GT ~ Cumul_CaZy$Taxa, col=c("#CD853F", "#BDB76B"), ylab="Average Relative Abundance%", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(GT ~ Taxa, vertical = TRUE, data = Cumul_CaZy, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off()
wilcox.test (Cumul_CaZy$GT ~ Cumul_CaZy$Taxa)
#W = 20, p-value = 0.1023

###################
#------ Figure 5f
###################

#----- CaZy Analysis ot get the significantly discriminating functions (For HeatMap)
#-- Assign functions to 
#------- To show the variables across samples
Bins_relative_abun_meta <- cbind(Total_bins_cazy_final_relab, meta$Taxa)
source('~/Work/Metagenomic_Analysis/Final_Analysis/GenomeReconstuction/Total_Prevo_SpiroBins/WilcoxPairwise.R')
Wilcoxon_pairwise(Bins_relative_abun_meta)

Selected_cazy <- read.csv(file="Selected_CaZy_annotated.txt", sep = "\t", row.names = 1, header =T)
Selected_cazy_matrix <- as.matrix(Selected_cazy)
colors = c(seq(-3,-2,length=100),seq(-2,0.5,length=100),seq(0.5,6,length=100))
my_palette <- colorRampPalette(c("gray55", "pink", "darkgreen"))(n = 299)

jpeg("CaZy_heatMap.jpg", height = 6, width = 6, units = 'in', res = 600)
aheatmap(sqrt(Selected_cazy_matrix), scale = "col", col=my_palette, distfun = "spearman", hclustfun = "complete", fontsize=8)
dev.off ()

jpeg("CaZy_heatMap_new.jpg", height = 8, width = 10, units = 'in', res = 300)
heatmap.2(Selected_cazy_matrix, scale = "column", col=my_palette, distfun = "spearman", hclustfun = "complete", fontsize=10)
dev.off ()

############################
#--- Supplementary Figures
############################
#-------------------------------------------------------------
# FigureS7 and S8. Genome Reconstruction and phylogenetic tree
#-------------------------------------------------------------
bins <- read.csv(file="Bins_gorilla_com50_cont10.txt", sep = "\t", row.names = 1, header =T)
bins[which(bins['Completeness'] >= 50 & abs(bins['Contamination']) <= 4 ),"group"] <- "Partial"
bins[which(bins['Completeness'] >= 70 & abs(bins['Contamination']) <= 10 ),"group"] <- "Medium-complete"
bins[which(bins['Completeness'] >= 90 & abs(bins['Contamination']) <= 5 ),"group"] <- "Near-complete"
bins_naremoved <- na.omit(bins)

jpeg("GenomeQualityfromGorilla.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(bins_naremoved, aes(x=Completeness, y=Contamination)) + geom_point(size = 0.5) + geom_point(aes(colour = group), size = 2.5) + scale_color_manual(values=Colors) + theme(axis.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), axis.title = element_text(size = 15, face = "bold")) + labs(x="Completeness (%)", y="Contamination (%)") + theme_bw() + ggtitle("")
dev.off ()
jpeg("GenomeQualityfromGorilla_scaffolds.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(bins_naremoved, aes(x=Completeness, y=Contigs)) + geom_point(size = 0.5) + geom_point(aes(colour = group), size = 2.5) + scale_color_manual(values=Colors) + theme(axis.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), axis.title = element_text(size = 15, face = "bold")) + labs(x="Completeness (%)", y="Number of scaffolds") + theme_bw() + ggtitle("")
dev.off ()

bins <- read.csv(file="Bins_human_comp50_cont10.txt", sep = "\t", row.names = 1, header =T)
bins[which(bins['Completeness'] >= 50 & abs(bins['Contamination']) <= 4 ),"group"] <- "Partial"
bins[which(bins['Completeness'] >= 70 & abs(bins['Contamination']) <= 10 ),"group"] <- "Medium-complete"
bins[which(bins['Completeness'] >= 90 & abs(bins['Contamination']) <= 5 ),"group"] <- "Near-complete"
bins_naremoved <- na.omit(bins)

jpeg("GenomeQualityfromHuman.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(bins_naremoved, aes(x=Completeness, y=Contamination)) + geom_point(size = 0.5) + geom_point(aes(colour = group), size = 2.5) + scale_color_manual(values=Colors) + theme(axis.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), axis.title = element_text(size = 15, face = "bold")) + labs(x="Completeness (%)", y="Contamination (%)") + theme_bw() + ggtitle("")
dev.off ()
jpeg("GenomeQualityfromHuman_scaffolds.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(bins_naremoved, aes(x=Completeness, y=Contigs)) + geom_point(size = 0.5) + geom_point(aes(colour = group), size = 2.5) + scale_color_manual(values=Colors) + theme(axis.text = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14, face = "bold"), legend.title = element_text(size = 16, face = "bold"), axis.title = element_text(size = 15, face = "bold")) + labs(x="Completeness (%)", y="Number of scaffolds") + theme_bw() + ggtitle("")
dev.off ()

# To make alignment for Prevotella bins
#~/Softwares/ezTree-master/ezTree -list Prevo_listSpeices -out genome.out -thread 120 -evalue 1e-5
# To make alignment for Spirochaetaceae bins
#~/Softwares/ezTree-master/ezTree -list list_species -out genome.out -thread 120 -evalue 1e-5

#Evolgenius was used to create the Phylogenetic Trees
#http://www.evolgenius.info/evolview/#mytrees/Meta/Trepo

#----------------------------------
# FigureS9 - Bins distribution
#----------------------------------
bins_distribution <- read.csv (file="Bins_depth_distribute.txt", sep="\t", row.names=1, header=T)
bins_distribution_t <- data.frame(t(bins_distribution))
bin_proportions <- bins_distribution_t/colSums(bins_distribution_t)[col(bins_distribution_t)]*100
bin_distribution_proportions <- data.frame(t(bin_proportions))
row.names (bin_distribution_proportions)
meta <- read.csv (file = "Final_Data/Metadata.txt", sep="\t", row.names=1, header=T)
row.names (meta)

jpeg("Bin148_new.jpg", height = 3, width = 2, units = 'in', res = 600)
boxplot(sqrt(bin_distribution_proportions$Bin148) ~ meta$Group1, col=c("skyblue4", "skyblue", "wheat4", "wheat"), ylab="Relative depth of Bin148", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8), ylim = c(0, 4.5))
stripchart(sqrt(bin_distribution_proportions$Bin148) ~ meta$Group1, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()

jpeg("Bin487_new2.jpg", height = 3, width = 2, units = 'in', res = 600)
boxplot(sqrt(sqrt(bin_distribution_proportions$Bin487)) ~ meta$Group1, col=c("skyblue4", "skyblue", "wheat4", "wheat"), ylab="Relative depth of Bin487", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8), ylim = c(0, 3))
stripchart(sqrt(sqrt(bin_distribution_proportions$Bin487)) ~ meta$Group1, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()

#-------- Final Distribution Graph which was Used
bin_distribution_proportions_meta <- cbind (bin_distribution_proportions, meta)
#write.table(bin_distribution_proportions_meta, file="bin_distribution_proportions_meta.txt", sep="\t")

bin_distribution_proportions_meta <- read.csv(file="Final_Data/bin_distribution_proportions_meta.txt",sep="\t",row.names=1,header=T)

Group1 <- as.factor(bin_distribution_proportions_meta$Group1)
jpeg("Bin148_distribution_New.jpg", height = 5, width = 7, units = 'in', res = 600)
barplot(sqrt(bin_distribution_proportions_meta$Bin148), col=c("skyblue4", "skyblue", "wheat4", "wheat")[Group1])
dev.off ()
jpeg("Bin487_distribution_New.jpg", height = 5, width = 7, units = 'in', res = 600)
barplot(sqrt(sqrt(bin_distribution_proportions_meta$Bin487)), col=c("skyblue4", "skyblue", "wheat4", "wheat")[Group1])
dev.off ()


#---- Pvalues for Gorilla and Human
bin_distribution_human <- subset(bin_distribution_proportions_meta, Group=="Human")
wilcox.test(bin_distribution_human$Bin148 ~ bin_distribution_human$Group1)
wilcox.test(bin_distribution_human$Bin487 ~ bin_distribution_human$Group1)

bin_distribution_gorilla <- subset(bin_distribution_proportions_meta, Group=="Gorilla")
wilcox.test(bin_distribution_gorilla$Bin148 ~ bin_distribution_gorilla$Group1)
wilcox.test(bin_distribution_gorilla$Bin487 ~ bin_distribution_gorilla$Group1)


barplot(bin_distribution_proportions_meta$Bin148, col=c("brown4", "orangered3", "green4", "limegreen")[Group1])

