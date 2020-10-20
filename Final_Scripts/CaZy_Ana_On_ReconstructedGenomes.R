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

jpeg("Observed.jpg", height = 4, width = 2, units = 'in', res = 600)
boxplot(total_div$N.obs ~ total_div$Taxa, col=c("darkred", "darkgreen"), ylab="Observed CaZy", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(N.obs ~ Taxa, vertical = TRUE, data = total_div, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()
wilcox.test(total_div$N.obs ~ meta$Taxa)

jpeg("Shnnon.jpg", height = 4, width = 2, units = 'in', res = 600)
boxplot(total_div$Shannon ~ total_div$Taxa, col=c("darkred", "darkgreen"), ylab="Shannon's H",  cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(Shannon ~ Taxa, vertical = TRUE, data = total_div, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off ()
wilcox.test(total_div$Shannon ~ meta$Taxa)

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
fviz_pca_ind(res.pca,col.ind=groups,palette=c("darkred","darkgreen"),addEllipses=TRUE,ellipse.type="confidence",legend.title="Groups",repel =TRUE)
dev.off ()
jpeg("PCA_variables.jpg", height = 6, width = 8, units = 'in', res = 600)
fviz_pca_var(res.pca, col.var = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
dev.off ()


#----- Cumultaive Graph
write.table(Total_bins_cazy_final_relab, file = "Cumulative_Relative_abundace.txt", sep = "\t")
Cumul_CaZy <- read.csv(file="Cumulative_Total_rel.txt", sep = "\t", row.names = 1, header = T)

Colors <- c("darkgreen", "darkred")

jpeg("CMB.jpg", height = 4, width = 2, units = 'in', res = 300)
boxplot(Cumul_CaZy$CBM ~ Cumul_CaZy$Taxa, col=c("darkred", "darkgreen"), ylab="Average Relative Abundance%", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(CBM ~ Taxa, vertical = TRUE, data = Cumul_CaZy, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off()
jpeg("CE.jpg", height = 4, width = 2, units = 'in', res = 300)
boxplot(Cumul_CaZy$CE ~ Cumul_CaZy$Taxa, col=c("darkred", "darkgreen"), ylab="Average Relative Abundance%", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(CE ~ Taxa, vertical = TRUE, data = Cumul_CaZy, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off()
jpeg("GH.jpg", height = 4, width = 2, units = 'in', res = 300)
boxplot(Cumul_CaZy$GH ~ Cumul_CaZy$Taxa, col=c("darkred", "darkgreen"), ylab="Average Relative Abundance%", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(GH ~ Taxa, vertical = TRUE, data = Cumul_CaZy, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off()
jpeg("GT.jpg", height = 4, width = 2, units = 'in', res = 300)
boxplot(Cumul_CaZy$GT ~ Cumul_CaZy$Taxa, col=c("darkred", "darkgreen"), ylab="Average Relative Abundance%", cex.axis=0.8, par(cex.lab=0.8), par(cex.axis=0.8))
stripchart(GT ~ Taxa, vertical = TRUE, data = Cumul_CaZy, method = "jitter", add = TRUE, pch = 20, col = 'black')
dev.off()


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

