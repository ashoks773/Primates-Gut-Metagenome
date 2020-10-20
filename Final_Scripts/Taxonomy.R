# Taxonomic analysis using NCBI and HMP (Steps I have to take from Vishnu)
setwd("~/Work/Metagenomic_Analysis/Final_Analysis/TaxonomicAnalysis")
taxa <- read.csv (file = "NCBI_HMP_TAXA_ABUNDANCE", sep = "\t", header = T, row.names = 1)
taxa_proportions <- taxa/colSums(taxa)[col(taxa)]
write.table(taxa_proportions, file ="NCBI_HMP_TAXA_proportions.txt", sep = "\t")
taxa1 <- data.frame(t(taxa_proportions))

metadata <- read.csv(file="../Metadata.txt", sep = "\t", row.names=1, header = T)

bray_dist <- vegdist(taxa1, "bray")
taxa_bray_pcoa <-pcoa(bray_dist)
Groups <- metadata[,1:3]
taxa_bray_pcoa$values[1:2,]
mds.var.per = round(taxa_bray_pcoa$values$Eigenvalues/sum(taxa_bray_pcoa$values$Eigenvalues)*100, 1)
Taxa_Bray_PCoA_MATRIX <- taxa_bray_pcoa$vectors[,1:2]
Taxa_Bray_PCoA_MATRIX <- data.frame(Taxa_Bray_PCoA_MATRIX)
Taxa_Bray_PCoA_MATRIX_New <- cbind(Taxa_Bray_PCoA_MATRIX, Groups)

Color <- c("darkgreen", "darkred")
Colors <- c("brown4", "orangered3", "green4", "limegreen")

jpeg("Taxa_Bray_PCoA_Group.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Taxa_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
jpeg("Taxa_Bray_PCoA_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Taxa_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCoA1 - ", mds.var.per[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
bray_dist <- as.matrix(bray_dist)
adonis (bray_dist ~ Group1, data=metadata)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Group1     3    4.2454 1.41515  20.359 0.56513  0.001 ***


#------------Gorilla only
metadata_gorilla <- read.csv (file = "../Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
taxa_gorilla <- taxa1[1:23,]
taxa_bray_pcoa_gorilla <-pcoa(vegdist(taxa_gorilla, "bray"))
taxa_bray_pcoa_gorilla$values[1:2,]
mds.var.per_g = round(taxa_bray_pcoa_gorilla$values$Eigenvalues/sum(taxa_bray_pcoa_gorilla$values$Eigenvalues)*100, 1)
Taxa_Bray_PCoA_MATRIX_G <- taxa_bray_pcoa_gorilla$vectors[,1:2]
Taxa_Bray_PCoA_MATRIX_G <- data.frame(Taxa_Bray_PCoA_MATRIX_G)
Taxa_Bray_PCoA_MATRIX_G_New <- cbind(Taxa_Bray_PCoA_MATRIX_G, Groups)

Color <- c("green4", "limegreen")
jpeg("Bray_PCoA_Gorilla_Group1.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(Taxa_Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=4) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
bray_dist_gorilla <-vegdist(taxa_gorilla, "bray")
adonis (bray_dist_gorilla ~ Group1, data=metadata_gorilla)
#Group1     1   0.12489 0.124886  3.8257 0.1541  0.001 ***
  
setwd("~/Work/Metagenomic_Analysis/Final_Analysis/TaxonomicAnalysis")
taxa <- read.csv (file = "../NCBI_HMP_TAXA_ABUNDANCE", sep = "\t", header = T, row.names = 1)
taxa_proportions <- taxa/colSums(taxa)[col(taxa)]
taxa1 <- data.frame(t(taxa_proportions))
metadata <- read.csv(file="Metadata_twoGroup.txt", sep = "\t", row.names=1, header = T)
Group <- metadata$Group3
taxa_proportions_Group <- cbind(taxa1, Group)

#-- Indicator species Labdsv
library (labdsv)
#taxa_formatted <- taxa_proportions_Group[ , which(!apply(taxa_proportions_Group==0,2,all))]
#OTU_formatted[is.na(OTU_formatted)] <- 0
iva <- indval(taxa_proportions_Group[,1:810], taxa_proportions_Group$Group)
gr <- iva$maxcls
iv <- iva$indcls
pv <- iva$pval
#fr <- apply(taxa_proportions_Group[,1:810]>0, 2, sum)[iva$pval<=0.05]
indvalsummary <- data.frame(group=gr, indval=iv, pvalue=pv)
#indvalsummary <- indvalsummary[order(indvalsummary$group, -indvalsummary$indval),]
write.table (indvalsummary, file="indvalsummary_New.txt", sep = "\t")

source('~/Work/Metagenomic_Analysis/Final_Analysis/TaxonomicAnalysis/TaxaSelectionTwoGroups/Odds_Ratio_Calculation.R')
Calculate_Odds_Ratio(taxa_proportions_Group)

#---- To calculate fold change in Human
taxa_human <- taxa1[24:51,]
metadata_human <- read.csv(file="../../Metadata_Human.txt", sep = "\t", row.names=1, header = T)
Group1 <- metadata_human$Group1
taxa_human_Group <- cbind (taxa_human, Group1)
source('~/Work/Metagenomic_Analysis/Final_Analysis/TaxonomicAnalysis/NewTaxaSelection/Human_Odds_Ratio_Calculation.R')
Calculate_Odds_Ratio(taxa_human_Group)

#--------------------------- To plot selected Microbes
#https://jkzorz.github.io/2019/06/05/Bubble-plots.html
taxa_Selected <- read.csv (file = "Sel_taxa_proportions_Group.txt", sep = "\t", header = T, row.names = 1)
Samples <- row.names(taxa_Selected)
Group <- taxa_Selected$Group

taxa_Sel_per <- taxa_Selected[,2:28]*100
taxa_Sel_Group <- cbind(Samples, Group, taxa_Sel_per)

pcm = melt(taxa_Sel_Group, id = c("Samples", "Group"))

pcm$Sample <- factor(pcm$Sample,levels=unique(pcm$Sample))

jpeg("Bubble_plot_larger.New.png", height = 6, width = 6, units = 'in', res = 600)
ggplot(pcm, aes(x = Sample, y = variable)) + 
  geom_point(aes(size = value, fill = Group), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 80), range = c(1,10), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "Group")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 4, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", size = 7, face = "bold"), 
        legend.text = element_text(size = 7, colour ="black"), 
        legend.title = element_text(size = 7, face = "bold"), panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +  
  scale_fill_manual(values = c("brown4", "orangered3", "green4", "limegreen"), guide = guide_legend(override.aes = list(size=5))) + 
  scale_y_discrete(limits = rev(levels(pcm$variable))) 
dev.off ()


#------- AUC plot to check the classification accuracy

library (randomForest)
RF1 <- randomForest(Group1~., data=taxa_Selected, ntree = 10000, proximity = TRUE, importance = TRUE, do.trace = 1000, cv.fold = 10, na.action = na.omit)
jpeg("RF_MDA_compareAll4Groups.jpg", height = 12, width = 6, units = 'in', res = 600)
varImpPlot(RF1, type=1)
dev.off ()
MDA <- importance(RF1, type=1)
write.table (MDA, file = "MDA_comparisionbetweenAllFour.txt", sep = "\t")

jpeg("mda_check.jpg", height = 12, width = 6, units = 'in', res = 600)
varImpPlot(RF1, type=1, pch=23, lcolor="black")
dev.off ()

OOB.votes <- RF1$votes

# ------------------------ Sort Way ------------------------------------------------------------------------------------------------#

png ('ROC_multiclass.png', height = 4, width = 4, units = 'in', res = 600)
plot.roc(taxa_Selected$Group1, OOB.votes[,1], col='brown4', main = "ROC Curve for MUlticlass", lwd=3, identity.lty=4) #, auc.polygon=TRUE) #, auc.polygon.col="blue") # Area under the curve value will be the same at all the cut-off values but the 
plot.roc(taxa_Selected$Group1, OOB.votes[,2], add=TRUE, col='orangered3', main = "", lwd=3, identity.lty=4) #, auc.polygon=TRUE) #, auc.polygon.col="red" )# # Area under the curve value will be the same at all the cut-off values but the 
plot.roc(taxa_Selected$Group1, OOB.votes[,3], add=TRUE, col='green4', main = "", lwd=3, identity.lty=4)   # Area under the curve value will be the same at all the cut-off values but the 
plot.roc(taxa_Selected$Group1, OOB.votes[,4], add=TRUE, col='limegreen', main = "", lwd=3, identity.lty=4)    # Area under the curve value will be the same at all the cut-off values but the 
legend("bottomright", legend=c("BaAka", "Bantu", "Dry", "Wet"), col=c("brown4", "orangered3", "green4", "limegreen"), lwd=1, cex = 0.7)
dev.off ()
#plot.roc(taxa_Selected$Group1, OOB.votes[,2], add=TRUE, col='orangered3', main = "", lwd=3, identity.lty=4) #, auc.polygon=TRUE) #, auc.polygon.col="red" )# # Area under the curve value will be the same at all the cut-off values but the 
multiclass.roc(taxa_Selected$Group1, OOB.votes)