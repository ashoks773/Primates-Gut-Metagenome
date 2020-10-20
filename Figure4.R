# Author: Ashok Kumar Sharma

# Taxonomic analysis using NCBI and HMP (Steps I have to take from Vishnu)
setwd("~/Work/Metagenomic_Analysis/Final_Analysis/TaxonomicAnalysis")
taxa <- read.csv (file = "Final_Data/NCBI_HMP_TAXA_ABUNDANCE", sep = "\t", header = T, row.names = 1)
taxa_proportions <- taxa/colSums(taxa)[col(taxa)]
#write.table(taxa_proportions, file ="NCBI_HMP_TAXA_proportions.txt", sep = "\t")
taxa1 <- data.frame(t(taxa_proportions))

metadata <- read.csv(file="Final_Data//Metadata.txt", sep = "\t", row.names=1, header = T)

bray_dist <- vegdist(taxa1, "bray")
taxa_bray_pcoa <-pcoa(bray_dist)
Groups <- metadata[,1:3]
taxa_bray_pcoa$values[1:2,]
mds.var.per1 <- round(taxa_bray_pcoa$values$Rel_corr_eig[1]*100, 2)
mds.var.per2 <- round(taxa_bray_pcoa$values$Rel_corr_eig[2]*100, 2)
Taxa_Bray_PCoA_MATRIX <- taxa_bray_pcoa$vectors[,1:2]
Taxa_Bray_PCoA_MATRIX <- data.frame(Taxa_Bray_PCoA_MATRIX)
Taxa_Bray_PCoA_MATRIX_New <- cbind(Taxa_Bray_PCoA_MATRIX, Groups)

Color <- c("skyblue4", "wheat4")
Colors <- c("skyblue4", "skyblue", "wheat4", "wheat")

jpeg("Taxa_Bray_PCoA_Group1.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(Taxa_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCoA1 - ", mds.var.per1, "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per2, "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
bray_dist <- as.matrix(bray_dist)
adonis (bray_dist ~ Group1, data=metadata)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#Group1     3    4.2454 1.41515  20.359 0.56513  0.001 ***

#------- Plot Axis1 and Axis2

pcoa2 <- Taxa_Bray_PCoA_MATRIX_New[,c(2,4)]
pcoa2_Melted <- melt(pcoa2, id.vars = "Group1")
jpeg("Bray_PCoA2_Distances.jpg", height = 4, width = 2, units = 'in', res = 600)
ggplot(data = pcoa2_Melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCo2") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()

kruskalmc(pcoa2$Axis.2 ~ pcoa2$Group1)


#------------Gorilla only
metadata_gorilla <- read.csv (file = "Final_Data/Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
taxa_gorilla <- taxa1[1:23,]
taxa_bray_pcoa_gorilla <-pcoa(vegdist(taxa_gorilla, "bray"))
taxa_bray_pcoa_gorilla$values[1:2,]
mds.var.per_g1 <- round(taxa_bray_pcoa_gorilla$values$Rel_corr_eig[1]*100, 2)
mds.var.per_g2 <- round(taxa_bray_pcoa_gorilla$values$Rel_corr_eig[2]*100, 2)
#mds.var.per_g = round(taxa_bray_pcoa_gorilla$values$Eigenvalues/sum(taxa_bray_pcoa_gorilla$values$Eigenvalues)*100, 1)
Taxa_Bray_PCoA_MATRIX_G <- taxa_bray_pcoa_gorilla$vectors[,1:2]
Taxa_Bray_PCoA_MATRIX_G <- data.frame(Taxa_Bray_PCoA_MATRIX_G)
Taxa_Bray_PCoA_MATRIX_G_New <- cbind(Taxa_Bray_PCoA_MATRIX_G, Groups)

Color <- c("wheat4", "wheat")
jpeg("Bray_PCoA_Gorilla_Group1.jpg", height = 2.5, width = 3.5, units = 'in', res = 600)
ggplot(Taxa_Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=4) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g1, "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g2, "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
bray_dist_gorilla <-vegdist(taxa_gorilla, "bray")
adonis (bray_dist_gorilla ~ Group1, data=metadata_gorilla)
#Group1     1   0.12489 0.124886  3.8257 0.1541  0.001 ***
  
#-------------------------- This whole step is now repetaed in NewTaxaSelection Folder
#----Important Taxa Identification
#RandomForest
Group <- metadata$Group
taxa_proportions_Group <- cbind(taxa1, Group)
library (randomForest)
RF1 <- randomForest(Group~., data=taxa_proportions_Group, ntree = 1000, proximity = TRUE, importance = TRUE, do.trace = 100, cv.fold = 10, na.action = na.omit)
jpeg("RF_MDA_gorillavsHuman.jpg", height = 12, width = 6, units = 'in', res = 600)
varImpPlot(RF1, type=1)
dev.off ()
MDA <- importance(RF1, type=1)
write.table (MDA, file = "MDA_GorillavsHuman.txt", sep = "\t")

Group1 <- metadata$Group1
taxa_proportions_Group <- cbind(taxa1, Group1)
library (randomForest)
RF1 <- randomForest(Group1~., data=taxa_proportions_Group, ntree = 1000, proximity = TRUE, importance = TRUE, do.trace = 100, cv.fold = 10, na.action = na.omit)
jpeg("RF_MDA_compareAll4Groups.jpg", height = 12, width = 6, units = 'in', res = 600)
varImpPlot(RF1, type=1)
dev.off ()
MDA <- importance(RF1, type=1)
write.table (MDA, file = "MDA_comparisionbetweenAllFour.txt", sep = "\t")

#ntree      OOB      1      2      3      4
#100:   7.84%  0.00% 14.29%  9.09%  8.33%
#200:   5.88%  0.00% 14.29%  0.00%  8.33%
#300:   3.92%  0.00%  7.14%  0.00%  8.33%
#400:   3.92%  0.00%  7.14%  0.00%  8.33%
#500:   3.92%  0.00%  7.14%  0.00%  8.33%
#600:   3.92%  0.00%  7.14%  0.00%  8.33%
#700:   3.92%  0.00%  7.14%  0.00%  8.33%
#800:   3.92%  0.00%  7.14%  0.00%  8.33%
#900:   3.92%  0.00%  7.14%  0.00%  8.33%
#1000:   3.92%  0.00%  7.14%  0.00%  8.33%

#---- ROC Curve
OOB.votes <- RF1$votes
write.table (OOB.votes, file = "OOB_pred", sep = "\t")
predictions <- read.csv (file = "OOB_pred", sep = "\t")
roc.multi <- multiclass.roc(taxa_proportions_Group$Group1, predictions$Bantu)
auc(roc.multi)
# AUC: 0.869
rs <- roc.multi[['rocs']]
plot.roc(rs[[1]])
sapply(2:length(rs),function(i) lines.roc(rs[[i]],col=c("green")))

jpeg("AUC.png", height = 4, width = 4, units = 'in', res = 600)
plot.roc(rs[[1]], print.auc = TRUE, grid=c(0.1, 0.2), lwd=3)
dev.off ()
#RF1$confusion
#BaAka Bantu Dry Wet class.error
#BaAka    14     0   0   0  0.00000000
#Bantu     1    13   0   0  0.07142857
#Dry       0     0  11   0  0.00000000
#Wet       0     0   1  11  0.08333333

#---- Make Taxa Bubble Plot
#Rscript Taxa_bubbleNew.R