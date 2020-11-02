# Author: Ashok Kumar Sharma

# Step 1: Take only lipid metabolism pathways - from total pathway relative abundance table
Pathway <- read.csv(file = "Lipid_metabolism//Lipid_metabolism_raw.txt", row.names=1, header = T, sep = "\t")
Pathway <- read.csv(file = "check.txt", row.names=1, header = T, sep = "\t")

Pathway_proportions <- Pathway/colSums(Pathway)[col(Pathway)]
Pathway1 <- data.frame(t(Pathway_proportions))
metadata <- read.csv (file = "../../Metadata.txt", row.names=1, header = T, sep = "\t")

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

jpeg("Lipid_Group1_try.jpg", height = 3, width = 4, units = 'in', res = 600)
ggplot(Pathway_Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCoA1 - ", mds.var.per1, "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per2, "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

distance_bray <-vegdist(Pathway1, "bray")
adonis(distance_bray ~ Group*Group1, data=metadata, permutations=999)
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
#Group      1  0.050977 0.050977  37.087 0.37857   0.01 **
# Group1     2  0.019079 0.009540   6.940 0.14169   0.01 **

#---Gorilla only
metadata_gorilla <- read.csv (file = "../../Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
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
jpeg("Lipid_Gorilla_Group1.jpg", height = 2.5, width = 3.5, units = 'in', res = 600)
ggplot(Pathway_Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g1, "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g2, "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

dist_bray_gorilla <-vegdist(pathway_gorilla, "bray")
adonis(dist_bray_gorilla ~ Group1, data=metadata_gorilla, permutations=99)
#Df SumsOfSqs    MeanSqs F.Model      R2 Pr(>F)  
#Group2         1 0.0029329 0.00293295  3.8866 0.14474   0.03 *
# Group1         1 0.0021860 0.00218600  2.8968 0.10788   0.05 *

#--- Only Human
metadata_human <- read.csv (file = "../../Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_human[,2:3]
pathway_human <- Pathway1[24:51,]
pathway_bray_pcoa_human <-pcoa(vegdist(pathway_human, "bray"))
pathway_bray_pcoa_human$values[1:2,]
mds.var.per_g = round(pathway_bray_pcoa_human$values$Eigenvalues/sum(pathway_bray_pcoa_human$values$Eigenvalues)*100, 1)
Pathway_Bray_PCoA_MATRIX_H <- pathway_bray_pcoa_human$vectors[,1:2]
Pathway_Bray_PCoA_MATRIX_H <- data.frame(Pathway_Bray_PCoA_MATRIX_H)
Pathway_Bray_PCoA_MATRIX_H_New <- cbind(Pathway_Bray_PCoA_MATRIX_H, Groups)

Color <- c("skyblue4", "skyblue")
jpeg("Lipid_Human_Group1.jpg", height = 2.5, width = 3.5, units = 'in', res = 600)
ggplot(Pathway_Bray_PCoA_MATRIX_H_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

dist_bray_human <-vegdist(pathway_human, "bray")
adonis(dist_bray_human ~ Group1, data=metadata_human, permutations=99)

#-------- Only Lipid metabolism pathways were fetched from total Pathway relative abundances
#------ Individual BoxPlots for Lipid Metabolism for Gorilla
lipid <- read.csv(file = "Lipid_metabolism/Lipid_metabolism.txt", row.names=1, header = T, sep = "\t")
lipid1 <- data.frame(t(lipid))
path2 <- sqrt(lipid1)
path3 <- path2*100

metadata_gorilla <- read.csv (file = "Final_Data/Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
pathway_gorilla <- path3[1:23,]
Color <- c( "wheat4", "wheat")
Pathways_sqrt_per_new <- cbind(pathway_gorilla, Groups)

#--
ketone <- Pathways_sqrt_per_new[,c(1,12)]
ketone_melt <- melt(ketone, id.vars = "Group1")
jpeg("G_KetoneBodies.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = ketone_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Synthesis.and.degradation.of.ketone.bodies") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Synthesis.and.degradation.of.ketone.bodies ~ Pathways_sqrt_per_new$Group1)
#W = 22, p-value = 0.007408

Glycerophospho <- Pathways_sqrt_per_new[,c(2,12)]
Glycerophospho_melt <- melt(Glycerophospho, id.vars = "Group1")
jpeg("G_Glycerophospho.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Glycerophospho_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Glycerophospholipid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Glycerophospholipid.metabolism ~ Pathways_sqrt_per_new$Group1)

Glycerolipid <- Pathways_sqrt_per_new[,c(3,12)]
Glycerolipid_melt <- melt(Glycerolipid, id.vars = "Group1")
jpeg("G_Glycerolipid.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Glycerolipid_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Glycerolipid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Glycerolipid.metabolism ~ Pathways_sqrt_per_new$Group1)

Sphingolipid <- Pathways_sqrt_per_new[,c(4,12)]
Sphingolipid_melt <- melt(Sphingolipid, id.vars = "Group1")
jpeg("G_Sphingolipid.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Sphingolipid_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Sphingolipid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Sphingolipid.metabolism ~ Pathways_sqrt_per_new$Group1)

Arachidonic <- Pathways_sqrt_per_new[,c(5,12)]
Arachidonic_melt <- melt(Arachidonic, id.vars = "Group1")
jpeg("G_Arachidonic.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Arachidonic_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Arachidonic.acid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Arachidonic.acid.metabolism ~ Pathways_sqrt_per_new$Group1)

Linoleic <- Pathways_sqrt_per_new[,c(6,12)]
Linoleic_melt <- melt(Linoleic, id.vars = "Group1")
jpeg("G_Linoleic.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Linoleic_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Linoleic.acid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Linoleic.acid.metabolism ~ Pathways_sqrt_per_new$Group1, method='FDR')
#W = 34, p-value = 0.05122

Ether.lipid <- Pathways_sqrt_per_new[,c(7,12)]
Ether.lipid_melt <- melt(Ether.lipid, id.vars = "Group1")
jpeg("G_Ether.lipid.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Ether.lipid_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Ether.lipid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Ether.lipid.metabolism ~ Pathways_sqrt_per_new$Group1)

Fatty.acid.degradation <- Pathways_sqrt_per_new[,c(8,12)]
Fatty.acid.degradation_melt <- melt(Fatty.acid.degradation, id.vars = "Group1")
jpeg("G_fATTY_Degradation.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Fatty.acid.degradation_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Fatty.acid.degradation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Fatty.acid.degradation ~ Pathways_sqrt_per_new$Group1)

Fatty.acid.elongation <- Pathways_sqrt_per_new[,c(10,12)]
Fatty.acid.elongation_melt <- melt(Fatty.acid.elongation, id.vars = "Group1")
jpeg("G_fATTY_elongation.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Fatty.acid.degradation_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Fatty.acid.elongation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Fatty.acid.elongation ~ Pathways_sqrt_per_new$Group1)

Fatty.acid.metabolism <- Pathways_sqrt_per_new[,c(10,12)]
Fatty.acid.metabolism_melt <- melt(Fatty.acid.metabolism, id.vars = "Group1")
jpeg("G_fATTY_metabolism.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Fatty.acid.metabolism_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Fatty.acid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Fatty.acid.metabolism ~ Pathways_sqrt_per_new$Group1)


#------ Individual BoxPlots for Lipid Metabolism for Human
#lipid <- read.csv(file = "Lipid_metabolism.txt", row.names=1, header = T, sep = "\t")
#lipid1 <- data.frame(t(lipid))
#path2 <- sqrt(lipid1)
#path3 <- path2*100

metadata_human <- read.csv (file = "Final_Data/Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_human[,2:3]
pathway_human <- path3[24:51,]
Color <- c("skyblue4", "skyblue")

Pathways_sqrt_per_new <- cbind(pathway_human, Groups)

#--
ketone <- Pathways_sqrt_per_new[,c(1,12)]
ketone_melt <- melt(ketone, id.vars = "Group1")
jpeg("H_KetoneBodies.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = ketone_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Synthesis.and.degradation.of.ketone.bodies") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Synthesis.and.degradation.of.ketone.bodies ~ Pathways_sqrt_per_new$Group1)
#W = 19, p-value = 0.0001015

Glycerophospho <- Pathways_sqrt_per_new[,c(2,12)]
Glycerophospho_melt <- melt(Glycerophospho, id.vars = "Group1")
jpeg("H_Glycerophospho.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Glycerophospho_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Glycerophospholipid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Glycerophospholipid.metabolism ~ Pathways_sqrt_per_new$Group1)

Glycerolipid <- Pathways_sqrt_per_new[,c(3,12)]
Glycerolipid_melt <- melt(Glycerolipid, id.vars = "Group1")
jpeg("H_Glycerolipid.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Glycerolipid_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Glycerolipid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Glycerolipid.metabolism ~ Pathways_sqrt_per_new$Group1)

Sphingolipid <- Pathways_sqrt_per_new[,c(4,12)]
Sphingolipid_melt <- melt(Sphingolipid, id.vars = "Group1")
jpeg("H_Sphingolipid.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Sphingolipid_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Sphingolipid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Sphingolipid.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 157, p-value = 0.005761

Arachidonic <- Pathways_sqrt_per_new[,c(5,12)]
Arachidonic_melt <- melt(Arachidonic, id.vars = "Group1")
jpeg("H_Arachidonic.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Arachidonic_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Arachidonic.acid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Arachidonic.acid.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 169, p-value = 0.0006462

Linoleic <- Pathways_sqrt_per_new[,c(6,12)]
Linoleic_melt <- melt(Linoleic, id.vars = "Group1")
jpeg("H_Linoleic.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Linoleic_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Linoleic.acid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Linoleic.acid.metabolism ~ Pathways_sqrt_per_new$Group1, method='FDR')
#W = 168, p-value = 0.0007944

Ether.lipid <- Pathways_sqrt_per_new[,c(7,12)]
Ether.lipid_melt <- melt(Ether.lipid, id.vars = "Group1")
jpeg("H_Ether.lipid.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Ether.lipid_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Ether.lipid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Ether.lipid.metabolism ~ Pathways_sqrt_per_new$Group1)

Fatty.acid.degradation <- Pathways_sqrt_per_new[,c(8,12)]
Fatty.acid.degradation_melt <- melt(Fatty.acid.degradation, id.vars = "Group1")
jpeg("H_fATTY_Degradation.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Fatty.acid.degradation_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Fatty.acid.degradation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Fatty.acid.degradation ~ Pathways_sqrt_per_new$Group1)

Fatty.acid.elongation <- Pathways_sqrt_per_new[,c(10,12)]
Fatty.acid.elongation_melt <- melt(Fatty.acid.elongation, id.vars = "Group1")
jpeg("H_fATTY_elongation.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Fatty.acid.degradation_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Fatty.acid.elongation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Fatty.acid.elongation ~ Pathways_sqrt_per_new$Group1)

Fatty.acid.metabolism <- Pathways_sqrt_per_new[,c(10,12)]
Fatty.acid.metabolism_melt <- melt(Fatty.acid.metabolism, id.vars = "Group1")
jpeg("H_fATTY_metabolism.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Fatty.acid.metabolism_melt,aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Fatty.acid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Fatty.acid.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 155, p-value = 0.007866



########################################
#--------- Total diversity Correlations
########################################
metadata <- read.csv (file = "Final_Data/Metadata.txt", row.names=1, header = T, sep = "\t")

Carbo <- read.csv(file = "Final_Data/Fam_ABUNDANCE", row.names=1, header = T, sep = "\t")
Carbo <- data.frame(t(Carbo))
Carbo_totalSum <- rowSums(Carbo)
Carbo_totalSum <- sqrt(Carbo_totalSum)*100
Carbo_totalSum <- data.frame(Carbo_totalSum)

Carbo_Shannon <- diversity(Carbo)
Carbo_Simp <- diversity(Carbo, "invsimpson")
div_carbo <- data.frame(Carbo_Shannon, Carbo_Simp)

Fat <- read.csv(file = "Lipid_metabolism//Lipid_metabolism_raw.txt", row.names=1, header = T, sep = "\t")
Fat <- data.frame(t(Fat))
Fat_totalSum <- rowSums(Fat)
Fat_totalSum <- sqrt(Fat_totalSum)*100
Fat_totalSum <- data.frame(Fat_totalSum)

Fat_Shannon <- diversity(Fat)
Fat_Simp <- diversity(Fat, "invsimpson")
div_Fat <- data.frame(Fat_Shannon, Fat_Simp)

#--- Merge Carbo-diversity, Fat-Divesity, Ketone Bodies abundance (from path3, sqrt %) and Metadata
Complete_Table <- data.frame(div_carbo, div_Fat, path3, Carbo_totalSum, Fat_totalSum, metadata)
Complete_Table_filtered <- Complete_Table[,c(1:5,16:19)]

Group1 <- as.factor (Complete_Table_filtered$Group1)

#####################
#---- Carbo Diversity
#####################
jpeg("Carbo_Keto_New.jpg", height = 2.5, width = 3.5, units = 'in', res = 600)
ggplot(data = Complete_Table_filtered, mapping = aes(x = Carbo_Shannon, y = Synthesis.and.degradation.of.ketone.bodies)) + 
  geom_point(mapping = aes(color = Group1)) + scale_color_manual(values=Colors) + 
  theme_bw() +
  xlab("Diversity (Carbohydrate metabolism)") + ylab("Sqrt abundance (keone bodies) %") +
  geom_smooth(method=lm, color="darkgray", fill="gray")
dev.off ()
Correlation <- corr.test(Complete_Table_filtered$Carbo_Shannon, Complete_Table_filtered$Synthesis.and.degradation.of.ketone.bodies, method="spearman", adjust="fdr", alpha =0.5, use = "pairwise")

##########################################
#---- Carbo Cumulative relative abundance
##########################################
jpeg("Carbo_Keto_Relab.jpg", height = 2.5, width = 3.5, units = 'in', res = 600)
ggplot(data = Complete_Table_filtered, mapping = aes(x = Carbo_totalSum, y = Synthesis.and.degradation.of.ketone.bodies)) + 
  geom_point(mapping = aes(color = Group1)) + scale_color_manual(values=Colors) + 
  theme_bw() +
  xlab("Sqrt abundance (carbo meta) %") + ylab("Sqrt abundance (keone bodies) %") +
  geom_smooth(method=lm, color="darkgray", fill="gray")
dev.off ()
Correlation <- corr.test(Complete_Table_filtered$Carbo_totalSum, Complete_Table_filtered$Synthesis.and.degradation.of.ketone.bodies, method="spearman", adjust="fdr", alpha =0.5, use = "pairwise")
