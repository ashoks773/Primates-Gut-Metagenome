#Step2: PCA Analysis using GENE PROPORTION TABLE
Gene <- read.csv(file = "../GENE_PROPORTION_TABLE", row.names=1, header = T, sep = "\t")
metadata <- read.csv (file = "../Metadata.txt", row.names=1, header = T, sep = "\t")
Gene1 <- data.frame(t(Gene))

#--- Method Bray
bray_dist <-vegdist(Gene1, "bray")
bray_dist <- as.matrix (bray_dist)
m2 <- melt(bray_dist)
write.table(m2, file="Bray_Distances.txt", sep = "\t")

bray_pcoa <-pcoa(vegdist(Gene1, "bray"))
Groups <- metadata[,1:3]
bray_pcoa$values[1:2,]
mds.var.per = round(bray_pcoa$values$Eigenvalues/sum(bray_pcoa$values$Eigenvalues)*100, 1)
Bray_PCoA_MATRIX <- bray_pcoa$vectors[,1:2]
Bray_PCoA_MATRIX <- data.frame(Bray_PCoA_MATRIX)
Bray_PCoA_MATRIX_New <- cbind(Bray_PCoA_MATRIX, Groups)

Color <- c("darkgreen", "darkred")
Colors <- c("brown4", "orangered3", "green4", "limegreen")

jpeg("Bray_PCoA_Group.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
jpeg("Bray_PCoA_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Colors) + xlab(paste("PCoA1 - ", mds.var.per[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

bray_dist <-vegdist(Gene1, "bray")
adonis(bray_dist ~ Group, data=metadata)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Group      1    5.6899  5.6899  20.816 0.29816   0.01 **
  
adonis(bray_dist ~ Group1, data=metadata)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Group1     3    7.3059 2.43530  9.7184 0.38284   0.01 **
#Residuals 47   11.7776 0.25059         0.61716  


pcoa1 <- Bray_PCoA_MATRIX_New[,c(1,4)]
pcoa1_Melted <- melt(pcoa1, id.vars = "Group1")
jpeg("Bray_PCoA1_Distances.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = pcoa1_Melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA1") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + coord_flip() + theme(legend.position='none')
dev.off ()
pcoa2 <- Bray_PCoA_MATRIX_New[,c(2,4)]
pcoa2_Melted <- melt(pcoa2, id.vars = "Group1")
jpeg("Bray_PCoA2_Distances.jpg", height = 4, width = 4, units = 'in', res = 600)
ggplot(data = pcoa2_Melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("Bray-Curtis Distances") + labs(x="",y="PCoA2") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()

#---Gorilla only
metadata_gorilla <- read.csv (file = "../Metadata_Gorilla.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_gorilla[,2:3]
gene_gorilla <- Gene1[1:23,]
bray_pcoa_gorilla <-pcoa(vegdist(gene_gorilla, "bray"))
bray_pcoa_gorilla$values[1:2,]
mds.var.per_g = round(bray_pcoa_gorilla$values$Eigenvalues/sum(bray_pcoa_gorilla$values$Eigenvalues)*100, 1)
Bray_PCoA_MATRIX_G <- bray_pcoa_gorilla$vectors[,1:2]
Bray_PCoA_MATRIX_G <- data.frame(Bray_PCoA_MATRIX_G)
Bray_PCoA_MATRIX_G_New <- cbind(Bray_PCoA_MATRIX_G, Groups)

Color <- c("green4", "limegreen")
jpeg("Bray_PCoA_Gorilla_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()
jpeg("Bray_PCoA_Gorilla_Group2.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_G_New, aes(x=Axis.1, y=Axis.2, colour=Group1, shape=Group2)) + geom_point(size=3) + geom_point(aes(shape=Group2)) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_g[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_g[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted") + geom_point(shape=18)
dev.off ()

bray_dist_gorilla <-vegdist(gene_gorilla, "bray")
adonis(bray_dist_gorilla ~ Group1 + Group2, data=metadata_gorilla, permutations=99)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Group1     1    0.7022 0.70225  3.8418 0.14601   0.01 **
#Group2     1    0.4517 0.45169  2.4711 0.09391   0.01 **

#---- Human Only
metadata_human <- read.csv (file = "../Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_human[,2:3]
gene_human <- Gene1[24:51,]
bray_pcoa_human <-pcoa(vegdist(gene_human, "bray"))
bray_pcoa_human$values[1:2,]
mds.var.per_h = round(bray_pcoa_human$values$Eigenvalues/sum(bray_pcoa_human$values$Eigenvalues)*100, 1)
Bray_PCoA_MATRIX_H <- bray_pcoa_human$vectors[,1:2]
Bray_PCoA_MATRIX_H <- data.frame(Bray_PCoA_MATRIX_H)
Bray_PCoA_MATRIX_H_New <- cbind(Bray_PCoA_MATRIX_H, Groups)

Color <- c("brown4", "orangered3")
jpeg("Bray_PCoA_Human_Group1.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(Bray_PCoA_MATRIX_H_New, aes(x=Axis.1, y=Axis.2, colour=Group1)) + geom_point(size=3) + scale_shape_manual(values =1:20) + scale_color_manual(values=Color) + xlab(paste("PCoA1 - ", mds.var.per_h[1], "%", sep="")) + ylab(paste("PCoA2 - ", mds.var.per_h[2], "%", sep="")) + ggtitle("") + theme(axis.text.x = element_text(size = 12, colour = "black", face = "bold"), axis.text.y = element_text(size = 12, colour = "black", face = "bold"), legend.text = element_text(size = 14, colour = "black"), legend.title = element_text(size = 16, face = "bold")) + theme_bw() + geom_vline(xintercept = 0, linetype="dotted") +  geom_hline(yintercept = 0, linetype="dotted")
dev.off ()

bray_dist_human <-vegdist(gene_human, "bray")
adonis2(bray_dist_human ~ Group1, data=metadata_human, permutations=99)
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
#Group1     1    0.9138 0.91377  3.0975 0.10645   0.01 **