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

#Digestion of Mannan, Galactomann ...  (GH113)
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

#--- Not Used-----
#---- Selected Cazy from all six groups
Selected_cazy <- rbind(Animal_poly, Plant_poly, Algal_poly, Fungal_poly, Lignin_degradation, Biosynthesis_poly)
Selected_cazy <- data.frame(t(Selected_cazy))
Selected_cazy <- as.matrix(Selected_cazy)
Selected_cazy_sel <- subset( Selected_cazy, select = -GH27 ) # GH27 dropped because it was creating a bias in the analysis
categories <- metadata[,1:2]
rownames(categories) <- rownames(Selected_cazy_sel)
annot.color.col <- list('Group'=c('darkgreen', 'darkred'),'Group1'=c("brown4", "orangered3", "green4", "limegreen"))
aheatmap(sqrt(Selected_cazy_sel)*100, color = "-RdYlBu2:100", breaks = 0,  main = "Association", distfun = "spearman", hclustfun = "complete", fontsize=10,  filename="DesHeatMap_Euclidian_Names.png", scale = "row", Rowv = NA, Colv = NA, border_color = NA, legend = TRUE, cexRow = .7, cexCol = .7, annRow = categories, annCol = NA, annColors = annot.color.col, annLegend = TRUE)
# Below command we have used to color a heatmap Groups
#---------
#Transpose plot
Selected_cazy_sel <- data.frame(t(Selected_cazy_sel))
Selected_cazy_sel<- as.matrix(Selected_cazy_sel)
rownames(categories) <- colnames(Selected_cazy_sel)
annot.color.col <- list('Group'=c('darkgreen', 'darkred'),'Group1'=c("brown4", "orangered3", "green4", "limegreen"))
aheatmap(sqrt(Selected_cazy_sel)*100, color = "-RdYlBu2:100", breaks = 0,  main = "Association", distfun = "spearman", hclustfun = "complete", fontsize=10,  filename="CaZy_tried.png", scale = "column", Rowv = TRUE, Colv = TRUE, border_color = NA, legend = TRUE, cexRow = .7, cexCol = .7, annRow = categories, annCol = NA, annColors = annot.color.col, annLegend = TRUE)

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
