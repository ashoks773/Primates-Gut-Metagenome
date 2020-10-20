#New Results

metadata_human <- read.csv (file = "../Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Colors <- c("brown4", "orangered3")

#Digestion of animal Polysaccharides (GH13_32, GH13_23, GH13_41, GH13_16, PL13, PL21_1)
Animal_poly <- read.csv(file="Animal_poly", sep = "\t", row.names = 1, header =T)
Animal <- colSums(Animal_poly)
Animal <- data.frame(Animal)
Animal <- Animal[24:51,]
Animal_Group <- cbind(Animal, metadata_human)
jpeg("Animal_h.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Animal_Group, aes(x=Group1, y=sqrt(Animal)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of animal polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Animal ~ Group1, data = Animal_Group)
#W = 106, p-value = 0.7345

#Digestion of plant Polysaccharides (GH27, GH43_8, PL4_2)
Plant_poly <- read.csv(file="plant_poly", sep = "\t", row.names = 1, header =T)
Plant <- colSums(Plant_poly)
Plant <- data.frame(Plant)
Plant <- Plant[24:51,]
Plant_Group <- cbind(Plant, metadata_human)
jpeg("Plant_h.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Plant_Group, aes(x=Group1, y=sqrt(Plant)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of plant polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Plant ~ Group1, data = Plant_Group)
#W = 183, p-value = 1.86e-05

#Digestion of fungal Polysaccharides (GH113)
Fungal_poly <- read.csv(file="Fungal_poly", sep = "\t", row.names = 1, header =T)
Fungal <- colSums(Fungal_poly)
Fungal <- data.frame(Fungal)
Fungal <- Fungal[24:51,]
Fungal_Group <- cbind(Fungal, metadata_human)
jpeg("Fungal_h.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Fungal_Group, aes(x=Group1, y=sqrt(Fungal)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of fungal polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Fungal ~ Group1, data = Fungal_Group)
#W = 99, p-value = 0.982

#Digestion of Algal Polysaccharides (GH96)
Algal_poly <- read.csv(file="Algal_poly", sep = "\t", row.names = 1, header =T)
Algal <- colSums(Algal_poly)
Algal <- data.frame(Algal)
Algal <- Algal[24:51,]
Algal_Group <- cbind(Algal, metadata_human)
jpeg("Algal_h.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Algal_Group, aes(x=Group1, y=sqrt(Algal)*100, fill=Group1)) + geom_boxplot() + ggtitle("Digetion of algal polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Algal ~ Group1, data = Algal_Group)
#W = 112, p-value = 0.1649

#-- Lignin degradation (AA1, AA3_2, AA7)
Lignin_degradation <- read.csv(file="Lignin_degradation", sep = "\t", row.names = 1, header =T)
Lignin <- colSums(Lignin_degradation)
Lignin <- data.frame(Lignin)
Lignin <- Lignin[24:51,]
Lignin_Group <- cbind(Lignin, metadata_human)
jpeg("Lignin_h.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Lignin_Group, aes(x=Group1, y=sqrt(Lignin)*100, fill=Group1)) + geom_boxplot() + ggtitle("Lignin degradation") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Lignin ~ Group1, data = Lignin_Group)
#W = 132, p-value = 0.1227

#--#Biosynthesis of Polysaccharides (GT24, GT88, GT39, GT64, GT10, GT21)
Biosynthesis_poly <- read.csv(file="Biosynthesis_poly", sep = "\t", row.names = 1, header =T)
Biosynthesis <- colSums(Biosynthesis_poly)
Biosynthesis <- data.frame(Biosynthesis)
Biosynthesis <- Biosynthesis[24:51,]
Biosynthesis_Group <- cbind(Biosynthesis, metadata_human)
jpeg("Biosynthesis_h.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Biosynthesis_Group, aes(x=Group1, y=sqrt(Biosynthesis)*100, fill=Group1)) + geom_boxplot() + ggtitle("Biosynthesis of polysaccharides") + labs(x="",y="sqrt (Cumulative abundance)%") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Biosynthesis ~ Group1, data = Biosynthesis_Group) 
#W = 65, p-value = 0.1371

