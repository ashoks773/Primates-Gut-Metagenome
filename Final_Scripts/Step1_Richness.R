#----Figure1: Gene Richness and diversity analsyis

#---Gene Count file used for the rarefaction and diversity analysis
Abundance <- read.csv(file = "../LargeFiles_Used/Gene_count_matrix.tab", row.names=1, header = T, sep = "\t")
colSums(Abundance)
new <- colSums(Abundance)
new <- data.frame(new)
View(new)
# The lowest value is 4498117 and hence, we will rarefy the data for 4.4 million (4400000) This is in one
# Sample named as 2011_60

#biom convert -i Gene_count_matrix.tab -o Gene_count_matrix.tab.biom --table-type="OTU table" --to-json
#multiple_rarefactions_even_depth.py -i Gene_count_matrix.tab.biom -o rarefied_files -d 4400000 -n 10
#alpha_diversity.py -i rarefied_files -m shannon,simpson,chao1,Observed -o alpha_diversity_out 

#perl calculate_mean.pl
#perl addLabels.pl ../../Metadata.txt Alpha_diversity_mean

Color <- c("darkgreen", "darkred")
Colors <- c("brown4", "orangered3", "green4", "limegreen")

setwd("~/Work/Metagenomic_Analysis/GeneRichness/Alpha_diversity_at_4.4mdepth_Using_Genecountsmatrix")

alpha_div_mean <- read.csv (file="Alpha_diversity_mean.Out", sep = "\t", header = TRUE, row.names = 1)

#alpha_div_obsreved <- alpha_div_mean[, c(4,6)]
#alpha_div_obsreved_melted <- melt(alpha_div_obsreved, id.vars = "Group")
#jpeg("observed.jpg", height = 4, width = 4, units = 'in', res = 600)
#ggplot(data = alpha_div_obsreved_melted, aes(x=reorder(Group, value, FUN=median), y=value, fill=Group)) + geom_boxplot() + ggtitle("Gene Richness") + labs(x="",y="Observed") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))
#dev.off ()
#wilcox.test(alpha_div_obsreved$Observed_Genes ~ alpha_div_obsreved$Group1)

alpha_div_obsreved1 <- alpha_div_mean[, c(4,7)]
alpha_div_obsreved1 <- melt(alpha_div_obsreved1, id.vars = "Group1")
jpeg("observed1.jpg", height = 3.5, width = 1.8, units = 'in', res = 600)
ggplot(data = alpha_div_obsreved1, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("") + labs(x="",y="Observed genes") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
kruskalmc(alpha_div_obsreved1$value ~ alpha_div_obsreved1$Group1)

#https://stats.idre.ucla.edu/r/faq/how-can-i-do-post-hoc-pairwise-comparisons-in-r/
#Pairwise T-test
pairwise.t.test(alpha_div_obsreved1$value, alpha_div_obsreved1$Group1, p.adj = "bonf")
#BaAka  Bantu  Dry   
#Bantu 0.4101 -      -     
#  Dry   0.6068 0.0078 -     
#  Wet   0.4959 0.0051 1.0000

#-- Trukey test
a1 <- aov(alpha_div_obsreved1$value ~ alpha_div_obsreved1$Group1) 
TukeyHSD(a1)
#$`alpha_div_obsreved1$Group1`
#diff        lwr       upr     p adj
#Bantu-BaAka -89736.686 -217847.26  38373.89 0.2566149
#Dry-BaAka    85743.184  -50823.15 222309.52 0.3495620
#Wet-BaAka    88782.035  -44559.68 222123.75 0.2987838
#Dry-Bantu   175479.869   38913.54 312046.20 0.0068633
#Wet-Bantu   178518.720   45177.00 311860.44 0.0045472
#Wet-Dry       3038.851 -138446.33 144524.03 0.9999318

#jpeg("observed1.jpg", height = 4, width = 4, units = 'in', res = 600)
#ggplot(data = alpha_div_obsreved1_melted, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Gene Richness") + labs(x="",y="Observed") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))
#dev.off ()
#stats <- kruskalmc(alpha_div_obsreved1$Observed_Genes~alpha_div_obsreved1$Group1, probs=0.05)
##Stats_p0.05 <- stats$dif.com
#write.table(Stats_p0.05, file = "StatsObservedGroup1_p0.05.txt", sep= "\t")

#--- Shannon
#alpha_div_shannon <- alpha_div_mean[, c(1,6)]
#alpha_div_shannon_melted <- melt(alpha_div_shannon, id.vars = "Group")
#jpeg("shannon_new.jpg", height = 3, width = 1, units = 'in', res = 600)
#ggplot(data = alpha_div_shannon_melted, aes(x=reorder(Group, value, FUN=median), y=value, fill=Group)) + geom_boxplot() + ggtitle("Gene Diversity") + labs(x="",y="Shannon") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
#d#ev.off ()
#Wilcoxon_pairwise(alpha_div_shannon)

alpha_div_shannon1 <- alpha_div_mean[, c(1,7)]
alpha_div_shannon1_melted <- melt(alpha_div_shannon1, id.vars = "Group1")

jpeg("shannon1_new.jpg", height = 3.5, width = 1.8, units = 'in', res = 600)
ggplot(data = alpha_div_shannon1_melted, aes(x=Group1, y=value, fill=Group1)) + geom_boxplot() + ggtitle("") + labs(x="",y="Shannon") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black")) + theme(legend.position='none') + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off ()
kruskalmc(alpha_div_shannon1_melted$value ~ alpha_div_shannon1_melted$Group1)

#jpeg("shannon1_new.jpg", height = 3, width = 1.5, units = 'in', res = 600)
#ggplot(data = alpha_div_shannon1_melted, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Gene Diversity") + labs(x="",y="Shannon") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
#dev.off ()
#stats <- kruskalmc(alpha_div_shannon1$Shannon~alpha_div_shannon1$Group1, probs=0.05)
#Stats_p0.05 <- stats$dif.com
#write.table(Stats_p0.05, file = "StatsShannonGroup1_p0.05.txt", sep= "\t")












#----
Colors <- c("green4", "limegreen")
alpha_div_gorilla_mean <- read.csv (file="Alpha_diversity_mean_gorilla.out", sep = "\t", header = TRUE, row.names = 1)

gorilla_obsreved <- alpha_div_gorilla_mean[, c(4,7,8)]
jpeg("observed3.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(aes(y = Observed_Genes, x = Group2, fill = Group1), data = gorilla_obsreved) + geom_boxplot() + ggtitle("Gene Richness") + labs(x="",y="Observed") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))
dev.off ()
gorilla_shannon <- alpha_div_gorilla_mean[, c(1,7,8)]
jpeg("shannon3.jpg", height = 4, width = 5, units = 'in', res = 600)
ggplot(aes(y = Shannon, x = Group2, fill = Group1), data = gorilla_shannon) + geom_boxplot() + ggtitle("Gene Diversity") + labs(x="",y="Shannon") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold")) + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1))
dev.off ()




#------ Only Human
alpha_div_mean <- read.csv (file="Alpha_diversity_mean.Out", sep = "\t", header = TRUE, row.names = 1)
alpha_div_gorilla <- subset(alpha_div_mean, Group == "Gorilla")
gorilla_shannon <- alpha_div_gorilla[, c(1,7)]
Colors <- c("green4", "limegreen")
jpeg("shannon_gorilla_only_Final.jpg", height = 3, width = 2, units = 'in', res = 600)
ggplot(aes(y = Shannon, x = Group1, fill = Group1), data = gorilla_shannon) + geom_boxplot() + ggtitle("Gene Diversity") + labs(x="",y="Shannon") + theme_classic() + scale_color_manual(values=Colors) + scale_fill_manual(values=Colors) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(axis.text.x = element_text(size = 8, colour = "black", face = "bold"), axis.text.y = element_text(size = 8, colour = "black", face = "bold"))
dev.off ()
