#--- Pathway plots New
#--New1
Pathways_new <- read.csv (file = "../../../../KEGG_Analysis/Pathway_plots/Check_pathways_new2.txt", sep = "\t", row.names = 1, header =T)
Pathways_new <- data.frame(t(Pathways_new))
path2 <- sqrt(Pathways_new)
path3 <- path2*100

metadata_human <- read.csv (file = "../../../../Metadata_Human.txt", row.names=1, header = T, sep = "\t")
Groups <- metadata_human[,2:3]
pathway_human <- path3[24:51,]
Color <- c("skyblue4", "skyblue")

Pathways_sqrt_per_new <- cbind(pathway_human, Groups)

Vitamin <- Pathways_sqrt_per_new[,c(1,38)]
Vitamin_melt <- melt(Vitamin, id.vars = "Group1")
jpeg("H_Vitamin B6 metabolism.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Vitamin_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Vitamin B6 metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold")) 
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Vitamin.B6.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 156, p-value = 0.006743

Methane <- Pathways_sqrt_per_new[,c(2,38)]
Methane_melt <- melt(Methane, id.vars = "Group1")
jpeg("H_Methane metabolism.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Methane_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Methane metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold")) 
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Methane.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 23, p-value = 0.0002687

Inositol <- Pathways_sqrt_per_new[,c(3,38)]
Inositol_melt <- melt(Inositol, id.vars = "Group1")
jpeg("H_InositolPhosphate metabolism.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Inositol_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Inositol phosphate metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold")) 
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Inositol.phosphate.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 56, p-value = 0.05562

MicrobialMetabolism <- Pathways_sqrt_per_new[,c(4,38)]
MicrobialMetabolism_melt <- melt(MicrobialMetabolism, id.vars = "Group1")
jpeg("H_MicrobialMetabolism.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = MicrobialMetabolism_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Microbial metabolism in diverse environments") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Microbial.metabolism.in.diverse.environments ~ Pathways_sqrt_per_new$Group1)
#W = 27, p-value = 0.0006462

Benzoatedegradation <- Pathways_sqrt_per_new[,c(5,38)]
Benzoatedegradation_melt <- melt(Benzoatedegradation, id.vars = "Group1")
jpeg("H_Benzoatedegradation.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Benzoatedegradation_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Benzoate degradation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Benzoate.degradation ~ Pathways_sqrt_per_new$Group1)
#W = 17, p-value = 5.973e-05

AminoBenzoatedegradation <- Pathways_sqrt_per_new[,c(6,38)]
AminoBenzoatedegradation_melt <- melt(AminoBenzoatedegradation, id.vars = "Group1")
jpeg("H_AminoBenzoatedegradation.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = AminoBenzoatedegradation_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Aminobenzoate degradation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Aminobenzoate.degradation ~ Pathways_sqrt_per_new$Group1)
#W = 51, p-value = 0.03103

GalactoseMetabolism <- Pathways_sqrt_per_new[,c(7,38)]
GalactoseMetabolism_melt <- melt(GalactoseMetabolism, id.vars = "Group1")
jpeg("H_GalactoseMetabolism.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = GalactoseMetabolism_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Galactose metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Galactose.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 74, p-value = 0.2852

ValineLeucineIsoLeucineDegradation <- Pathways_sqrt_per_new[,c(8,38)]
ValineLeucineIsoLeucineDegradation_melt <- melt(ValineLeucineIsoLeucineDegradation, id.vars = "Group1")
jpeg("H_ValineLeucineIsoLeucineDegradation.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = ValineLeucineIsoLeucineDegradation_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Valine leucine and isoleucine.degradation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Valine..leucine.and.isoleucine.degradation ~ Pathways_sqrt_per_new$Group1)
#W = 0, p-value = 4.985e-08

Ascorbate.and.aldarate.metabolism <- Pathways_sqrt_per_new[,c(9,38)]
Ascorbate.and.aldarate.metabolism_melt <- melt(Ascorbate.and.aldarate.metabolism, id.vars = "Group1")
jpeg("H_Ascorbate.and.aldarate.metabolism.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Ascorbate.and.aldarate.metabolism_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Ascorbate and aldarate metabolism_melt") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Ascorbate.and.aldarate.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 105, p-value = 0.7688

Fatty.acid.degradation <- Pathways_sqrt_per_new[,c(10,38)]
Fatty.acid.degradation_melt <- melt(Fatty.acid.degradation, id.vars = "Group1")
jpeg("H_Fatty.acid.degradation.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Fatty.acid.degradation_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Fatty acid degradation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Fatty.acid.degradation ~ Pathways_sqrt_per_new$Group1)
#W = 83, p-value = 0.5112

Lysine.degradation <- Pathways_sqrt_per_new[,c(11,38)]
Lysine.degradation_melt <- melt(Lysine.degradation, id.vars = "Group1")
jpeg("H_Lysine.degradation.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Lysine.degradation_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Lysine.degradation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Lysine.degradation ~ Pathways_sqrt_per_new$Group1)
#W = 26, p-value = 0.0005233

Tyrosinemetabolism <- Pathways_sqrt_per_new[,c(12,38)]
Tyrosinemetabolism_melt <- melt(Tyrosinemetabolism, id.vars = "Group1")
jpeg("H_Tyrosinemetabolism.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Tyrosinemetabolism_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Tyrosine metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Tyrosine.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 62, p-value = 0.1035

beta.Alanine.metabolism <- Pathways_sqrt_per_new[,c(13,38)]
beta.Alanine.metabolism_melt <- melt(beta.Alanine.metabolism, id.vars = "Group1")
jpeg("H_beta.Alanine.metabolism.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = beta.Alanine.metabolism_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("beta-Alanine metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$beta.Alanine.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 110, p-value = 0.6027

Styrene.degradation <- Pathways_sqrt_per_new[,c(14,38)]
Styrene.degradation_melt <- melt(Styrene.degradation, id.vars = "Group1")
jpeg("H_Styrene.degradation.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Styrene.degradation_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Styrene.degradation_melt") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Styrene.degradation ~ Pathways_sqrt_per_new$Group1)
#W = 34, p-value = 0.002486

ValineLeucineIsoLeucineBiosynthesis <- Pathways_sqrt_per_new[,c(15,38)]
ValineLeucineIsoLeucineBiosynthesis_melt <- melt(ValineLeucineIsoLeucineBiosynthesis, id.vars = "Group1")
jpeg("H_ValineLeucineIsoLeucineBiosynthesis.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = ValineLeucineIsoLeucineBiosynthesis_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Valine leucine and isoleucine biosynthesis") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Valine..leucine.and.isoleucine.biosynthesis ~ Pathways_sqrt_per_new$Group1)
#W = 95, p-value = 0.91

Nicotinate.and.nicotinamide.metabolism <- Pathways_sqrt_per_new[,c(16,38)]
Nicotinate.and.nicotinamide.metabolism_melt <- melt(Nicotinate.and.nicotinamide.metabolism, id.vars = "Group1")
jpeg("H_Nicotinate.and.nicotinamide.metabolism.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Nicotinate.and.nicotinamide.metabolism_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Nicotinate and nicotinamide metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Nicotinate.and.nicotinamide.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 118, p-value = 0.3761

Phenylalanine..tyrosine.and.tryptophan.biosynthesis <- Pathways_sqrt_per_new[,c(17,38)]
Phenylalanine..tyrosine.and.tryptophan.biosynthesis_melt <- melt(Phenylalanine..tyrosine.and.tryptophan.biosynthesis, id.vars = "Group1")
jpeg("H_Phenylalanine..tyrosine.and.tryptophan.biosynthesis.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Phenylalanine..tyrosine.and.tryptophan.biosynthesis_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Phenylalanine..tyrosine.and.tryptophan.biosynthesis") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Phenylalanine..tyrosine.and.tryptophan.biosynthesis ~ Pathways_sqrt_per_new$Group1)
#W = 143, p-value = 0.03948

Chlorocyclohexane.and.chlorobenzene.degradation <- Pathways_sqrt_per_new[,c(18,38)]
Chlorocyclohexane.and.chlorobenzene.degradation_melt <- melt(Chlorocyclohexane.and.chlorobenzene.degradation, id.vars = "Group1")
jpeg("H_Chlorocyclohexane.and.chlorobenzene.degradation.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Chlorocyclohexane.and.chlorobenzene.degradation_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Chlorocyclohexane.and.chlorobenzene.degradation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Chlorocyclohexane.and.chlorobenzene.degradation ~ Pathways_sqrt_per_new$Group1)
#W = 71, p-value = 0.2273

Alanine..aspartate.and.glutamate.metabolism <- Pathways_sqrt_per_new[,c(19,38)]
Alanine..aspartate.and.glutamate.metabolism_melt <- melt(Alanine..aspartate.and.glutamate.metabolism, id.vars = "Group1")
jpeg("H_Alanine..aspartate.and.glutamate.metabolism.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Alanine..aspartate.and.glutamate.metabolism_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Alanine..aspartate.and.glutamate.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Alanine..aspartate.and.glutamate.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 84, p-value = 0.5409

Taurine.and.hypotaurine.metabolism <- Pathways_sqrt_per_new[,c(20,38)]
Taurine.and.hypotaurine.metabolism_melt <- melt(Taurine.and.hypotaurine.metabolism, id.vars = "Group1")
jpeg("H_Taurine.and.hypotaurine.metabolism.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Taurine.and.hypotaurine.metabolism_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Taurine.and.hypotaurine.metabolism_melt") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Taurine.and.hypotaurine.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 131, p-value = 0.1371

Tropane..piperidine.and.pyridine.alkaloid.biosynthesis <- Pathways_sqrt_per_new[,c(22,38)]
Tropane..piperidine.and.pyridine.alkaloid.biosynthesis_melt <- melt(Tropane..piperidine.and.pyridine.alkaloid.biosynthesis, id.vars = "Group1")
jpeg("H_Tropane..piperidine.and.pyridine.alkaloid.biosynthesis.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Tropane..piperidine.and.pyridine.alkaloid.biosynthesis_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Tropane..piperidine.and.pyridine.alkaloid.biosynthesis") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Tropane..piperidine.and.pyridine.alkaloid.biosynthesis ~ Pathways_sqrt_per_new$Group1)
#W = 78, p-value = 0.3761

Nitrotoluene.degradation <- Pathways_sqrt_per_new[,c(25,38)]
Nitrotoluene.degradation_melt <- melt(Nitrotoluene.degradation, id.vars = "Group1")
jpeg("H_Nitrotoluene.degradation.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Nitrotoluene.degradation_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Nitrotoluene.degradation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Nitrotoluene.degradation ~ Pathways_sqrt_per_new$Group1)
#W = 42, p-value = 0.009146

Glycerolipid.metabolism <- Pathways_sqrt_per_new[,c(26,38)]
Glycerolipid.metabolism_melt <- melt(Glycerolipid.metabolism, id.vars = "Group1")
jpeg("H_Glycerolipid.metabolism.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Glycerolipid.metabolism_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Glycerolipid.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Glycerolipid.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 48, p-value = 0.02119

Pyrimidine.metabolism <- Pathways_sqrt_per_new[,c(27,38)]
Pyrimidine.metabolism_melt <- melt(Pyrimidine.metabolism, id.vars = "Group1")
jpeg("H_Pyrimidine.metabolism.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Pyrimidine.metabolism_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Pyrimidine.metabolism") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Pyrimidine.metabolism ~ Pathways_sqrt_per_new$Group1)
#W = 24, p-value = 0.0003375

Biosynthesis.of.amino.acids <- Pathways_sqrt_per_new[,c(28,38)]
Biosynthesis.of.amino.acids_melt <- melt(Biosynthesis.of.amino.acids, id.vars = "Group1")
jpeg("H_Biosynthesis.of.amino.acids.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Biosynthesis.of.amino.acids_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Biosynthesis.of.amino.acids") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Biosynthesis.of.amino.acids ~ Pathways_sqrt_per_new$Group1)
#W = 81, p-value = 0.4544

Degradation.of.aromatic.compounds <- Pathways_sqrt_per_new[,c(29,38)]
Degradation.of.aromatic.compounds_melt <- melt(Degradation.of.aromatic.compounds, id.vars = "Group1")
jpeg("H_Degradation.of.aromatic.compounds_melt.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Degradation.of.aromatic.compounds_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Degradation.of.aromatic.compounds") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Degradation.of.aromatic.compounds ~ Pathways_sqrt_per_new$Group1)
#W = 48, p-value = 0.02119

Oxidative.phosphorylation <- Pathways_sqrt_per_new[,c(30,38)]
Oxidative.phosphorylation_melt <- melt(Oxidative.phosphorylation, id.vars = "Group1")
jpeg("H_Oxidative.phosphorylation.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Oxidative.phosphorylation_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Oxidative.phosphorylation") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Oxidative.phosphorylation ~ Pathways_sqrt_per_new$Group1)
#W = 54, p-value = 0.04437

Biosynthesis.of.secondary.metabolites <- Pathways_sqrt_per_new[,c(31,38)]
Biosynthesis.of.secondary.metabolites_melt <- melt(Biosynthesis.of.secondary.metabolites, id.vars = "Group1")
jpeg("H_Biosynthesis.of.secondary.metabolites.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Biosynthesis.of.secondary.metabolites_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Biosynthesis.of.secondary.metabolites") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Biosynthesis.of.secondary.metabolites ~ Pathways_sqrt_per_new$Group1)
#W = 97, p-value = 0.982

Lipopolysaccharide.biosynthesis <- Pathways_sqrt_per_new[,c(32,38)]
Lipopolysaccharide.biosynthesis_melt <- melt(Lipopolysaccharide.biosynthesis, id.vars = "Group1")
jpeg("H_Lipopolysaccharide.biosynthesis.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = Lipopolysaccharide.biosynthesis_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Lipopolysaccharide.biosynthesis") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Lipopolysaccharide.biosynthesis ~ Pathways_sqrt_per_new$Group1)
#W = 164, p-value = 0.001729

One.carbon.pool.by.folate <- Pathways_sqrt_per_new[,c(33,38)]
One.carbon.pool.by.folate_melt <- melt(One.carbon.pool.by.folate, id.vars = "Group1")
jpeg("H_One.carbon.pool.by.folate.jpg", height = 3, width = 1, units = 'in', res = 600)
ggplot(data = One.carbon.pool.by.folate_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("One.carbon.pool.by.folate") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$One.carbon.pool.by.folate ~ Pathways_sqrt_per_new$Group1)
#W = 37, p-value = 0.004162

Streptomycin.biosynthesis <- Pathways_sqrt_per_new[,c(34,38)]
Streptomycin.biosynthesis_melt <- melt(Streptomycin.biosynthesis, id.vars = "Group1")
jpeg("H_Streptomycin.biosynthesis_new.jpg", height = 2.5, width = 1, units = 'in', res = 600)
#ggplot(data = Streptomycin.biosynthesis_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Streptomycin.biosynthesis") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
ggplot(data = Streptomycin.biosynthesis, aes(x=Group1, y=Streptomycin.biosynthesis, fill=Group1)) + geom_boxplot() + ggtitle("Streptomycin.biosynthesis") + labs(x="",y="sqrt (Relative abundance)%") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Streptomycin.biosynthesis ~ Pathways_sqrt_per_new$Group1)
#W = 179, p-value = 5.973e-05

Two.component.system <- Pathways_sqrt_per_new[,c(35,38)]
Two.component.system_melt <- melt(Two.component.system, id.vars = "Group1")
jpeg("H_Two.component.system.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Two.component.system_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Two.component.system") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Two.component.system ~ Pathways_sqrt_per_new$Group1)
#W = 40, p-value = 0.006743

ABC.transporters <- Pathways_sqrt_per_new[,c(36,38)]
ABC.transporters_melt <- melt(ABC.transporters, id.vars = "Group1")
jpeg("H_ABC.transporters.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = ABC.transporters_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("ABC.transporters") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$ABC.transporters ~ Pathways_sqrt_per_new$Group1)
#W = 12, p-value = 1.356e-05

Phosphotransferase.system..PTS <- Pathways_sqrt_per_new[,c(37,38)]
Phosphotransferase.system..PTS_melt <- melt(Phosphotransferase.system..PTS, id.vars = "Group1")
jpeg("H_Phosphotransferase.system..PTS.jpg", height = 2.5, width = 1, units = 'in', res = 600)
ggplot(data = Phosphotransferase.system..PTS_melt, aes(x=reorder(Group1, value, FUN=median), y=value, fill=Group1)) + geom_boxplot() + ggtitle("Phosphotransferase.system..PTS") + labs(x="",y="sqrt(Relative abundance) %") + theme_classic() + scale_color_manual(values=Color) + scale_fill_manual(values=Color) + theme(plot.title = element_text(size = 6, face = "bold")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.text = element_text(colour="black", size=8)) + theme(legend.title = element_blank()) + theme(legend.position='none') + theme(axis.text.x = element_text(size = 10, colour = "black", face = "bold"), axis.text.y = element_text(size = 10, colour = "black", face = "bold"))
dev.off ()
wilcox.test(Pathways_sqrt_per_new$Phosphotransferase.system..PTS. ~ Pathways_sqrt_per_new$Group1)
#W = 45, p-value = 0.01411

