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

