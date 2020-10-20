taxa_Selected <- read.csv (file = "Sel_taxa_proportions_Group.txt", sep = "\t", header = T, row.names = 1)
Samples <- row.names(taxa_Selected)
Group <- taxa_Selected$Group

#----To plot selected Microbes Total 16 (14 from first method and 2 additional (indval>0.7) from second method)
taxa_Sel_per <- sqrt(taxa_Selected[,2:17])*100
# Try a cube root
#taxa_Sel_per <- 100*taxa_Selected[,2:17]^(1/3)
taxa_Sel_Group <- cbind(Samples, Group, taxa_Sel_per)
pcm = melt(taxa_Sel_Group, id = c("Samples", "Group"))
pcm$Sample <- factor(pcm$Sample,levels=unique(pcm$Sample))

jpeg("Bubble_plot_Taxa.png", height = 6, width = 6, units = 'in', res = 600)
ggplot(pcm, aes(x = Sample, y = variable)) + 
  geom_point(aes(size = value, fill = Group), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Sqrt (Relative Abundance)%", fill = "Group")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 4, angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", size = 7, face = "bold"), 
        legend.text = element_text(size = 7, colour ="black"), 
        legend.title = element_text(size = 7, face = "bold"), panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 1.1), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +  
  scale_fill_manual(values = c("skyblue4", "skyblue", "wheat4", "wheat"), guide = guide_legend(override.aes = list(size=5))) + 
  scale_y_discrete(limits = rev(levels(pcm$variable))) 
dev.off ()
