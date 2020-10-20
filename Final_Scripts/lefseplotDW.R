Pathways <- read.csv (file="DW_GSA_Stats_pathways_Dry_Wet.txt", sep = "\t", header = T)
Pathways_stats <- Pathways[,c(1,2,4)]

library("dplyr")
library("ggplot2")
df <- Pathways_stats %>%
  mutate(
    Pathways = factor(Pathways, levels = Pathways[order(Stat_dist.dir.up, decreasing = TRUE)]),
    label_y = ifelse(Stat_dist.dir.up < 0, 0.02, -0.02),
    label_hjust = ifelse(Stat_dist.dir.up < 0, 0, 1)
  )

my_plot <- ggplot(df, aes(x = Pathways, y = Stat_dist.dir.up, fill = Group)) +
  geom_bar(stat = "identity", col = "black") +
  geom_text(aes(y = label_y, label = Pathways, hjust = label_hjust)) +
  coord_flip() +
  scale_fill_manual(values = c(Dry = "wheat4", Wet = "wheat")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top",
        legend.justification = 0.9,
        plot.title = element_text(size = 24, face = "bold"),
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(expression("Stat_dist.dir.up"),
                     breaks = -6:6, limits = c(-9.5, 6))

print(my_plot)
ggsave("Pathways_Dry_vs_Wet.png", width = 10, height = 8, dpi = 600, my_plot)