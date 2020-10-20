Xeno <- read.csv (file="Dry_Wet_Lyfse_Plot_WithNames.txt", sep = "\t", header = T)
Xeno_stats <- Xeno[,c(1,2,3)]

library("dplyr")
library("ggplot2")
df <- Xeno_stats %>%
  mutate(
    Xeno = factor(Xeno, levels = Xeno[order(Log_Odds_Ratio, decreasing = TRUE)]),
    label_y = ifelse(Log_Odds_Ratio < 0, 0.02, -0.02),
    label_hjust = ifelse(Log_Odds_Ratio < 0, 0, 1)
)

my_plot <- ggplot(df, aes(x = Xeno, y = Log_Odds_Ratio, fill = Group)) +
  geom_bar(stat = "identity", col = "black") +
  geom_text(aes(y = label_y, label = Xeno, hjust = label_hjust)) +
  coord_flip() +
  scale_fill_manual(values = c(Dry = "wheat4", Wet = "wheat")) +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top",
        legend.justification = 0.9,
        plot.title = element_text(size = 20, face = "bold"),
        legend.title = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey80", linetype = "dashed"),
        panel.grid.minor.x = element_blank()) +
  scale_y_continuous(expression("Log_Odds_Ratio"),
                     breaks = -4.7:3.3)
print(my_plot)
ggsave("Xeno_Dry_Wet.Names.png", width = 12, height = 7, dpi = 600, my_plot)
