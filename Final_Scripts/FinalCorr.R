#https://cran.r-project.org/web/packages/minerva/minerva.pdf

Sel_taxa <- read.csv(file="Total_Sel_Taxa.txt", sep="\t", row.names=1, header = T)
Sel_fun <- read.csv(file="Total_Sel_Functions.txt", sep="\t", row.names=1, header = T)
library (minerva)

Sel_taxa <- as.matrix(Sel_taxa)
Sel_fun <- as.matrix(Sel_fun)

cstats(Sel_taxa, Sel_fun, alpha = 0.6, C = 15, est = "mic_approx")

mictools(Sel_taxa, alpha = 9, C = 5, seed = 0, nperm = 2e+05, p.adjust.method = "BH")


Correlations <- mine(Sel_taxa, y=Sel_fun[,1], master = NULL, alpha = 0.6, C = 15, n.cores = 1, var.thr = 1e-05, eps = NULL, est = "mic_approx",na.rm = FALSE, use = "all.obs", normalization = FALSE)

#--- Finally to calculate MIC Maximal Information Coefficient
mine_stat(Sel_taxa[,1], Sel_fun[,1], alpha = 0.6, C = 15, est = "mic_approx", measure = "mic", eps = NA_real_, p = -1, norm = FALSE)




#------- CCREPE
# Using ccrepe to correlate Metab_F and Microbiota
ccrepe_taxa_fun <- ccrepe(x = Sel_taxa, y = Sel_fun, sim.score = nc.score, iterations = 50, min.subj = 5)

# Extract the correlation coefficients and q values from ccrepe_Microb_Metab_F
r_ccrepe_taxa_fun <- ccrepe_taxa_fun$sim.score
q_ccrepe_taxa_fun <- ccrepe_taxa_fun$q.values

r_ccrepe_taxa_fun <- melt(r_ccrepe_taxa_fun)
q_ccrepe_taxa_fun <- melt (q_ccrepe_taxa_fun)

q_r_ccrepe_taxa_fun <- cbind(r_ccrepe_taxa_fun,q_ccrepe_taxa_fun)
q_r_ccrepe_taxa_fun_Final <- q_r_ccrepe_taxa_fun[,c(1:3,6)]


colnames(q_r_ccrepe_taxa_fun_Final) <- c("Microbe", "Function", "r_value", "q_value")
q_r_ccrepe_taxa_fun_Final_formatted <- q_r_ccrepe_taxa_fun_Final %>% filter(!is.na(q_value))
write.table (q_r_ccrepe_taxa_fun_Final_formatted, file="Microbe_Fun_Corr.txt", sep = "\t")

#-- To make a plot
graph_taxa_path_Spear_corsTT <- q_r_ccrepe_taxa_fun_Final_formatted %>%
  filter(abs(r_value) > .4) %>%
  filter(abs(q_value) < 0.05) %>%
  graph_from_data_frame(directed = FALSE)

#--- This figure is replotted using Cytoscape

jpeg("Microbe_Fun_Corr_Spear.jpg", height = 12, width = 12, units = 'in', res = 600)
ggraph(graph_taxa_path_Spear_corsTT) +
  geom_edge_link(aes(edge_alpha = abs(r_value), edge_width = abs(r_value), color = r_value)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("darkred", "forestgreen")) +
  scale_edge_width(range = c(0.4, 1)) +
  geom_node_point(color = "darkgray", size = 3) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_graph() +
  labs(title = "Correlations between taxa and function")
dev.off ()

