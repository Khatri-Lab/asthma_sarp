###############################################################
library(tidyverse)
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/plot_themes.R")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conorm_data_umap_obj.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_pheno_cluster.rds")
###############################################################


###############################################################
plot_cluster_data <- cbind(conorm_umap, combined_pheno)
names(cluster_colors) <- c("1", "2", "3", "4")
plot_cluster <- ggplot(plot_cluster_data, aes(x = V1, y = V2)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  geom_point(aes(colour = as.factor(Cluster)), size = point_size) +
  theme(legend.position = "none") +
  scale_colour_manual(values = cluster_colors)
plot_cluster_marginal <- ggMarginal(plot_cluster, groupColour = TRUE, groupFill = TRUE)
###############################################################


###############################################################
plot_severity <- ggplot(plot_cluster_data, aes(x = V1, y = V2)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  geom_point(aes(colour = Severity), size = point_size) +
  theme(legend.position = "none") +
  scale_colour_manual(values = severity_colours)
plot_severity_marginal <- ggMarginal(plot_severity, groupColour = TRUE, groupFill = TRUE)
###############################################################


###############################################################
cluster_summary <- plot_cluster_data %>% group_by(Cluster, Severity) %>% summarize(Count = n()) %>% mutate(Prop = Count / sum(Count))
cluster_hist <- ggplot(cluster_summary %>% group_by(Cluster) %>% summarize(Total = sum(Count)), aes(x = as.factor(Cluster), y = Total, fill = as.factor(Cluster))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  ylab("Count") + xlab("Cluster") +
  labs(fill = "Cluster") + guides(fill = guide_legend(ncol = 1, override.aes = list(size = text_size_labels * 6 / 9))) +
  theme(panel.grid.major.y = element_line(linewidth = line_size_text_repel_mm))
###############################################################


###############################################################
severity_hist <- ggplot(cluster_summary, aes(x = as.factor(Cluster), y = Prop)) +
  geom_bar(aes(fill = Severity), stat = "identity") +
  scale_fill_manual(values = severity_colours) +
  ylab("Proportion") + xlab("Cluster") +
  labs(fill = "Severity") + guides(fill = guide_legend(ncol = 1, override.aes = list(size = text_size_labels * 6 / 9)))
###############################################################


###############################################################
figures_folder <- "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/"
umap_combined <- plot_grid(plot_cluster_marginal, plot_severity_marginal, nrow = 1, align = "hv", axis = "tblr")
bars_combined <- plot_grid(cluster_hist, severity_hist, nrow = 1, align = "hv", axis = "tblr")

cairo_pdf(paste0(figures_folder, "fig_2_UMAP_clustering.pdf"), width = 7.2, height = 5, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))

print(umap_combined, vp = viewport(layout.pos.row = 2:60, layout.pos.col = 2:99))
print(bars_combined, vp = viewport(layout.pos.row = 61:100, layout.pos.col = 2:99))

grid.text(x = unit(0.02, "npc"), y = unit(0.98, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.515, "npc"), y = unit(0.98, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.02, "npc"), y = unit(0.385, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.515, "npc"), y = unit(0.385, "npc"), label = "D", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))

dev.off()
###############################################################

