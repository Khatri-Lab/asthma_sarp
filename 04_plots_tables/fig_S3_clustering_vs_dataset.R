###############################################################
library(tidyverse)
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/plot_themes.R")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/raw_data_umap_obj.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conorm_data_umap_obj.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_pheno_cluster.rds")
###############################################################


###############################################################
plot_conorm_data <- cbind(conorm_umap, combined_pheno)
names(cluster_colors) <- c("1", "2", "3", "4")
plot_conorm <- ggplot(plot_conorm_data, aes(x = V1, y = V2)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  geom_point(aes(colour = Dataset), size = point_size) +
  theme(legend.position = "none") +
  scale_colour_manual(values = dataset_colours)
plot_conorm_marginal <- ggMarginal(plot_conorm, groupColour = TRUE, groupFill = TRUE)
###############################################################


###############################################################
raw_umap$dataset[raw_umap$dataset == "SARP UCSF"] <- "SARP 3"
plot_raw <- ggplot(raw_umap, aes(x = V1, y = V2)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  geom_point(aes(colour = dataset), size = point_size) +
  theme(legend.position = "none") +
  scale_colour_manual(values = dataset_colours)
plot_raw_marginal <- ggMarginal(plot_raw, groupColour = TRUE, groupFill = TRUE)
###############################################################


###############################################################
cluster_summary <- plot_conorm_data %>% group_by(Cluster, Dataset) %>% summarize(Count = n()) %>% mutate(Prop = Count / sum(Count))
dataset_hist <- ggplot(cluster_summary, aes(x = as.factor(Cluster), y = Prop)) +
  geom_bar(aes(fill = Dataset), stat = "identity") +
  scale_fill_manual(values = dataset_colours) +
  ylab("Proportion") + xlab("Cluster") + theme(legend.position = "left") +
  labs(fill = "Dataset") + guides(fill = guide_legend(ncol = 1, override.aes = list(size = text_size_labels * 6 / 9)))
###############################################################


###############################################################
figures_folder <- "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/"
umap_combined <- plot_grid(plot_conorm_marginal, plot_raw_marginal, nrow = 1, align = "hv", axis = "tblr")

cairo_pdf(paste0(figures_folder, "fig_S3_UMAP_dataset.pdf"), width = 7.2, height = 2.5, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))

print(umap_combined, vp = viewport(layout.pos.row = 2:100, layout.pos.col = 1:65))
print(dataset_hist, vp = viewport(layout.pos.row = 2:100, layout.pos.col = 66:99))

grid.text(x = unit(0.02, "npc"), y = unit(0.965, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.32, "npc"), y = unit(0.965, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.80, "npc"), y = unit(0.965, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))

dev.off()
###############################################################

