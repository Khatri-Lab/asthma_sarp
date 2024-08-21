###############################################################
library(readr)
library(data.table)
library(MetaIntegrator)
library(tidyverse)
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/plot_themes.R")
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/pooled_ROC_plot.R")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_expr.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/discovery_validation_bayes_meta_objs.rds")
###############################################################


###############################################################
compute_sample_scores <- function(meta_obj, split_name)
{
  pheno_score_frame <- Reduce(rbind, lapply(meta_obj$originalData, function(x) data.frame(bec_score = calculateScore(sig_genes, x, TRUE), dataset = x$formattedName, pheno = x$class)))
  pheno_score_frame$Split <- split_name
  return(pheno_score_frame)
}
combined_bec_scores <- rbind(compute_sample_scores(discovery_meta_obj, "Discovery"), compute_sample_scores(validation_meta_obj, "Validation"))
combined_bec_scores$pheno <- ifelse(combined_bec_scores$pheno == 0, "Healthy", "Asthma")
combined_bec_scores$pheno <- ordered(combined_bec_scores$pheno, levels = c("Healthy", "Asthma"))
score_plot <- ggplot(combined_bec_scores, aes(x = pheno, y = bec_score)) +
  theme(axis.title.x = element_blank()) +
  ylab("BEC asthma (z score)") +
  geom_beeswarm(cex = cex_val * 1.2, size = point_size, colour = black_text_colour) +
  geom_boxplot(fill = NA, lwd = line_size_main_mm, width = 0.4, outlier.shape = NA, colour = black_text_colour) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  stat_compare_means(method = "wilcox.test",
    comparisons = list(c("Healthy", "Asthma")),
    size = text_size_labels, colour = black_text_colour,
    bracket.size = line_size_text_repel_mm) +
  facet_wrap(vars(dataset), nrow = 2, scales = "free")
###############################################################


###############################################################
figures_folder <- "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/"

cairo_pdf(paste0(figures_folder, "fig_S2_BEC_score_by_dataset.pdf"), width = 7.2, height = 4, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))

print(score_plot, vp = viewport(layout.pos.row = 2:99, layout.pos.col = 2:99))

grid.text(x = unit(0.05, "npc"), y = unit(0.97, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.24, "npc"), y = unit(0.97, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.43, "npc"), y = unit(0.97, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.62, "npc"), y = unit(0.97, "npc"), label = "D", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.81, "npc"), y = unit(0.97, "npc"), label = "E", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.05, "npc"), y = unit(0.47, "npc"), label = "F", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.24, "npc"), y = unit(0.47, "npc"), label = "G", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.43, "npc"), y = unit(0.47, "npc"), label = "H", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.62, "npc"), y = unit(0.47, "npc"), label = "I", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.81, "npc"), y = unit(0.47, "npc"), label = "J", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))

dev.off()
###############################################################
