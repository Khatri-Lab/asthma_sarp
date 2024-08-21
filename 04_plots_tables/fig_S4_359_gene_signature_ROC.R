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
sig_genes_subset <- sig_genes
sig_genes_subset$posGeneNames <- sig_genes_subset$posGeneNames[sig_genes_subset$posGeneNames %in% rownames(combined_expr)]
sig_genes_subset$negGeneNames <- sig_genes_subset$negGeneNames[sig_genes_subset$negGeneNames %in% rownames(combined_expr)]
ROC_discovery <- IL_pooledROCPlot(discovery_meta_obj, sig_genes_subset, custom_colors = discovery_data_colours, title = "Discovery, 359 gene subset") + theme(legend.background = element_rect(fill = "#FFFFFFEE", linewidth = 0, colour = "transparent"), panel.grid.major = element_line(linewidth = line_size_text_repel_mm))
ROC_validation <- IL_pooledROCPlot(validation_meta_obj, sig_genes_subset, custom_colors = validation_data_colours, title = "Validation, 359 gene subset") + theme(legend.background = element_rect(fill = "#FFFFFFEE", linewidth = 0, colour = "transparent"), panel.grid.major = element_line(linewidth = line_size_text_repel_mm))
ROC_combined <- plot_grid(ROC_discovery, ROC_validation, nrow = 1, align = "hv", axis = "tblr")
###############################################################


###############################################################
compute_sample_scores <- function(meta_obj, split_name)
{
  pheno_score_frame <- Reduce(rbind, lapply(meta_obj$originalData, function(x) data.frame(score_505 = calculateScore(sig_genes, x, TRUE), score_359 = calculateScore(sig_genes_subset, x, TRUE))))
  pheno_score_frame$Split <- split_name
  return(pheno_score_frame)
}
combined_bec_scores <- rbind(compute_sample_scores(discovery_meta_obj, "Discovery"), compute_sample_scores(validation_meta_obj, "Validation"))
score_corr_plot <- ggplot(combined_bec_scores, aes(x = score_505, y = score_359)) +
  theme(legend.position = "none", panel.grid.major = element_line(linewidth = line_size_text_repel_mm)) +
  xlab("505 gene signature score") + ylab("359 gene subset score") +
  geom_point(size = point_size, color = black_line_colour) +
  geom_smooth(method = "lm", formula = y~x, linewidth = line_size_main_mm, color = black_text_colour) +
  stat_cor(method = "pearson", label.x = -1, label.y = 2, size = text_size_labels, family = base_font_family, colour = black_text_colour) +
  facet_wrap(vars(Split), nrow = 1, scales = "free")
###############################################################


###############################################################
figures_folder <- "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/"

cairo_pdf(paste0(figures_folder, "fig_S4_359_gene_signature_ROC.pdf"), width = 7.2, height = 5.5, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))

print(ROC_combined, vp = viewport(layout.pos.row = 2:60, layout.pos.col = 2:99))
print(score_corr_plot, vp = viewport(layout.pos.row = 61:100, layout.pos.col = 2:99))

grid.text(x = unit(0.02, "npc"), y = unit(0.97, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.52, "npc"), y = unit(0.97, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.02, "npc"), y = unit(0.375, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.52, "npc"), y = unit(0.375, "npc"), label = "D", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))

dev.off()
###############################################################

