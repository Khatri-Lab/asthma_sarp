###############################################################
library(readr)
library(tidyverse)
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/plot_themes.R")
###############################################################


###############################################################
tile_height <- unit(base_text_size, "pt")
tile_width <- unit(base_text_size, "pt")
ht_opt("simple_anno_size" = tile_height)
ht_opt("heatmap_row_names_gp" = text_setting_small_it)
figures_folder <- "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/"
###############################################################


###############################################################
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/GSE115770_baseline_expr_pheno.rds")
GSE115770_pheno_cluster <- read_csv("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/GSE115770_baseline_pheno_clusters.csv")
sarp_pheno_expr <- read_csv("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/SARP_cluster_pheno_expr.csv")
gene_list <- c("POSTN", "SERPINB2", "CLCA1", "CCL26", "CDH26", "ITLN1", "MS4A2", "CPA3", "ALOX15", "IL18R1", "GATA2", "CXCL1", "SCGB3A1", "SCGB1A1", "MUC5B", "MUC12", "CEACAM5", "NOTCH2", "RBPJ", "HEY1", "FKBP5", "HSD11B2", "CST1", "TOP2A")
###############################################################


###############################################################
baseline_expr <- as.data.frame(baseline_expr[GSE115770_pheno_cluster$library.sampleId...1, ])
baseline_expr$Cluster <- GSE115770_pheno_cluster$pred_cluster_raw_sarp
###############################################################


###############################################################
cluster_summary <- GSE115770_pheno_cluster %>% group_by(pred_cluster_raw_sarp) %>% summarize(count = n())
cluster_plot <- ggplot(cluster_summary, aes(y = as.factor(pred_cluster_raw_sarp), x = count, fill = as.factor(pred_cluster_raw_sarp))) +
  xlab("# Subjects (baseline)") + ylab("Cluster") +
  geom_bar(stat = "identity", width = 0.8) +
  theme(legend.position = "none") +
  scale_fill_manual(values = cluster_colors) +
  scale_x_continuous(breaks = c(6, 28, 55))
###############################################################


###############################################################
min_expr <- -1.1547
max_expr <- 1.1547
expr_matrix_colour_pts <- c(min_expr, min_expr/2, 0, max_expr/2, max_expr)
#expr_cols <- colorRamp2(expr_matrix_colour_pts, rev(c("#762A83", "#C2A5CF", "#F7F7F7", "#ACD39E", "#1B7837")))
expr_cols <- colorRamp2(expr_matrix_colour_pts, scico(5, palette = "vik"))
expr_matrix_legend_pts <- c(min_expr, max_expr)
expr_matrix_legend_labels <- c(round(min_expr, 2), round(max_expr, 2))
expr_legend <- Legend(col_fun = expr_cols,
  at = expr_matrix_legend_pts,
  labels = expr_matrix_legend_labels,
  title = "Row-scaled\nexpression", title_position = "topcenter",
  grid_height = unit(base_text_size / 3, "pt"),
  legend_width = 4 * unit(base_text_size, "pt"),
  background = "transparent",
  title_gp = text_setting_small,
  title_gap = gap_size, gap = gap_size,
  labels_gp = text_setting_sub_small,
  direction = "horizontal")
###############################################################


###############################################################
plot_mean_heatmap <- function(data_matrix, plot_title)
{
  data_matrix[, gene_list] <- log2(1 + data_matrix[, gene_list])
  data_matrix <- subset(data_matrix, select = c(gene_list, "Cluster"))
  data_matrix <- data_matrix[data_matrix$Cluster != 2, ]
  heatmap_matrix <- as.data.frame(data_matrix %>% group_by(Cluster) %>% summarise_all(mean))
  rownames(heatmap_matrix) <- paste0("Cluster ", heatmap_matrix$Cluster)
  heatmap_matrix <- subset(heatmap_matrix, select = -c(Cluster))
  heatmap_matrix <- t(scale(heatmap_matrix))
  ht_mp <- Heatmap(heatmap_matrix,
      col = expr_cols,
      width = ncol(heatmap_matrix) * tile_width,
      name = plot_title,
      show_heatmap_legend = FALSE,
      row_dend_gp = dend_lines_gp,
      row_dend_width = unit(2 * base_text_size, "pt"),
      cluster_columns = FALSE,
      gap = gap_size,
      column_title = plot_title)
  return(ht_mp)
}
altmann_htmp <- plot_mean_heatmap(baseline_expr, "Altman")
sarp_12_htmp <- plot_mean_heatmap(sarp_pheno_expr[sarp_pheno_expr$dataset == "sarp_1_2", ], "SARP 1-2")
sarp_3_htmp <- plot_mean_heatmap(sarp_pheno_expr[sarp_pheno_expr$dataset == "sarp_3", ], "SARP 3")
mean_expr_ht_mp <- altmann_htmp + sarp_12_htmp + sarp_3_htmp
###############################################################


###############################################################
GSE115770_pheno_cluster$cluster_factor <- as.factor(GSE115770_pheno_cluster$pred_cluster_raw_sarp)
clinical_var_summary <- GSE115770_pheno_cluster[GSE115770_pheno_cluster$pred_cluster_raw_sarp != 1, c("FEV1.Percent.Predicted.at.Visit", "FEV1..FVC.at.Visit", "Exacerbations.total", "FeNO.at.Visit", "Blood.Eosinophil.Count", "Nasal.Eosinophil.Count", "cluster_factor")]
colnames(clinical_var_summary) <- c("FEV1 (%pred)", "FEV1/FVC ratio", "Exacerbations", "FeNO (ppb)", "Blood eosinophil count", "Nasal eosinophil count", "cluster")
clinical_var_summary <- clinical_var_summary %>% pivot_longer(!cluster, names_to = "clinical_var", values_to = "clinical_values")
clinical_var_summary$clinical_var <- ordered(clinical_var_summary$clinical_var, levels = c("FEV1 (%pred)", "FEV1/FVC ratio", "Exacerbations", "FeNO (ppb)", "Blood eosinophil count", "Nasal eosinophil count"))
clinical_var_summary_plot <- ggplot(clinical_var_summary, aes(x = cluster, y = clinical_values, fill = cluster)) +
  xlab("Cluster") +
  geom_boxplot(lwd = line_size_main_mm, width = 0.4, colour = black_text_colour, outlier.shape = NA) +
  geom_beeswarm(cex = cex_val, size = point_size, colour = black_text_colour) +
  stat_compare_means(method = "wilcox.test",
    comparisons = list(c("3", "4")),
    size = text_size_labels, colour = black_text_colour,
    bracket.size = line_size_text_repel_mm) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  scale_fill_manual(values = cluster_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.075, 0.125))) +
  facet_wrap(vars(clinical_var), nrow = 3, scales = "free_y")
###############################################################


###############################################################
cairo_pdf(paste0(figures_folder, "fig_6_nasal_lavage.pdf"), width = 7.2, height = 5.2, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))

print(cluster_plot, vp = viewport(layout.pos.row = 1:18, layout.pos.col = 4:43))

pushViewport(viewport(layout.pos.row = 19:98, layout.pos.col = 2:44))
draw(mean_expr_ht_mp, height = length(gene_list) * tile_height, gap = 3 * gap_size, background = "transparent", newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 90:99, layout.pos.col = 38:43))
draw(expr_legend)
upViewport(1)

print(clinical_var_summary_plot, vp = viewport(layout.pos.row = 2:99, layout.pos.col = 48:99))

grid.text(x = unit(0.02, "npc"), y = unit(0.98, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.02, "npc"), y = unit(0.81, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.49, "npc"), y = unit(0.97, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.74, "npc"), y = unit(0.97, "npc"), label = "D", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.49, "npc"), y = unit(0.655, "npc"), label = "E", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.74, "npc"), y = unit(0.655, "npc"), label = "F", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.49, "npc"), y = unit(0.35, "npc"), label = "G", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.74, "npc"), y = unit(0.35, "npc"), label = "H", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))

dev.off()
###############################################################

