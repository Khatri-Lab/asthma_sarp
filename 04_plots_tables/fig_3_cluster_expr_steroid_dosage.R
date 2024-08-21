###############################################################
library(tidyverse)
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/plot_themes.R")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_expr.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_pheno_cluster.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/SARP_meta_objs.rds")
###############################################################


###############################################################
tile_height <- unit(base_text_size, "pt")
tile_width <- unit(base_text_size, "pt")
ht_opt("simple_anno_size" = tile_height)
ht_opt("heatmap_row_names_gp" = text_setting_small_it)
figures_folder <- "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/"
names(cluster_colors) <- c(1:4)
###############################################################


###############################################################
SARP_3_cluster_pheno <- combined_pheno[combined_pheno$Dataset == "SARP 3", ]
genes_of_interest <- c("HEY1", "NPNT", "FHDC1", "CXCL1", "MUC5B", "C3", "SCGB3A1", "NOTCH2", "RBPJ", "SCGB1A1", "ATP13A4", "HEG1", "DAPK1", "TOP2A", "MUC12", "NDUFA8", "ZMAT2", "SYAP1", "GATA2", "CPA3", "MS4A2", "CLCA1", "ITLN1", "CCL26", "CST1", "SERPINB2", "ALOX15", "CDH26", "POSTN", "CEACAM5", "IL18R1", "KRT23", "HSD11B2", "FKBP5")
combined_data <- cbind(combined_pheno, t(combined_expr[genes_of_interest, ]))
combined_data <- combined_data %>% pivot_longer(genes_of_interest, names_to = "genes", values_to = "expr") %>%
  group_by(Cluster, Sex, Severity, genes) %>% summarize(mean_expr = mean(expr)) %>%
  pivot_wider(names_from = "genes", values_from = "mean_expr")
combined_data <- na.omit(combined_data[combined_data$Severity != "Unspecified", ])
combined_data$Cluster <- ordered(combined_data$Cluster, levels = c("1", "2", "3", "4"))
combined_data$Severity <- ordered(combined_data$Severity, levels = c("Healthy", "Mild/moderate", "Severe"))
combined_data$Sex <- ordered(combined_data$Sex, levels = c("female", "male"))
combined_data <- combined_data[with(combined_data, order(Cluster, Sex, Severity)), ]
combined_pheno <- combined_data[, c("Cluster", "Sex", "Severity")]
combined_expr_cont <- combined_data[, genes_of_interest]
combined_expr_cont <- as.matrix(t(scale(combined_expr_cont)))
###############################################################


###############################################################
combined_expr <- combined_expr_cont
combined_expr[combined_expr_cont > 2] <- "(2, 3.03]"
combined_expr[combined_expr_cont <= 2 & combined_expr_cont > 1] <- "(1, 2]"
combined_expr[combined_expr_cont <= 1 & combined_expr_cont > 0] <- "(0, 1]"
combined_expr[combined_expr_cont > -1 & combined_expr_cont <= 0] <- "(-1, 0]"
combined_expr[combined_expr_cont > -2 & combined_expr_cont <= -1] <- "(-2, -1]"
combined_expr[combined_expr_cont <= -2] <- "[-4.01, -2]"
###############################################################


###############################################################
expr_cols <- scico(6, palette = "berlin")
names(expr_cols) <- rev(c("(2, 3.03]", "(1, 2]", "(0, 1]", "(-1, 0]", "(-2, -1]", "[-4.01, -2]"))

create_discrete_legend <- function(legend_colours, legend_title, num_col)
{
  discrete_legend <- Legend(labels = names(legend_colours),
    legend_gp = gpar(fill = legend_colours, color = "transparent"),
    title = legend_title, title_position = "topleft",
    grid_height = unit(base_text_size, "pt"),
    grid_width = unit(base_text_size, "pt"),
    background = "transparent",
    title_gp = text_setting_small,
    title_gap = gap_size, gap = gap_size,
    labels_gp = text_setting_sub_small,
    direction = "horizontal", ncol = num_col, by_row = TRUE)
}

cluster_legend <- create_discrete_legend(cluster_colors, "Cluster", 4)
severity_legend <- create_discrete_legend(severity_colours[1:3], "Severity", 3)
sex_legend <- create_discrete_legend(sex_colors, "Sex", 2)
expr_legend <- create_discrete_legend(expr_cols, "Row-scaled expression", 3)

ht_mp_annotation <- HeatmapAnnotation(
  Cluster = combined_pheno$Cluster, Sex = combined_pheno$Sex, Severity = combined_pheno$Severity,
  col = list(Cluster = cluster_colors, Sex = sex_colors, Severity = severity_colours),
  show_legend = FALSE,
  annotation_name_gp = text_setting_small,
  annotation_name_side = "left")

ht_mp_expr <- Heatmap(combined_expr,
  col = expr_cols, na_col = "#FFEE99",
  show_column_names = FALSE,
  rect_gp = dend_lines_gp,
  width = ncol(combined_expr) * tile_width,
  name = "expr_matrix", column_title = "Mean expression in group",
  column_title_side = "bottom",
  row_names_side = "left",
  show_heatmap_legend = FALSE,
  column_split = combined_pheno$Cluster,
  gap = gap_size, top_annotation = ht_mp_annotation)
###############################################################


###############################################################
modified_subj_id <- SARP_3_cluster_pheno %>% rownames() %>% sub(pattern = "X", replacement = "") %>% sub(pattern = "\\.", replacement = "-") %>% sub(pattern = "\\.", replacement = "-")
SARP_3_cluster_pheno$participant <- as.character(modified_subj_id)
SARP_3_cluster_pheno <- merge(SARP_3_cluster_pheno, SARP_3$pheno, by = "participant")
SARP_3_cluster_pheno <- subset(SARP_3_cluster_pheno, select = c('rxhx_1120', 'rxhx_1232', 'Cluster'))
colnames(SARP_3_cluster_pheno) <- c('Currently on ICS', 'High dose ICS', 'Cluster')
###############################################################


###############################################################
SARP_3_cluster_pheno_long <- na.omit(SARP_3_cluster_pheno %>% pivot_longer(!Cluster, names_to = "dosage_name", values_to = "dosage_value"))
SARP3_pheno_summary <- as.data.frame(SARP_3_cluster_pheno_long %>% group_by(Cluster, dosage_name) %>% summarize(count = sum(dosage_value)))
SARP3_pheno_summary$frac <- SARP3_pheno_summary$count / dim(SARP_3_cluster_pheno)[1]

# Create contingency tables for chi square test
SARP_3_cluster_pheno_long$cluster_severe <- "1-2"
SARP_3_cluster_pheno_long$cluster_severe[SARP_3_cluster_pheno_long$Cluster %in% c(3, 4)] <- "3-4"
SARP3_pheno_contingency <- as.data.frame(na.omit(SARP_3_cluster_pheno_long) %>% group_by(cluster_severe, dosage_name) %>% summarize(Yes = sum(dosage_value), No = sum(1- dosage_value)))
###############################################################


###############################################################
compute_chi_sq_plot_frac <- function(dosage_to_plot, y_label_pos)
{
  SARP3_pheno_contingency <- SARP3_pheno_contingency[SARP3_pheno_contingency$dosage_name == dosage_to_plot, c("cluster_severe", "Yes", "No")]
  rownames(SARP3_pheno_contingency) <- SARP3_pheno_contingency$cluster_severe
  SARP3_pheno_contingency <- subset(SARP3_pheno_contingency, select = -c(cluster_severe))
  ICS_p_val <- pairwise_prop_test(t(SARP3_pheno_contingency))
  p_val <- paste0("p = ", round(ICS_p_val$p, 4))
  ICS_summary <- ggplot(SARP3_pheno_summary[SARP3_pheno_summary$dosage_name == dosage_to_plot, ], aes(x = Cluster, y = frac, fill = Cluster)) +
    xlab("Cluster") + ylab("Fraction of subjects") + ggtitle(dosage_to_plot) +
    geom_bar(stat = "identity", width = 0.8) +
    theme(legend.position = "none") +
    annotate(geom = "text", label = p_val, x = 1.5, y = y_label_pos,
      size = text_size_labels, colour = black_text_colour) +
    scale_fill_manual(values = cluster_colors)
  return(ICS_summary)
}
ICS_plot <- plot_grid(compute_chi_sq_plot_frac("Currently on ICS", 0.2), compute_chi_sq_plot_frac("High dose ICS", 0.15), align = "hv", axis = "tblr", ncol = 1)
###############################################################


###############################################################
cairo_pdf(paste0(figures_folder, "fig_3_conorm_expr_heatmap_dosage_effect.pdf"), width = 7.2, height = 5.5, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 100, ncol = 100)))

pushViewport(viewport(layout.pos.row = 2:99, layout.pos.col = 1:62))
draw(ht_mp_expr, height = length(genes_of_interest) * tile_height, gap = gap_size, background = "transparent", newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 15:23, layout.pos.col = 63:93))
draw(expr_legend)
upViewport(1)

pushViewport(viewport(layout.pos.row = 2:6, layout.pos.col = 63:79))
draw(cluster_legend)
upViewport(1)

pushViewport(viewport(layout.pos.row = 8:13, layout.pos.col = 63:94))
draw(severity_legend)
upViewport(1)

pushViewport(viewport(layout.pos.row = 2:6, layout.pos.col = 80:96))
draw(sex_legend)
upViewport(1)

print(ICS_plot, vp = viewport(layout.pos.row = 28:99, layout.pos.col = 64:99))

grid.text(x = unit(0.02, "npc"), y = unit(0.98, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.635, "npc"), y = unit(0.72, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.635, "npc"), y = unit(0.36, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))

dev.off()
###############################################################
