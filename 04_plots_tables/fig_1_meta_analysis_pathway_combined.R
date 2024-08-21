###############################################################
library(readr)
library(data.table)
library(MetaIntegrator)
library(tidyverse)
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/plot_themes.R")
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/JT_test_function.R")
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/pooled_ROC_plot.R")
###############################################################


###############################################################
tile_height <- unit(base_text_size, "pt")
tile_width <- unit(0.8, "pt")
ht_opt("simple_anno_size" = tile_height)
figures_folder <- "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/"
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/discovery_validation_bayes_meta_objs.rds")
###############################################################


###############################################################
min_eff_size <- -3.4
max_eff_size <- 3.4
eff_size_matrix_colour_pts <- c(min_eff_size, min_eff_size/2, 0, max_eff_size/2, max_eff_size)
eff_size_cols <- colorRamp2(eff_size_matrix_colour_pts, rev(scico(5, palette = "bam")))
eff_size_matrix_legend_pts <- c(min_eff_size, max_eff_size)
eff_size_matrix_legend_labels <- c(round(min_eff_size, 2), round(max_eff_size, 2))
eff_size_legend <- Legend(col_fun = eff_size_cols,
  at = eff_size_matrix_legend_pts,
  labels = eff_size_matrix_legend_labels,
  title = "Effect\nsize", title_position = "topcenter",
  grid_width = unit(base_text_size / 3, "pt"),
  legend_height = 4 * unit(base_text_size, "pt"),
  background = "transparent",
  title_gp = text_setting_small,
  title_gap = gap_size, gap = gap_size,
  labels_gp = text_setting_sub_small,
  direction = "vertical")
gene_order_frame <- discovery_meta_obj$bayesianMeta$finalResults[discovery_meta_obj$bayesianMeta$finalResults$Gene %in% sig_genes$gene_sig, c("Gene", "ES")]
gene_order <- gene_order_frame$Gene[order(gene_order_frame$ES)]
gene_order <- gsub("-", ".", gene_order)
plot_eff_size_heatmap <- function(meta_obj, plot_title)
{
  eff_size_matrix <- merge(meta_obj$bayesianMeta$datasetEffectSizes[meta_obj$bayesianMeta$datasetEffectSizes$Gene %in% sig_genes$gene_sig, colnames(meta_obj$bayesianMeta$datasetEffectSizes) != "nStudies"], meta_obj$bayesianMeta$finalResults[meta_obj$bayesianMeta$finalResults$Gene %in% sig_genes$gene_sig, c("Gene", "ES")], by = "Gene", all = TRUE)
  colnames(eff_size_matrix)[ncol(eff_size_matrix)] <- "Pooled"
  rownames(eff_size_matrix) <- eff_size_matrix$Gene
  eff_size_matrix <- data.frame(t(eff_size_matrix[, colnames(eff_size_matrix) != "Gene"]))
  if(sum(!(gene_order %in% colnames(eff_size_matrix))) > 0)
  {
    eff_size_matrix[, gene_order[!(gene_order %in% colnames(eff_size_matrix))]] <- NA
  }
  eff_size_matrix <- as.matrix(eff_size_matrix[, gene_order])
  ht_mp_eff_size <- Heatmap(eff_size_matrix,
    col = eff_size_cols, na_col = "#9A9A9A",
    show_column_names = FALSE,
    width = ncol(eff_size_matrix) * tile_width,
    name = plot_title, column_title = plot_title,
    column_title_side = "top", column_title_gp = text_setting_small,
    row_names_side = "right",
    show_heatmap_legend = FALSE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    gap = gap_size)
  return(ht_mp_eff_size)
}
discovery_ht_mp <- plot_eff_size_heatmap(discovery_meta_obj, "Differentially expressed genes, Discovery")
validation_ht_mp <- plot_eff_size_heatmap(validation_meta_obj, "Differentially expressed genes, Validation")
###############################################################


###############################################################
#names(validation_meta_obj$originalData) <- names(validation_data_colours)
ROC_discovery <- IL_pooledROCPlot(discovery_meta_obj, sig_genes, custom_colors = discovery_data_colours, title = "Discovery") + theme(legend.background = element_rect(fill = "#FFFFFFEE", linewidth = 0, colour = "transparent"), panel.grid.major = element_line(linewidth = line_size_text_repel_mm))
ROC_validation <- IL_pooledROCPlot(validation_meta_obj, sig_genes, custom_colors = validation_data_colours, title = "Validation") + theme(legend.background = element_rect(fill = "#FFFFFFEE", linewidth = 0, colour = "transparent"), panel.grid.major = element_line(linewidth = line_size_text_repel_mm))
ROC_combined <- plot_grid(ROC_discovery, ROC_validation, nrow = 1, align = "hv", axis = "tblr")
###############################################################


###############################################################
compute_sample_scores <- function(meta_obj, split_name)
{
  pheno_score_frame <- Reduce(rbind, lapply(meta_obj$originalData, function(x) data.frame(Score = calculateScore(sig_genes, x, TRUE), Severity = x$pheno$severity)))
  pheno_score_frame$Severity <- ifelse(grepl("Mild|moderate|Moderate", pheno_score_frame$Severity), "Mild/moderate", pheno_score_frame$Severity)
  pheno_score_frame$Severity <- ifelse(grepl("Severe", pheno_score_frame$Severity), "Severe", pheno_score_frame$Severity)
  pheno_score_frame$Severity <- ifelse(grepl("healthy", pheno_score_frame$Severity), "Healthy", pheno_score_frame$Severity)
  pheno_score_frame$Severity <- ifelse(grepl("Unspecified", pheno_score_frame$Severity), "Unspecified", pheno_score_frame$Severity)
  pheno_score_frame$Severity <- ifelse(grepl("rhinitis", pheno_score_frame$Severity), "Healthy", pheno_score_frame$Severity)
  pheno_score_frame[pheno_score_frame$Severity == "Unspecified", ] <- NA
  pheno_score_frame <- na.omit(pheno_score_frame)
  pheno_score_frame$Severity <- ordered(pheno_score_frame$Severity, levels = c("Healthy", "Mild/moderate", "Severe"))
  pheno_score_frame$Split <- split_name
  jt_stat <- JT.test(data = pheno_score_frame$Score, class = pheno_score_frame$Severity, labs = c("Healthy", "Mild/moderate", "Severe"))
  p_value_print <- paste0("p(JT) = ", format(jt_stat$p.value, digits = 2, scientific = TRUE))
  if(jt_stat$p.value == 0)
  {
    p_value_print <- "p(JT) < 2.2e-16"
  }
  pheno_score_frame$jt_p <- p_value_print
  return(pheno_score_frame)
}
combined_bec_scores <- rbind(compute_sample_scores(discovery_meta_obj, "Discovery"), compute_sample_scores(validation_meta_obj, "Validation"))
p_value_print <- combined_bec_scores %>% group_by(Split) %>% summarise(jt_p = first(jt_p), Severity = "Healthy")
bec_score_summary <- ggplot(combined_bec_scores, aes(x = Severity, y = Score, fill = Severity)) +
  ylab("Scaled BEC score") +
  geom_violin(lwd = line_size_text_repel_mm, colour = black_text_colour) +
  geom_beeswarm(cex = cex_val * 1.1, size = point_size, colour = black_text_colour, show.legend = FALSE) +
  geom_boxplot(fill = NA, lwd = line_size_main_mm, width = 0.4, outlier.shape = NA, colour = black_text_colour) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
    legend.position = "bottom", legend.direction = "horizontal", legend.margin = margin(1, 1, 1, 1), legend.box.margin = margin(-10, 5, 3, 3),
    panel.grid.major.y = element_line(linewidth = line_size_text_repel_mm)) +
  geom_text(data = p_value_print, aes(label = jt_p),
    x = -Inf, y = max(combined_bec_scores$Score),
    size = text_size_labels, hjust = -0.25, vjust = 0.75, colour = black_text_colour) +
  scale_fill_manual(values = severity_colours) +
  labs(fill = NULL) + guides(fill = guide_legend(nrow = 1, override.aes = list(size = text_size_labels * 6 / 9))) +
  facet_wrap(vars(Split), nrow = 1)
###############################################################


###############################################################
gene_set_enrichment <- read_csv("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/ArtZEnrichment.csv")
gene_set_enrichment <- gene_set_enrichment[order(gene_set_enrichment$direction, decreasing = TRUE), ]
gene_set_enrichment$figure.set.name <- gsub("_", " ", gene_set_enrichment$set.name.short)
gene_set_enrichment <- gene_set_enrichment[gene_set_enrichment$padj < 0.05, ]
###############################################################


###############################################################
pathway_barplot <- ggplot(gene_set_enrichment, aes(x = figure.set.name, y = direction, fill = padj)) +
  geom_bar(stat = "identity") +
  #scale_fill_gradient(low = "#2C194C", high = "#CAC6F4") +
  scale_fill_gradient(low = "#14385A", high = "#A2B0CA") +
  scale_x_discrete(limits = rev(c(gene_set_enrichment$figure.set.name))) +
  ylab("Direction") +
  theme(axis.title.y = element_blank(), legend.position = c(.8, 0.25),
    panel.grid.major = element_line(linewidth = line_size_text_repel_mm)) +
  labs(fill = "Adjusted\np-value") +
  coord_flip()
###############################################################


###############################################################
cairo_pdf(paste0(figures_folder, "fig_1_meta_analysis_pathway.pdf"), width = 7.2, height = 9.7, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 200, ncol = 100)))

pushViewport(viewport(layout.pos.row = 1:27, layout.pos.col = 2:93))
draw(discovery_ht_mp, height = (length(discovery_data_colours) + 1) * tile_height, gap = gap_size, background = "transparent", newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 28:44, layout.pos.col = 1:93))
draw(validation_ht_mp, height = (length(validation_data_colours) + 1) * tile_height, gap = gap_size, background = "transparent", newpage = FALSE)
upViewport(1)

pushViewport(viewport(layout.pos.row = 1:27, layout.pos.col = 94:99))
draw(eff_size_legend)
upViewport(1)

print(ROC_combined, vp = viewport(layout.pos.row = 45:95, layout.pos.col = 2:99))
print(bec_score_summary, vp = viewport(layout.pos.row = 96:139, layout.pos.col = 2:99))
print(pathway_barplot, vp = viewport(layout.pos.row = 140:200, layout.pos.col = 2:99))

grid.text(x = unit(0.02, "npc"), y = unit(0.99, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.02, "npc"), y = unit(0.85, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.02, "npc"), y = unit(0.765, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.53, "npc"), y = unit(0.765, "npc"), label = "D", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.02, "npc"), y = unit(0.505, "npc"), label = "E", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.53, "npc"), y = unit(0.505, "npc"), label = "F", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.02, "npc"), y = unit(0.305, "npc"), label = "G", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))

dev.off()
###############################################################

