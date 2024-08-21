###############################################################
library(tidyverse)
library(stringr)
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/plot_themes.R")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/discovery_validation_bayes_meta_objs.rds")
###############################################################


###############################################################
tsai_signature <- read_csv("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/tsai_signature.csv")
tsai_up_genes <- tsai_signature$`Gene Symbol`[tsai_signature$`Hedges g` > 0] %>% sub(pattern = " ", replacement = "") %>% strsplit("///") %>% unlist()
tsai_down_genes <- tsai_signature$`Gene Symbol`[tsai_signature$`Hedges g` < 0] %>% sub(pattern = " ", replacement = "") %>% strsplit("///") %>% unlist()
length(intersect(tsai_up_genes, sig_genes$posGeneNames))#153
length(intersect(tsai_up_genes, sig_genes$negGeneNames))#0
length(intersect(tsai_down_genes, sig_genes$posGeneNames))#0
length(intersect(tsai_down_genes, sig_genes$negGeneNames))#101
###############################################################


schofield_eosniophil <-


###############################################################
names(cluster_colors) <- c("1", "2", "3", "4")
plot_cluster <- ggplot(combined_pheno, aes(x = UMAP_1, y = UMAP_2)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  geom_point(aes(colour = as.factor(Cluster)), size = point_size) +
  theme(legend.position = "none") +
  scale_colour_manual(values = cluster_colors)
plot_cluster_marginal <- ggMarginal(plot_cluster, groupColour = TRUE, groupFill = TRUE)
###############################################################


###############################################################
plot_severity <- ggplot(combined_pheno, aes(x = UMAP_1, y = UMAP_2)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  geom_point(aes(colour = Severity), size = point_size) +
  theme(legend.position = "none") +
  scale_colour_manual(values = severity_colours)
plot_severity_marginal <- ggMarginal(plot_severity, groupColour = TRUE, groupFill = TRUE)
###############################################################


###############################################################
cluster_summary <- combined_pheno %>% group_by(Cluster, Severity) %>% summarize(Count = n()) %>% mutate(Prop = Count / sum(Count))
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
# Test for differences across clusters for the following clinical measurements from SARP
severity_var <- c("preFEV1_Quanpct_pred", "max_fev1fvc", "WG_derived_num_exac12", "rxhx_1040", "acq_1")
names(severity_var) <- c("FEV1 (% pred)", "FEV1/FVC ratio", "Exacerbations",  "Weekly SABA", "Nighttime awakenings")

diagnostic_var <- c("pc20_comb", "mrv_max_chg")
names(diagnostic_var) <- c("Methacholine challenge", "Max. Bronchodilator")

biomarker_var <- c("feno", "ige_log", "Th2GM")
names(biomarker_var) <- c("FeNO (ppb)", "IgE (log10(kU/L))", "Sputum Th2GM")

eosinophil_var <- c("eos_count", "sptr_1090", "bal_1040")
names(eosinophil_var) <- c("Blood eosinophil (cells/\U00B5L)", "Sputum eosinophil (%)", "BAL eosinophil (%)")

demographic_var <- c("age_v2", "age_onset", "bmi")
names(demographic_var) <- c("Age (yrs)", "Age Asthma Onset (yrs)", "BMI")

all_variables <- c(Cluster = "Cluster", severity_var, diagnostic_var, biomarker_var, eosinophil_var, demographic_var)
###############################################################


###############################################################
SARP_3_cluster_pheno <- combined_pheno[combined_pheno$Dataset == "SARP 3", ]
modified_subj_id <- SARP_3_cluster_pheno %>% rownames() %>% sub(pattern = "X", replacement = "") %>% sub(pattern = "\\.", replacement = "-") %>% sub(pattern = "\\.", replacement = "-")
SARP_3_cluster_pheno$participant <- as.character(modified_subj_id)
SARP_3_cluster_pheno_full <- merge(SARP_3_cluster_pheno, SARP_3$pheno, by = "participant")
SARP_3_cluster_pheno <- SARP_3_cluster_pheno_full[, all_variables]
colnames(SARP_3_cluster_pheno) <- names(all_variables)
SARP_3_cluster_pheno_plot <- SARP_3_cluster_pheno %>% pivot_longer(!Cluster, names_to = "pheno_name", values_to = "pheno_value")
SARP_3_cluster_pheno_plot$pheno_value[SARP_3_cluster_pheno_plot$pheno_name == "BAL eosinophil (%)" & SARP_3_cluster_pheno_plot$pheno_value > 10] <- NA
###############################################################


###############################################################
adjusted_kruskal_test <- SARP_3_cluster_pheno_plot %>% group_by(pheno_name) %>%
  kruskal_test(pheno_value ~ as.factor(Cluster))
y_pos_label <- SARP_3_cluster_pheno_plot %>% group_by(pheno_name) %>%
  summarize(y = max(pheno_value, na.rm = TRUE), x = 2.5, Cluster = "2")
adjusted_kruskal_test <- merge(adjusted_kruskal_test, y_pos_label, by = "pheno_name")
adjusted_kruskal_test$fdr <- p.adjust(adjusted_kruskal_test$p, method = "fdr")
adjusted_kruskal_test$fdr_print <- paste0("fdr (Kruskal) = ", format(adjusted_kruskal_test$fdr, digits = 2))
###############################################################


###############################################################
SARP_3_cluster_pheno_plot <- transform(SARP_3_cluster_pheno_plot,
  pheno_name = factor(pheno_name, levels = c("FEV1 (% pred)", "FEV1/FVC ratio", "Exacerbations",  "Weekly SABA", "Nighttime awakenings", "Methacholine challenge", "Max. Bronchodilator", "FeNO (ppb)", "IgE (log10(kU/L))", "Blood eosinophil (cells/\U00B5L)", "Sputum eosinophil (%)", "BAL eosinophil (%)", "Sputum Th2GM", "Age (yrs)", "Age Asthma Onset (yrs)", "BMI")))
names(cluster_colors) <- c("1", "2", "3", "4")
SARP_3_cluster_pheno_plot$Cluster <- as.factor(SARP_3_cluster_pheno_plot$Cluster)
clinical_var_plot <- ggplot(SARP_3_cluster_pheno_plot, aes(x = Cluster, y = pheno_value, fill = Cluster)) +
  xlab("Cluster") +
  geom_boxplot(lwd = line_size_main_mm, width = 0.5, colour = black_text_colour, outlier.shape = NA) +
  geom_beeswarm(cex = cex_val, size = point_size, colour = black_text_colour) +
  scale_fill_manual(values = cluster_colors) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  geom_text(data = adjusted_kruskal_test, mapping = aes(label = fdr_print, x = x, y = y),# x = "2", y = y_pos_label),
    size = text_size_labels, colour = black_text_colour,
    vjust = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  facet_wrap(pheno_name ~ ., nrow = 4, scale = "free_y")
###############################################################


###############################################################
SARP_3_cluster_pheno <- subset(SARP_3_cluster_pheno_full, select = c('rxhx_1120', 'rxhx_1232', 'Cluster'))
colnames(SARP_3_cluster_pheno) <- c('Currently on ICS', 'High dose ICS', 'Cluster')
SARP_3_cluster_pheno_long <- na.omit(SARP_3_cluster_pheno %>% pivot_longer(!Cluster, names_to = "dosage_name", values_to = "dosage_value"))
SARP3_pheno_summary <- as.data.frame(SARP_3_cluster_pheno_long %>% group_by(Cluster, dosage_name) %>% summarize(count = sum(dosage_value)))
SARP3_pheno_summary$frac <- SARP3_pheno_summary$count / dim(SARP_3_cluster_pheno)[1]

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
  ICS_summary <- ggplot(SARP3_pheno_summary[SARP3_pheno_summary$dosage_name == dosage_to_plot, ], aes(x = as.factor(Cluster), y = frac, fill = as.factor(Cluster))) +
    xlab("Cluster") + ylab("Fraction of subjects") + ggtitle(dosage_to_plot) +
    geom_bar(stat = "identity", width = 0.8) +
    theme(legend.position = "none") +
    # geom_text(label = p_val, x = 1.5, y = y_label_pos,
    #   size = text_size_labels, colour = black_text_colour) +
    scale_fill_manual(values = cluster_colors)
  return(ICS_summary)
}
ICS_plot <- plot_grid(compute_chi_sq_plot_frac("Currently on ICS", 0.2), compute_chi_sq_plot_frac("High dose ICS", 0.15), align = "hv", axis = "tblr", nrow = 1)
###############################################################


###############################################################
figures_folder <- "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures_reviewer/"
umap_combined <- plot_grid(plot_cluster_marginal, plot_severity_marginal, nrow = 1, align = "hv", axis = "tblr")
bars_combined <- plot_grid(cluster_hist, severity_hist, nrow = 1, align = "hv", axis = "tblr")

cairo_pdf(paste0(figures_folder, "clustering_by_all_non_sig_genes.pdf"), width = 7.2, height = 15.2, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 152, ncol = 100)))

print(umap_combined, vp = viewport(layout.pos.row = 2:30, layout.pos.col = 2:99))
print(bars_combined, vp = viewport(layout.pos.row = 31:50, layout.pos.col = 2:99))
print(clinical_var_plot, vp = viewport(layout.pos.row = 51:122, layout.pos.col = 2:99))
print(ICS_plot, vp = viewport(layout.pos.row = 123:152, layout.pos.col = 2:99))

# grid.text(x = unit(0.02, "npc"), y = unit(0.98, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
# grid.text(x = unit(0.515, "npc"), y = unit(0.98, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
# grid.text(x = unit(0.02, "npc"), y = unit(0.385, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
# grid.text(x = unit(0.515, "npc"), y = unit(0.385, "npc"), label = "D", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))

dev.off()
###############################################################

