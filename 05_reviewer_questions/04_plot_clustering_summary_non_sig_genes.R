###############################################################
library(tidyverse)
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/plot_themes.R")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_with_all_non_sig_genes_pheno_with_umap_clusters.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/SARP_meta_objs.rds")
###############################################################


###############################################################
combined_pheno$Cluster <- 5 - combined_pheno$Cluster
figures_folder <- "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures_reviewer/"
###############################################################


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
severity_var <- c("preFEV1_Quanpct_pred", "max_fev1fvc", "WG_derived_num_exac12", "rxhx_1040", "acq_1")
names(severity_var) <- c("FEV1 (% pred)", "FEV1/FVC ratio", "Exacerbations",  "Weekly SABA", "Nighttime awakenings")

diagnostic_var <- c("pc20_comb", "mrv_max_chg", "scr2_1030")
names(diagnostic_var) <- c("Methacholine challenge", "Max. Bronchodilator reversibility", "Bronchodilator response (%)")

biomarker_var <- c("feno", "ige_log", "one_specIgE", "num_specIgE", "Th2GM")
names(biomarker_var) <- c("FeNO (ppb)", "IgE (log10(kU/L))", "Any specific IgE test", "Number of specific IgE tests", "Sputum Th2 GM score")

demographic_var <- c("age_v2", "age_onset", "bmi", "Sex")
names(demographic_var) <- c("Age (yrs)", "Age Asthma Onset (yrs)", "BMI", "Sex (% male)")

bal_var <- c("bal_1030", "bal_1010", "bal_1020", "bal_1040")
names(bal_var) <- c("BAL Neutrophils (%)", "BAL Lymphocytes (%)", "Alveolar macrophages (%)", "BAL Eosinophils (%)")

blood_count_var <- c("neut_count", "lymph_count", "mono_count", "eos_count")
names(blood_count_var) <- c("Blood Neutrophils (cells/\U00B5L)", "Blood Lymphocytes (cells/\U00B5L)", "Blood Monocytes (cells/\U00B5L)", "Blood Eosinophils (cells/\U00B5L)")

blood_var <- c("lab_1010", "lab_1020", "lab_1030", "lab_1040")
names(blood_var) <- c("Blood Neutrophils (%)", "Blood Lymphocytes (%)", "Blood Monocytes (%)", "Blood Eosinophils (%)")

sputum_var <- c("sptr_1080", "sptr_1100","sptr_1070", "sptr_1090")
names(sputum_var) <- c("Sputum Neutrophil (%)", "Sputum Lymphocytes (%)",  "Sputum Macrophage (%)", "Sputum Eosinophils (%)")

allergic_var <- c("allergies", "allergic_envallergen", "allergic_foods", "eczema_ever", "asthma_mother", "asthma_father", "asthma_child", "asthma_sibling", "atopic_mother", "atopic_father", "atopic_child", "atopic_sibling")
names(allergic_var) <- c("Allergies diagnosis", "Environmental allergy", "Food allergy", "History of eczema", "Maternal asthma", "Paternal asthma", "Child with asthma", "Sibling with asthma", "Atopic mother", "Atopic father", "Atopic child", "Atopic sibling")

provoke_var <- c("provoke_exercise", "provoke_stairs", "provoke_nsaids", "provoke_uri", "provoke_irritants", "provoke_weather", "provoke_cold", "provoke_emotional", "provoke_tobacco", "provoke_preserv", "provoke_envallergen")
names(provoke_var) <- c("Exercise (provoking)", "Stairs (provoking)", "NSAIDs (provoking)", "URI (provoking)", "Irritants (provoking)", "Weather (provoking)", "Cold (provoking)", "Emotions (provoking)", "Tobacco (provoking)", "Preservatives (provoking)", "Environmental allergens (provoking)")
#provoke_var <- provoke_var[c(1,3:9,11)] #selecting for sepcific provoking factors

smoke_var <- c("smoke2_gestation", "smoke2_home", "smoke2_day" )
names(smoke_var) <- c("Smoke in utero", "Smoke in childhood", "Current smoke exposure")

corticosteroids_var <- c("rxhx_1120", "rxhx_1232")
names(corticosteroids_var) <- c("Currently on ICS", "High dose ICS")

all_variables <- c(Cluster = "Cluster", demographic_var, severity_var, diagnostic_var, biomarker_var, bal_var, blood_count_var, blood_var, sputum_var, allergic_var, provoke_var, smoke_var, corticosteroids_var)
###############################################################


###############################################################
SARP_3_cluster_pheno <- combined_pheno[combined_pheno$Dataset == "SARP 3", ]
modified_subj_id <- SARP_3_cluster_pheno %>% rownames() %>% sub(pattern = "X", replacement = "") %>% sub(pattern = "\\.", replacement = "-") %>% sub(pattern = "\\.", replacement = "-")
SARP_3_cluster_pheno$participant <- as.character(modified_subj_id)
SARP_3_cluster_pheno_full <- merge(SARP_3_cluster_pheno, SARP_3$pheno, by = "participant")
SARP_3_cluster_pheno <- SARP_3_cluster_pheno_full[, all_variables]
colnames(SARP_3_cluster_pheno) <- names(all_variables)
###############################################################


###############################################################
SARP_3_cluster_pheno$`Sex (% male)` <- ifelse(SARP_3_cluster_pheno$`Sex (% male)` == "male", 1, 0)
SARP_3_cluster_pheno$`Sex (% male)` <- as.numeric(SARP_3_cluster_pheno$`Sex (% male)`) * 100
SARP_3_cluster_pheno$`Bronchodilator response (%)` <- SARP_3_cluster_pheno$`Bronchodilator response (%)` * 100
SARP_3_cluster_pheno_plot <- SARP_3_cluster_pheno %>% pivot_longer(!Cluster, names_to = "pheno_name", values_to = "pheno_value")
###############################################################


###############################################################
adjusted_kruskal_test <- SARP_3_cluster_pheno_plot %>% group_by(pheno_name) %>%
  kruskal_test(pheno_value ~ as.factor(Cluster))
y_pos_label <- SARP_3_cluster_pheno_plot %>% group_by(pheno_name) %>%
  summarize(y = max(pheno_value, na.rm = TRUE) * 1.05, x = 2.5, Cluster = "2")
adjusted_kruskal_test <- merge(adjusted_kruskal_test, y_pos_label, by = "pheno_name")
adjusted_kruskal_test$fdr <- p.adjust(adjusted_kruskal_test$p, method = "fdr")
adjusted_kruskal_test$fdr_print <- paste0("fdr (KW) = ", format(adjusted_kruskal_test$fdr, digits = 2))
###############################################################


###############################################################
variables_to_plot <- c("FEV1 (% pred)", "FEV1/FVC ratio", "Exacerbations", "Weekly SABA", "Nighttime awakenings", "Max. Bronchodilator reversibility", "Methacholine challenge", "FeNO (ppb)", "IgE (log10(kU/L))", "Sputum Th2 GM score", "Sputum Eosinophils (%)", "Blood Eosinophils (%)")
SARP_3_cluster_pheno_plot_selected <- transform(SARP_3_cluster_pheno_plot[SARP_3_cluster_pheno_plot$pheno_name %in% variables_to_plot, ],
  pheno_name = factor(pheno_name, levels = variables_to_plot))
names(cluster_colors) <- c("1", "2", "3", "4")
SARP_3_cluster_pheno_plot_selected$Cluster <- as.factor(SARP_3_cluster_pheno_plot_selected$Cluster)
clinical_var_plot <- ggplot(SARP_3_cluster_pheno_plot_selected, aes(x = Cluster, y = pheno_value, fill = Cluster)) +
  xlab("Cluster") +
  geom_boxplot(lwd = line_size_main_mm, width = 0.5, colour = black_text_colour, outlier.shape = NA) +
  geom_beeswarm(cex = cex_val, size = point_size, colour = black_text_colour) +
  scale_fill_manual(values = cluster_colors) +
  theme(legend.position = "none", axis.title.y = element_blank()) +
  geom_text(data = adjusted_kruskal_test[adjusted_kruskal_test$pheno_name %in% variables_to_plot, ], mapping = aes(label = fdr_print, x = x, y = y),
    size = text_size_labels, colour = black_text_colour,
    vjust = 0) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  facet_wrap(as.factor(pheno_name) ~ ., nrow = 4, scale = "free_y")
ggsave(paste0(figures_folder, "fig_R2_clinical_variables_summary_by_non_sig_genes.pdf"), clinical_var_plot, width = 7.2, height = 7.2, bg = "transparent", device = cairo_pdf)
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
umap_combined <- plot_grid(plot_cluster_marginal, plot_severity_marginal, nrow = 1, align = "hv", axis = "tblr")
bars_combined <- plot_grid(cluster_hist, severity_hist, nrow = 1, align = "hv", axis = "tblr")

cairo_pdf(paste0(figures_folder, "fig_R1_clustering_by_all_non_sig_genes.pdf"), width = 7.2, height = 7.5, bg = "transparent")

pushViewport(viewport(layout = grid.layout(nrow = 75, ncol = 100)))

print(umap_combined, vp = viewport(layout.pos.row = 2:30, layout.pos.col = 2:99))
print(bars_combined, vp = viewport(layout.pos.row = 31:50, layout.pos.col = 2:99))
print(ICS_plot, vp = viewport(layout.pos.row = 51:75, layout.pos.col = 2:99))

grid.text(x = unit(0.02, "npc"), y = unit(0.98, "npc"), label = "A", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.515, "npc"), y = unit(0.98, "npc"), label = "B", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.02, "npc"), y = unit(0.62, "npc"), label = "C", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.515, "npc"), y = unit(0.62, "npc"), label = "D", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.02, "npc"), y = unit(0.33, "npc"), label = "E", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))
grid.text(x = unit(0.515, "npc"), y = unit(0.33, "npc"), label = "F", gp = gpar(fontsize = base_text_size + 2, fontface = "bold", col = black_text_colour))

dev.off()
###############################################################
