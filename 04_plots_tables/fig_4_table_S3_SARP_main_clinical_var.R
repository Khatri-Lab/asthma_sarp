###############################################################
library(tidyverse)
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/plot_themes.R")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_pheno_cluster.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/SARP_meta_objs.rds")
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
SARP_3_cluster_pheno <- merge(SARP_3_cluster_pheno, SARP_3$pheno, by = "participant")
SARP_3_cluster_pheno <- SARP_3_cluster_pheno[, all_variables]
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
pheno_value_cluster_summary <- SARP_3_cluster_pheno_plot %>% group_by(pheno_name, Cluster) %>%
  summarize(mean = mean(pheno_value, na.rm = TRUE), sd = sd(pheno_value, na.rm = TRUE), na = sum(is.na(pheno_value)))
adjusted_kruskal_test$fdr <- p.adjust(adjusted_kruskal_test$p, method = "fdr")
na_summary <- pheno_value_cluster_summary %>% group_by(pheno_name) %>% summarize(na = sum(na))
pheno_value_cluster_summary$cluster_summary <- paste0(round(pheno_value_cluster_summary$mean), " (", round(pheno_value_cluster_summary$sd), ")")
pheno_value_cluster_summary$Cluster <- paste0("Cluster ", pheno_value_cluster_summary$Cluster)
pheno_value_cluster_summary_wide <- pheno_value_cluster_summary[, c("pheno_name", "Cluster", "cluster_summary")] %>% pivot_wider(names_from = "Cluster", values_from = "cluster_summary")
combined_cluster_summary_for_printing <- merge(pheno_value_cluster_summary_wide, na_summary, by = "pheno_name")
combined_cluster_summary_for_printing <- merge(combined_cluster_summary_for_printing, adjusted_kruskal_test[, c("pheno_name", "p", "fdr")], by = "pheno_name")
###############################################################


###############################################################
variables_to_remove <- union(adjusted_kruskal_test$pheno_name[adjusted_kruskal_test$fdr > 0.1], pheno_value_cluster_summary$pheno_name[pheno_value_cluster_summary$sd == 0])
adjusted_wilcox_test <- SARP_3_cluster_pheno_plot[!(SARP_3_cluster_pheno_plot$pheno_name %in% variables_to_remove), ] %>% group_by(pheno_name) %>%
  wilcox_test(pheno_value ~ Cluster, comparisons = list(c("1", "2"), c("1", "3"), c("1", "4"), c("2", "3"), c("2", "4"), c("3", "4")))
adjusted_wilcox_test <- adjusted_wilcox_test %>% group_by(pheno_name) %>% mutate(fdr = p.adjust(p, method = "fdr"))
adjusted_wilcox_test$fdr <- format(adjusted_wilcox_test$fdr, digits = 3, scientific = TRUE)
adjusted_wilcox_test$comparison <- paste0(adjusted_wilcox_test$group1, " vs. ", adjusted_wilcox_test$group2, " (FDR)")
adjusted_wilcox_test <- adjusted_wilcox_test[, c("pheno_name", "comparison", "fdr")] %>% pivot_wider(names_from = "comparison", values_from = "fdr")
combined_cluster_summary_for_printing <- merge(combined_cluster_summary_for_printing, adjusted_wilcox_test, by = "pheno_name", all.x = TRUE)
###############################################################


###############################################################
combined_cluster_summary_for_printing$p <- format(combined_cluster_summary_for_printing$p, digits = 3, scientific = TRUE)
combined_cluster_summary_for_printing$fdr <- format(combined_cluster_summary_for_printing$fdr, digits = 3, scientific = TRUE)
colnames(combined_cluster_summary_for_printing) <- c("Variable", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "# NA", "p (Kruskal-Wallis)", "FDR (Kruskal-Wallis)", "1 vs. 2 (FDR)", "1 vs. 3 (FDR)", "1 vs. 4 (FDR)", "2 vs. 3 (FDR)", "2 vs. 4 (FDR)", "3 vs. 4 (FDR)")
combined_cluster_summary_for_printing <- combined_cluster_summary_for_printing %>% arrange(match(Variable, names(all_variables)))
write_csv(combined_cluster_summary_for_printing, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/table_S3_SARP_clinical_variable_summary.csv")
###############################################################


###############################################################
y_pos_label <- SARP_3_cluster_pheno_plot %>% group_by(pheno_name) %>%
  summarize(y = max(pheno_value, na.rm = TRUE) * 1.05, x = 2.5, Cluster = "2")
adjusted_kruskal_test <- merge(adjusted_kruskal_test, y_pos_label, by = "pheno_name")
adjusted_kruskal_test$fdr_print <- paste0("fdr (KW) = ", format(adjusted_kruskal_test$fdr, digits = 2))
variables_to_plot <- c("FEV1 (% pred)", "FEV1/FVC ratio", "Exacerbations", "Weekly SABA", "Nighttime awakenings", "Max. Bronchodilator reversibility", "Methacholine challenge", "FeNO (ppb)", "IgE (log10(kU/L))", "Sputum Th2 GM score", "Sputum Eosinophils (%)", "Blood Eosinophils (%)")
SARP_3_cluster_pheno_plot_selected <- transform(SARP_3_cluster_pheno_plot[SARP_3_cluster_pheno_plot$pheno_name %in% variables_to_plot, ],
  pheno_name = factor(pheno_name, levels = variables_to_plot))
names(cluster_colors) <- c("1", "2", "3", "4")
SARP_3_cluster_pheno_plot$Cluster <- as.factor(SARP_3_cluster_pheno_plot$Cluster)
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
ggsave("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/fig_4_clinical_vars.pdf", clinical_var_plot, width = 7.2, height = 7.2, device = cairo_pdf)
###############################################################