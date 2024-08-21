###############################################################
library(tidyverse)
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_pheno_cluster.rds")
###############################################################


###############################################################
combined_summary <- combined_pheno %>% group_by(Cluster) %>% summarize(overall_total = n(), Female = sum(Sex == "female", na.rm = TRUE), Asthma = sum(Severity != "Healthy", na.rm = TRUE)) %>% pivot_longer(c("Female", "Asthma"), names_to = "Demographic", values_to = "overall_value")
SARP_3_summary <- combined_pheno[combined_pheno$Dataset == "SARP 3", ] %>% group_by(Cluster) %>% summarize(sarp_total = n(), Female = sum(Sex == "female", na.rm = TRUE), Asthma = sum(Severity != "Healthy", na.rm = TRUE)) %>% pivot_longer(c("Female", "Asthma"), names_to = "Demographic", values_to = "sarp_value")
fisher_matrix <- merge(combined_summary, SARP_3_summary)
fisher_matrix <- fisher_matrix[, c("Cluster", "Demographic", "overall_total", "sarp_total", "overall_value", "sarp_value")]
fisher_p <- apply(fisher_matrix, 1, function(x) fisher.test(matrix(as.integer(x[3:6]), nrow = 2))$p.value)
fisher_matrix$p <- format(fisher_p, digits = 4)

fisher_matrix$overall_value <- paste0(fisher_matrix$overall_value, " (", format(100 * fisher_matrix$overall_value / fisher_matrix$overall_total, digits = 3), "%)")
fisher_matrix$sarp_value <- paste0(fisher_matrix$sarp_value, " (", format(100 * fisher_matrix$sarp_value / fisher_matrix$sarp_total, digits = 3), "%)")

female_summary <- fisher_matrix[fisher_matrix$Demographic == "Female", ]
colnames(female_summary) <- c("Cluster", "Demographic", "# Subjects (overall)", "# Subjects (SARP 3)", "# Female (overall)", "# Female (SARP 3)", "Fisher p-value (Female)")
female_summary <- subset(female_summary, select = -c(Demographic))
asthma_summary <- fisher_matrix[fisher_matrix$Demographic == "Asthma", ]
colnames(asthma_summary) <- c("Cluster", "Demographic", "# Subjects (overall)", "# Subjects (SARP 3)", "# Asthma (overall)", "# Asthma (SARP 3)", "Fisher p-value (Asthma)")
asthma_summary <- subset(asthma_summary, select = -c(Demographic))

cluster_demographic_summary <- merge(female_summary, asthma_summary)
write_csv(cluster_demographic_summary, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/table_S2_clustering_demographic_summary.csv")
###############################################################
