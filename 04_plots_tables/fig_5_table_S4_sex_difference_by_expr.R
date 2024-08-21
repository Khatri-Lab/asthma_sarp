###############################################################
library(tidyverse)
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/plot_themes.R")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_expr.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_pheno_cluster.rds")
###############################################################


###############################################################
combined_data <- cbind(t(combined_expr), combined_pheno)
combined_data <- combined_data[combined_data$Severity != "Healthy", ]
combined_data_long <- combined_data %>% pivot_longer(rownames(combined_expr), names_to = "genes", values_to = "expr")

wilcox_test_asthma_only <- combined_data_long %>% group_by(genes) %>% wilcox_test(expr ~ Sex)
eff_size_asthma_only <- combined_data_long %>% group_by(genes) %>% cohens_d(expr ~ Sex, hedges.correction = TRUE)
wilcox_test_asthma_only$fdr <- p.adjust(wilcox_test_asthma_only$p, method = "fdr")

sex_diff_genes <- wilcox_test_asthma_only$genes[wilcox_test_asthma_only$fdr < 0.1]
wilcox_test_asthma_only_cluster <- combined_data_long[combined_data_long$genes %in% sex_diff_genes, ] %>% group_by(Cluster, genes) %>% wilcox_test(expr ~ Sex)
eff_size_asthma_only_cluster <- combined_data_long[combined_data_long$genes %in% sex_diff_genes, ] %>% group_by(Cluster, genes) %>% cohens_d(expr ~ Sex, hedges.correction = TRUE)
wilcox_test_asthma_only_cluster$fdr <- p.adjust(wilcox_test_asthma_only_cluster$p, method = "fdr")
###############################################################


###############################################################
plot_gene_in_cluster <- function(gene_name, cluster_name)
{
  plot_data <- na.omit(combined_data_long[combined_data_long$genes == gene_name & combined_data_long$Cluster == cluster_name, c("expr", "Sex")])
  expr_plot <- ggplot(plot_data, aes(x = Sex, y = expr, fill = Sex)) +
    ylab(gene_name) + ggtitle(paste0("Cluster ", cluster_name)) +
    geom_boxplot(lwd = line_size_main_mm, width = 0.5, colour = black_text_colour, outlier.shape = NA) +
    geom_beeswarm(cex = cex_val, size = point_size, colour = black_text_colour) +
    scale_fill_manual(values = sex_colors) +
    theme(legend.position = "none", axis.title.x = element_blank()) +
    annotate("text", label = paste0("FDR = ", format(wilcox_test_asthma_only_cluster$fdr[wilcox_test_asthma_only_cluster$genes == gene_name & wilcox_test_asthma_only_cluster$Cluster == cluster_name], digits = 2)),
      x = 1.5, y = Inf,
      size = text_size_labels, colour = black_text_colour,
      vjust = 1.5) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  expr_plot
  return(expr_plot)
}
num_plots <- 0
list_of_plots <- list()
for(i in 1:nrow(wilcox_test_asthma_only_cluster))
{
  if(wilcox_test_asthma_only_cluster$fdr[i] < 0.05)
  {
    num_plots <- num_plots + 1
    list_of_plots[[num_plots]] <- plot_gene_in_cluster(wilcox_test_asthma_only_cluster$genes[i], wilcox_test_asthma_only_cluster$Cluster[i])
  }
}
combined_sex_plot <- plot_grid(plotlist = list_of_plots, nrow = 3, align = "hv", axis = "tblr")
ggsave("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/fig_5_sex_differences_by_cluster.pdf", combined_sex_plot, width = 7.2, height = 7.2, device = cairo_pdf)
###############################################################


###############################################################
combined_eff_size_fdr <- merge(eff_size_asthma_only[, c("genes", "effsize")], wilcox_test_asthma_only[, c("genes", "fdr")])
colnames(combined_eff_size_fdr) <- c("genes", "Hedges G (all clusters)", "FDR (all clusters)")
wilcox_test_asthma_only_cluster$Cluster <- paste0("FDR (cluster ", wilcox_test_asthma_only_cluster$Cluster, ")")
wilcox_test_asthma_only_cluster_wide <- wilcox_test_asthma_only_cluster[, c("genes", "fdr", "Cluster")] %>% pivot_wider(names_from = "Cluster", values_from = "fdr")
eff_size_asthma_only_cluster$Cluster <- paste0("Hedges G (cluster ", eff_size_asthma_only_cluster$Cluster, ")")
eff_size_asthma_only_cluster_wide <- eff_size_asthma_only_cluster[, c("genes", "effsize", "Cluster")] %>% pivot_wider(names_from = "Cluster", values_from = "effsize")
combined_eff_size_fdr <- merge(combined_eff_size_fdr, eff_size_asthma_only_cluster_wide, all.x = TRUE)
combined_eff_size_fdr <- merge(combined_eff_size_fdr, wilcox_test_asthma_only_cluster_wide, all.x = TRUE)
combined_eff_size_fdr <- combined_eff_size_fdr[, c(1:4, 8, 5, 9, 6, 10, 7, 11)]
for(i in 2:11)
{
  combined_eff_size_fdr[, i] <- format(combined_eff_size_fdr[, i], digits = 2)
}
write_csv(combined_eff_size_fdr, "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/table_S4_sex_differences_by_cluster.csv")
###############################################################
