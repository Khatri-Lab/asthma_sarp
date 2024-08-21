###############################################################
library(tidyverse)
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/discovery_validation_bayes_meta_objs.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_expr.rds")
###############################################################


###############################################################
effect_size_frame <- discovery_meta_obj$bayesianMeta$finalResults[discovery_meta_obj$bayesianMeta$finalResults$Gene %in% sig_genes$gene_sig, c("Gene", "ES", "ES_SD", "Tau", "TauSD", "Pr0", "nStudies")]
colnames(effect_size_frame) <- c("Gene", "Effect size", "Standard deviation (effect size)", "Tau (heterogeneity)", "Standard deviation (Tau)", "Probability (effect size = 0)", "# studies")
for(col_nums in c(2:6))
{
  effect_size_frame[, col_nums] <- format(effect_size_frame[, col_nums], digits = 3)
}
effect_size_frame$`Gene present after COCONUT conormalization` <- "No"
effect_size_frame$`Gene present after COCONUT conormalization`[effect_size_frame$Gene %in% rownames(combined_expr)] <- "Yes"
gene_order <- effect_size_frame$Gene[order(as.numeric(effect_size_frame$`Effect size`))]
effect_size_frame_for_printing <- effect_size_frame %>% arrange(match(Gene, gene_order))
write_csv(effect_size_frame_for_printing, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/table_S1_bayesian_effect_size_summary.csv")
###############################################################
