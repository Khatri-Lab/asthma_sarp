###############################################################
library(MetaIntegrator)
# Load the processed meta objects
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/discovery_validation_bayes_meta_objs.rds")
###############################################################


###############################################################
# Keep all genes and conormalize- we still exclude GSE18965 to keep consistent comparison
combined_meta_obj <- list(originalData = list(
  GSE41861 = discovery_meta_obj$originalData$GSE41861,
  GSE64913 = discovery_meta_obj$originalData$GSE64913,
  GSE89809 = discovery_meta_obj$originalData$GSE89809,
  GSE4302 = discovery_meta_obj$originalData$GSE4302,
  GSE104468 = discovery_meta_obj$originalData$GSE104468,
  GSE44037 = discovery_meta_obj$originalData$GSE44037,
  SARP.1.2 = validation_meta_obj$originalData$SARP.1.2,
  SARP.3 = validation_meta_obj$originalData$SARP.3,
  GSE67472 = validation_meta_obj$originalData$GSE67472))
conormalized_results <- coconutMetaIntegrator(combined_meta_obj)
save(conormalized_results, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_meta_objs_all_genes.rds")
###############################################################
