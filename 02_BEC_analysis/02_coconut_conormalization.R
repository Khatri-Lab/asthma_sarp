###############################################################
library(MetaIntegrator)
# Load the processed meta objects
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/discovery_validation_bayes_meta_objs.rds")
###############################################################


###############################################################
# Remove GSE18965 and keep only the 505 significant DEGs for conormalization
keep_only_significant_genes <- function(dataset)
{
  dataset$expr <- dataset$expr[dataset$keys %in% sig_genes$gene_sig, ]
  dataset$keys <- dataset$keys[dataset$keys %in% sig_genes$gene_sig]
  return(dataset)
}
combined_meta_obj <- list(originalData = list(
  GSE41861 = keep_only_significant_genes(discovery_meta_obj$originalData$GSE41861),
  GSE64913 = keep_only_significant_genes(discovery_meta_obj$originalData$GSE64913),
  GSE89809 = keep_only_significant_genes(discovery_meta_obj$originalData$GSE89809),
  GSE4302 = keep_only_significant_genes(discovery_meta_obj$originalData$GSE4302),
  GSE104468 = keep_only_significant_genes(discovery_meta_obj$originalData$GSE104468),
  GSE44037 = keep_only_significant_genes(discovery_meta_obj$originalData$GSE44037),
  SARP.1.2 = keep_only_significant_genes(validation_meta_obj$originalData$SARP.1.2),
  SARP.3 = keep_only_significant_genes(validation_meta_obj$originalData$SARP.3),
  GSE67472 = keep_only_significant_genes(validation_meta_obj$originalData$GSE67472)
))
conormalized_results <- coconutMetaIntegrator(combined_meta_obj)
save(conormalized_results, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_meta_objs.rds")
###############################################################


###############################################################
# Combine co-normalized gene expression into a data frame
cases_expr <- Reduce(cbind, lapply(conormalized_results$COCONUTList, function(x) x$genes))
control_expr <- Reduce(cbind, lapply(conormalized_results$controlList$GSEs, function(x) x$genes))
combined_expr <- cbind(cases_expr, control_expr)
###############################################################


###############################################################
# Combine phenotype into a data frame
cases_pheno <- Reduce(rbind, lapply(conormalized_results$COCONUTList, function(x) data.frame(SampleID = rownames(x$pheno), Severity = as.character(x$pheno$severity), Dataset = as.character(x$pheno$dataSet), Sex = as.character(x$pheno$sex))))
control_pheno <- Reduce(rbind, lapply(conormalized_results$controlList$GSEs, function(x) data.frame(SampleID = rownames(x$pheno), Severity = as.character(x$pheno$severity), Dataset = as.character(x$pheno$dataSet), Sex = as.character(x$pheno$sex))))
combined_pheno <- rbind(cases_pheno, control_pheno)
rownames(combined_pheno) <- combined_pheno$SampleID
combined_pheno <- subset(combined_pheno, select = -c(SampleID))

combined_pheno$Severity <- ifelse(grepl("Mild|moderate|Moderate", combined_pheno$Severity), "Mild/moderate", combined_pheno$Severity)
combined_pheno$Severity <- ifelse(grepl("Severe", combined_pheno$Severity), "Severe", combined_pheno$Severity)
combined_pheno$Severity <- ifelse(grepl("healthy", combined_pheno$Severity), "Healthy", combined_pheno$Severity)
combined_pheno$Severity <- ifelse(grepl("Unspecified", combined_pheno$Severity), "Unspecified", combined_pheno$Severity)
combined_pheno$Severity <- ifelse(grepl("rhinitis", combined_pheno$Severity), "Healthy", combined_pheno$Severity)

combined_pheno$Dataset <- ifelse(grepl("SARP UCSF", combined_pheno$Dataset), "SARP 3", combined_pheno$Dataset)
###############################################################


###############################################################
# Remove one sample with extreme values- this affects subsequent UMAP and clustering
combined_expr <- subset(combined_expr, select = -c(GSM2800880))
combined_pheno <-  combined_pheno[rownames(combined_pheno) != "GSM2800880", ]
save(combined_expr, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_expr.rds")
save(combined_pheno, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_pheno.rds")
###############################################################
