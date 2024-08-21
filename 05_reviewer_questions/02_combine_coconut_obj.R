###############################################################
library(umap)
library(Rphenoannoy)
# Load the processed coconut object
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_meta_objs_all_genes.rds")
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
save(combined_expr, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_with_all_genes_expr.rds")
save(combined_pheno, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_with_all_genes_pheno.rds")
###############################################################
