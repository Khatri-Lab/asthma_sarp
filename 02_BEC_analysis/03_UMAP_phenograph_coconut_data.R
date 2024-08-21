###############################################################
library(umap)
library(bayesMetaintegrator)
library(Rphenoannoy)
# Load the coconut-combined expression
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_expr.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_pheno.rds")
# Load the processed meta objects for performing UMAP on raw data for supplemental figure
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/discovery_validation_bayes_meta_objs.rds")
###############################################################


###############################################################
# Convert probe names to gene names for the 359 DEG genes in coconut
list_of_common_genes <- rownames(combined_expr)
discovery_meta_obj <- probeToGene(discovery_meta_obj, list_of_common_genes)
validation_meta_obj <- probeToGene(validation_meta_obj, list_of_common_genes)
###############################################################


###############################################################
# Remove GSE18965 and merge the expression from others
discovery_meta_obj$bayesianMeta$originalData[[which(names(discovery_meta_obj$bayesianMeta$originalData) == "GSE18965")]] <- NULL
discovery_expr <- Reduce(merge, sapply(discovery_meta_obj$bayesianMeta$originalData, function(x) x$gene_expr))
validation_expr <- Reduce(merge, sapply(validation_meta_obj$bayesianMeta$originalData, function(x) x$gene_expr))
combined_raw_expr <- merge(discovery_expr, validation_expr)
rownames(combined_raw_expr) <- combined_raw_expr$gene
combined_raw_expr <- subset(combined_raw_expr, select = -c(gene, GSM2800880))
discovery_dataset <- Reduce(c, sapply(discovery_meta_obj$bayesianMeta$originalData, function(x) as.character(x$pheno$dataSet[rownames(x$pheno) != "GSM2800880"])))
validation_dataset <- Reduce(c, sapply(validation_meta_obj$bayesianMeta$originalData, function(x) as.character(x$pheno$dataSet[rownames(x$pheno) != "GSM2800880"])))
combined_raw_dataset <- c(discovery_dataset, validation_dataset)
###############################################################


###############################################################
# Compute UMAP for raw and conormalized data
compute_umap <- function(expr_matrix)
{
  umap_obj <- umap(t(expr_matrix), n_neighbors = 40, min_dist = 0.4)
  umap_data <- as.data.frame(umap_obj$layout)
  return(umap_data)
}
raw_umap <- compute_umap(combined_raw_expr)
raw_umap$dataset <- combined_raw_dataset
conorm_umap <- compute_umap(combined_expr)
###############################################################


###############################################################
# Compute phenograph clusters for conormalized data
phenograph_clusters <- Rphenoannoy::Rphenoannoy(data = t(combined_expr), k = 40)
cluster <- as.integer(membership(phenograph_clusters[[2]]))
combined_pheno$Cluster <- cluster
###############################################################


###############################################################
save(raw_umap, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/raw_data_umap_obj.rds")
save(conorm_umap, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conorm_data_umap_obj.rds")
save(combined_pheno, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_pheno_cluster.rds")
###############################################################
