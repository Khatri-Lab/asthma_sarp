###############################################################
library(umap)
library(Rphenoannoy)
# Load the processed coconut object
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/discovery_validation_bayes_meta_objs.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_with_all_genes_expr.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_with_all_genes_pheno.rds")
###############################################################


###############################################################
set.seed(42)
umap_obj <- umap(t(combined_expr), n_neighbors = 40, min_dist = 0.4)
umap_data <- as.data.frame(umap_obj$layout)
colnames(umap_data) <- c("UMAP_1", "UMAP_2")
set.seed(42)
phenograph_clusters <- Rphenoannoy::Rphenoannoy(data = t(combined_expr), k = 100)
cluster <- as.integer(membership(phenograph_clusters[[2]]))
combined_pheno$Cluster <- cluster
combined_pheno <- cbind(combined_pheno, umap_data)
save(combined_pheno, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_with_all_genes_pheno_with_umap_clusters.rds")
###############################################################


###############################################################
set.seed(1008)
umap_obj <- umap(t(combined_expr[!(rownames(combined_expr) %in% sig_genes$gene_sig), ]), n_neighbors = 40, min_dist = 0.4)
umap_data <- as.data.frame(umap_obj$layout)
colnames(umap_data) <- c("UMAP_1", "UMAP_2")
set.seed(1008)
phenograph_clusters <- Rphenoannoy::Rphenoannoy(data = t(combined_expr), k = 100)
cluster <- as.integer(membership(phenograph_clusters[[2]]))
combined_pheno$Cluster <- cluster
combined_pheno <- cbind(combined_pheno, umap_data)
save(combined_pheno, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_with_all_non_sig_genes_pheno_with_umap_clusters.rds")
###############################################################
