###############################################################
library(tidyverse)

library(bayesMetaintegrator)
# Load the processed meta objects
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/discovery_validation_meta_objs.rds")
###############################################################


###############################################################
# Convert probe names to gene names for the genes that are in the discovery objects
list_of_common_genes <- rownames(discovery_meta_obj$metaAnalysis$pooledResults[discovery_meta_obj$metaAnalysis$pooledResults$numStudies > 1, ])
discovery_meta_obj <- probeToGene(discovery_meta_obj, list_of_common_genes)
validation_meta_obj <- probeToGene(validation_meta_obj, list_of_common_genes)
###############################################################


###############################################################
# Get dataset-wise effect sizes and filter out genes that did not converge
discovery_meta_obj <- getDatasetEffectSizes(discovery = discovery_meta_obj, cores = 8, steps = 2000)
validation_meta_obj <- getDatasetEffectSizes(discovery = validation_meta_obj, cores = 8, steps = 2000)
###############################################################


###############################################################
# Combine effect sizes: meta analysis
discovery_meta_obj <- combineDatasets(discovery = discovery_meta_obj)
validation_meta_obj <- combineDatasets(discovery = validation_meta_obj)
# Filter by effect size and bayesian equivalent of p-value for significantly differently expressed genes
sig_genes <- filterGenes(discovery_meta_obj, probabilityThresh = 0.05, effectSizeThresh = 0.5, numberStudiesThresh = 2)
sig_genes <- sig_genes$filterResults$bayesian_PR_0.05ES_0.5nStudies_2
sig_genes$gene_sig <- unlist(c(sig_genes$posGeneNames, sig_genes$negGeneNames))
save(discovery_meta_obj, validation_meta_obj, sig_genes, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/discovery_validation_bayes_meta_objs.rds")
###############################################################
