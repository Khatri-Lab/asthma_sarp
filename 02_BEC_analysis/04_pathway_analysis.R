###############################################################
# Load the processed meta objects for finding effect sizes
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/discovery_validation_bayes_meta_objs.rds")
###############################################################


###############################################################
pathway_input_frame <- discovery_meta_obj$bayesianMeta$finalResults[discovery_meta_obj$bayesianMeta$finalResults$Gene %in% sig_genes$gene_sig, c("Gene", "ES")]
colnames(pathway_input_frame) <- c("gene", "ES")
pathway_input_frame$direction <- ifelse(pathway_input_frame$ES > 0, "up", "down")
write.csv(pathway_input_frame[, c("gene", "direction", "ES")], file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/pathway_analysis_input.csv", row.names = FALSE)
###############################################################


###############################################################
# pathway_analysis_input.csv goes into 
# https://rstudio-connect.khatrilab.stanford.edu/connect/#/apps/228/access
# and returns
# ArtZEnrichment.csv
###############################################################


###############################################################
pathway_output <- read.csv(file = paste0(out_folder, "ArtZEnrichment.csv"))
###############################################################