###############################################################
library(ConsensusClusterPlus)
source("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/code/plot_themes.R")
# Load conormalized expression matrix
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_expr.rds")
###############################################################


###############################################################
consensus_clustering_results <- ConsensusClusterPlus(as.matrix(combined_expr),
  maxK = 6, reps = 200, pItem = 0.8, pFeature = 1,
  title = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/figures/",
  clusterAlg = "km", distance = "euclidean", seed = 42, plot = "pdf")
###############################################################
