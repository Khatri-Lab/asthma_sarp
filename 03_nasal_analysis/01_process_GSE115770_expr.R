###############################################################
library(data.table)
library(readr)
library(tidyverse)
library('biomaRt')
Sys.setenv(BIOMART_CACHE = "/labs/khatrilab/ananthg/")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/conormalized_expr.rds")
###############################################################


###############################################################
# Data link: https://github.com/BenaroyaResearch/Peds_Asthma_Modules/tree/master/data
GSE115770_expr <- fread("/labs/khatrilab/ananthg/asthma_sarp/GSE115770_github_data_ananlyses/github_repo/data/raw_counts_nasal_muppits523.txt")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gene_names <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"), values = GSE115770_expr$V1, mart = mart)

colnames(GSE115770_expr)[1] <- colnames(gene_names)[1]
GSE115770_expr <- merge(GSE115770_expr, gene_names, by = colnames(gene_names)[1])
GSE115770_expr <- GSE115770_expr[, -1]

gene_counts <- table(GSE115770_expr$hgnc_symbol)
multiple_entry_matrix <- GSE115770_expr[GSE115770_expr$hgnc_symbol %in% names(gene_counts[gene_counts > 1]), ]
single_entry_matrix <- GSE115770_expr[GSE115770_expr$hgnc_symbol %in% names(gene_counts[gene_counts == 1]), ]
pruned_expr_matrix <- multiple_entry_matrix %>%
  pivot_longer(colnames(multiple_entry_matrix)[1:ncol(multiple_entry_matrix) - 1], names_to = "share_id", values_to = "expr") %>%
  group_by(hgnc_symbol, share_id) %>%
  summarize(median_expr = median(expr)) %>%
  pivot_wider(names_from = share_id, values_from = median_expr)
expr_matrix <- rbind(single_entry_matrix, pruned_expr_matrix[, colnames(single_entry_matrix)])

colnames(expr_matrix) <- sapply(strsplit(colnames(expr_matrix), split = "_"), function(x) x[[1]])
colnames(expr_matrix)[length(colnames(expr_matrix))] <- "hgnc_symbol"
###############################################################


###############################################################
GSE115770_pheno <- read_csv("/labs/khatrilab/ananthg/asthma_sarp/GSE115770_github_data_ananlyses/github_repo/data/totalNasalDesign_update_clinical variable groups.csv")
GSE115770_pheno <- GSE115770_pheno[, -1]
colnames(GSE115770_pheno) <- GSE115770_pheno[4, ]
GSE115770_pheno <- GSE115770_pheno[-c(1:4), ]
###############################################################


###############################################################
asthma_genes <- rownames(combined_expr)
baseline_subjects <- GSE115770_pheno$library.sampleId[GSE115770_pheno$Analysis.Visit == "Visit 0"]
baseline_pheno <- GSE115770_pheno[GSE115770_pheno$Analysis.Visit == "Visit 0", ]

baseline_expr <- as.data.frame(expr_matrix[expr_matrix$hgnc_symbol %in% asthma_genes, ])
rownames(baseline_expr) <- baseline_expr$hgnc_symbol
baseline_expr <- t(baseline_expr[, colnames(baseline_expr) %in% baseline_subjects])

save(baseline_expr, baseline_pheno, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/GSE115770_baseline_expr_pheno.rds")
write.csv(data.frame(baseline_expr), file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/GSE115770_baseline_expr.csv")
write_csv(baseline_pheno, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/GSE115770_baseline_pheno.csv")
###############################################################


###############################################################
sarp_1_2 <- read_csv("/labs/khatrilab/ilee/SARP/SARP Data/SARP 1-2/BEC_155subj_mRNA.csv")
sarp_3 <- read_csv("/labs/khatrilab/ilee/SARP/SARP Data/SARP 3/RawCountData_SARPEpithelialRNASeqData_Feb2021.csv")

SARP_cluster_pheno <- data.frame(t(read.csv("/labs/khatrilab/ilee/SARP/AvH_Paper_2022/Analysis - including GSE67472/Figures_7_26_2022/figure_3/coconut_expr.csv", row.names = 1)))
modified_subj_id <- SARP_cluster_pheno %>% rownames() %>% sub(pattern = "X", replacement = "") %>% sub(pattern = "\\.", replacement = "-") %>% sub(pattern = "\\.", replacement = "-")

sarp_1_2 <- as.data.frame(sarp_1_2[sarp_1_2$GeneSymbol %in% asthma_genes, ])
rownames(sarp_1_2) <- sarp_1_2$GeneSymbol
sarp_1_2 <- sarp_1_2[, -c(1, 2, 3)]
sarp_1_2 <- as.data.frame(t(sarp_1_2))
sarp_1_2$dataset <- "sarp_1_2"

sarp_3 <- as.data.frame(sarp_3[sarp_3$gene_symbol %in% asthma_genes, ])
rownames(sarp_3) <- sarp_3$gene_symbol
sarp_3 <- sarp_3[, -c(1, 2)]
sarp_3 <- as.data.frame(t(sarp_3))
sarp_3$dataset <- "sarp_3"

sarp_data <- rbind(sarp_1_2[, c(asthma_genes, "dataset")], sarp_3[, c(asthma_genes, "dataset")])
sarp_data$subject_id <- rownames(sarp_data)

sarp_data <- merge(sarp_data, data.frame(subject_id = modified_subj_id, Cluster = SARP_cluster_pheno$cluster), by = "subject_id")
write_csv(sarp_data, "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/SARP_cluster_pheno_expr.csv")
###############################################################
