###############################################################
# Load SARP 1&2 and SARP 3 data into MetaIntegrator format
###############################################################


###############################################################
library(MetaIntegrator)
library(tidyverse)
###############################################################


###############################################################
#SARP 1&2
# Read data as given by SARP
SARP_1_2_expr <- read.csv("/labs/khatrilab/ilee/SARP/SARP Data/SARP 1-2/BEC_155subj_mRNA.csv")
SARP_1_2_pheno <- read.csv("/labs/khatrilab/ilee/SARP/SARP Data/SARP 1-2/BEC_155subj_pheno.csv")
##########
# Some data switched around and FVC_pre_liter and FVC_pre_per columns switched
SARP_1_2_pheno[SARP_1_2_pheno$SampleID == "bec0112F", "FEV1_ratio"] <- 0.72500
SARP_1_2_pheno[SARP_1_2_pheno$SampleID == "bec0112F", "FVC_pre_liter"] <- 79
SARP_1_2_pheno[SARP_1_2_pheno$SampleID == "bec1325H", "FVC_pre_per"] <- 3.21000
SARP_1_2_pheno[SARP_1_2_pheno$SampleID == "bec1325H", "FVC_pre_liter"] <- 104
names(SARP_1_2_pheno)[names(SARP_1_2_pheno)== "FVC_pre_liter"] <- "A"
names(SARP_1_2_pheno)[names(SARP_1_2_pheno)== "FVC_pre_per"] <- "B"
names(SARP_1_2_pheno)[names(SARP_1_2_pheno)== "A"] <- "FVC_pre_per"
names(SARP_1_2_pheno)[names(SARP_1_2_pheno)== "B"] <- "FVC_pre_liter"
##########
# Remove BEC from sample names
sample_name <- unlist(strsplit(as.character(SARP_1_2_pheno$SampleID), split = "bec"))
sample_name <- sample_name[!sample_name == ""]
SARP_1_2_pheno$sample_name <- sample_name
rownames(SARP_1_2_pheno) <- sample_name
##########
# Change colnames of expr data to remove an X and non-subject columns
rownames(SARP_1_2_expr) <- SARP_1_2_expr$GeneSymbol
SARP_1_2_expr <- SARP_1_2_expr[, !colnames(SARP_1_2_expr) %in% c("ProbeName", "Entrez", "GeneSymbol")]
# Remove X from column names
expr_names <- SARP_1_2_expr %>% colnames() %>% sub(pattern = "X", replacement = "")
colnames(SARP_1_2_expr) <- expr_names
# Remove 2042A from pheno because there isn't expr data for it
SARP_1_2_pheno <- SARP_1_2_pheno[SARP_1_2_pheno$sample_name %in% colnames(SARP_1_2_expr), ]
##########
# Add class, condition, and severity like the GEO datasets
SARP_1_2_pheno$class <- SARP_1_2_pheno$Asthma
SARP_1_2_pheno$condition <- "asthma"
SARP_1_2_pheno$condition[SARP_1_2_pheno$Asthma == 0] <- "healthy"
SARP_1_2_pheno$severity <- "Mild/moderate"
SARP_1_2_pheno$severity[SARP_1_2_pheno$Severe == 1] <- "Severe"
SARP_1_2_pheno$severity[SARP_1_2_pheno$Asthma == 0] <- "Healthy"
##########
# Make meta object and check
keys <- rownames(SARP_1_2_expr)
names(keys) <- rownames(SARP_1_2_expr)
SARP_1_2 <- list(expr = as.matrix(SARP_1_2_expr), pheno = SARP_1_2_pheno, class = SARP_1_2_pheno$class, keys = keys, formattedName = "SARP_1_2")
checkDataObject(SARP_1_2, "Dataset")
###############################################################


###############################################################
#SARP 3
# Read data as given by SARP
SARP_3_expr <- read.csv("/labs/khatrilab/ilee/SARP/SARP Data/SARP 3/RLENormalized_BatchCorrected_SARPEpithelialRNASeqData_Feb2021.csv")
SARP_3_pheno <- read.csv("/labs/khatrilab/ilee/SARP/SARP Data/SARP 3/SARP3_20220118l1_SARPtoV9a_v3_1.csv")
##########
# Set HGNC symbol as rownames and remove non-expr columns
rownames(SARP_3_expr) <- SARP_3_expr$gene_symbol
SARP_3_expr <- SARP_3_expr[, !colnames(SARP_3_expr) %in% c("X", "gene_symbol", "ensembl_gene_id")]
# Remove X and . from column names
expr_names <- SARP_3_expr %>% colnames() %>% sub(pattern = "X", replacement = "") %>% sub(pattern = "\\.", replacement = "") %>% sub(pattern = "\\.", replacement = "")
colnames(SARP_3_expr) <- as.character(expr_names)
##########
# Keep only those subjects with expression in the pheno
SARP_3_pheno <- SARP_3_pheno[SARP_3_pheno$part %in% colnames(SARP_3_expr), ]
rownames(SARP_3_pheno) <- as.character(SARP_3_pheno$part)
SARP_3_pheno <- SARP_3_pheno[colnames(SARP_3_expr), ]
# Add class, condition, and severity like the GEO datasets
SARP_3_pheno$class <- 1 - SARP_3_pheno$healthy_control
SARP_3_pheno$condition <- "asthma"
SARP_3_pheno$condition[SARP_3_pheno$healthy_control == 1] <- "healthy"
SARP_3_pheno$severity_number <- SARP_3_pheno$severity
SARP_3_pheno$severity <- "Mild/moderate"
SARP_3_pheno$severity[SARP_3_pheno$severity_number == 0] <- "Healthy"
SARP_3_pheno$severity[SARP_3_pheno$severity_number == 3] <- "Severe"
##########
# Make meta object and check
keys <- rownames(SARP_3_expr)
names(keys) <- rownames(SARP_3_expr)
SARP_3 <- list(expr = as.matrix(SARP_3_expr), pheno = SARP_3_pheno, class = SARP_3_pheno$class, keys = keys, formattedName = "SARP_3")
checkDataObject(SARP_3, "Dataset")
###############################################################


###############################################################
save(SARP_1_2, SARP_3, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/SARP_meta_objs.rds")
###############################################################
