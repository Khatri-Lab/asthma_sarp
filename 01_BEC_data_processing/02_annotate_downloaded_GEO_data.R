###############################################################
# Annotate the following BEC datasets downloaded from NCBI-GEO to file:
# GSE4302, GSE18965, GSE44037, GSE41861, GSE64913, GSE89809, GSE104468, GSE67472.
###############################################################


###############################################################
library(MetaIntegrator)
# Load the datasets downloaded and annotated from GEO
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/downloaded_geo_data_unannotated.rds")
###############################################################


###############################################################
#GSE4302
# 42 nonsmoking subjects with mild/moderate asthma, 28 nonsmoking healthy controls, and 16 current smokers without asthma; bronchoscopy at baseline then 1 week after.
disease_status_orig_column <- as.character(GSE4302$pheno$`sample type:ch1`)
class <- rep(2, length(disease_status_orig_column))
class <- ifelse(grepl("Healthy control", disease_status_orig_column), 0, class)
class <- ifelse(grepl("Asthmatic at baseline", disease_status_orig_column), 1, class) #for now classifying all asthma as a case
class <- ifelse(grepl("Asthmatic after Placebo", disease_status_orig_column), 1, class)
class <- ifelse(grepl("Asthmatic after Flovent", disease_status_orig_column), 1, class)
# Smokers are removed
GSE4302$pheno$class <- class
GSE4302$class <- class
condition <- rep(NA, length(disease_status_orig_column))
condition <- ifelse(grepl("Healthy control", disease_status_orig_column), "healthy", condition)
condition <- ifelse(grepl("Asthmatic", disease_status_orig_column), "asthma", condition)
condition <- ifelse(grepl("Smoker", disease_status_orig_column), "smoker", condition)
GSE4302$pheno$condition <- condition
severity <- rep(NA, length(disease_status_orig_column)) #unlabeled "2's" are smokers
severity <- ifelse(grepl("Healthy control", disease_status_orig_column), "Healthy", severity)
severity <- ifelse(grepl("Asthmatic", disease_status_orig_column), "Mild/moderate", severity)
severity <- ifelse(grepl("Smoker", disease_status_orig_column), "Smoker", severity)
GSE4302$pheno$severity <- severity
GSE4302$pheno$tissue <- rep("bronchial epithelial cell", length(disease_status_orig_column)) 
GSE4302$pheno$allergic <- "Unknown"
GSE4302$pheno$age <- "24-61"
GSE4302$pheno$sex <- imputeSex(GSE4302)
###############################################################


###############################################################
#GSE18965
# Study of 112 healthy, healthy atopic, and atopic mild asthma patients; only 16 subjects Healthy (9) and Atopic asthma (7) had microarray done. No healthy atopics were included
# classified in "title" under HN = healthy normal or AA = atopic asthma
disease_status_orig_column <- as.character(GSE18965$pheno$title)
class <- rep(2, length(disease_status_orig_column))
class <- ifelse(grepl("HN", disease_status_orig_column), 0, class)
class <- ifelse(grepl("AA", disease_status_orig_column), 1, class)
GSE18965$pheno$class <- class
GSE18965$class <- class
#annotations for GSE18965
condition <- rep(2, length(disease_status_orig_column))
condition <- ifelse(grepl("HN", disease_status_orig_column), "healthy", condition)
condition <- ifelse(grepl("AA", disease_status_orig_column), "asthma", condition)
GSE18965$pheno$condition <- condition
GSE18965$pheno$severity <- ifelse(grepl("AA", disease_status_orig_column), "Mild", "Healthy")
GSE18965$pheno$tissue <- rep("bronchial epithelial cell", length(disease_status_orig_column))
GSE18965$pheno$allergic <- ifelse(grepl("AA", disease_status_orig_column), "allergic", "no")
GSE18965$pheno$age <- "2–17"
GSE18965$pheno$sex <- imputeSex(GSE18965)
###############################################################


###############################################################
#GSE41861
# upper and lower airway epithelial cells; serum IgE, periostin, and eosinophils as well; 54 asthmatic, 30 healthy. No associated paper. 
# $`severity:ch1` - Healthy, Mild, Moderate, Severe; $characteristics_ch1.1 Sex; $characteristics_ch1.4 - upper vs lower airway
# reported as allergic asthma for all asthmatics
disease_status_orig_column <- as.character(GSE41861$pheno$`severity:ch1`)
tissue_column <- as.character(GSE41861$pheno$characteristics_ch1.4)
class <- rep(2, length(disease_status_orig_column))
class <- ifelse(grepl("Healthy", disease_status_orig_column), 0, class)
class <- ifelse(grepl("Mild|Moderate|Severe", disease_status_orig_column), 1, class)
# Keep only the BEC samples and remove the Nasal samples
class <- ifelse(grepl("Nasal", tissue_column), 2, class)
GSE41861$pheno$class <- class
GSE41861$class <- class
GSE41861$pheno$condition <- ifelse(grepl("Mild|Moderate|Severe", disease_status_orig_column), "asthma", "healthy")
severity <- rep(NA, length(disease_status_orig_column))
severity <- ifelse("Healthy" == disease_status_orig_column, "Healthy", severity)
severity <- ifelse("Mild" == disease_status_orig_column, "Mild", severity)
severity <- ifelse("Moderate" == disease_status_orig_column, "Moderate", severity)
severity <- ifelse("Severe" == disease_status_orig_column, "Severe", severity)
GSE41861$pheno$severity <- severity
GSE41861$pheno$tissue <- rep(NA, length(tissue_column))
GSE41861$pheno$tissue <- ifelse(grepl("Nasal", tissue_column), "nasal epithelial cell", "bronchial epithelial cell")
GSE41861$pheno$age <- as.numeric(GSE41861$pheno$`age:ch1`)
GSE41861$pheno$allergic <- rep("allergic",length(disease_status_orig_column))
# All patients with asthma are described as allergic in the paper
GSE41861$pheno$allergic <- ifelse(grepl("asthma", GSE41861$pheno$condition), "allergic", "no")
GSE41861$pheno$sex <- GSE41861$pheno$sex_kl
###############################################################


###############################################################
#GSE44037
# 17 subjects; 6 allergic asthma, 5 allergic rhinitis, 6 healthy; nasal and BECs paired in subjects; studied off of medications
# asthma subjects were intermittent - interpretted as mild
disease_status_orig_column <- GSE44037$pheno$characteristics_ch1.4
tissue <- rep(2, length(disease_status_orig_column))
tissue <- ifelse(grepl("bronchial", GSE44037$pheno$`tissue:ch1`), "bronchial epithelial cell", tissue)
tissue <- ifelse(grepl("nasal", GSE44037$pheno$`tissue:ch1`), "nasal epithelial cell", tissue)
GSE44037$pheno$tissue <- tissue
class <- rep(2, length(disease_status_orig_column))
class <- ifelse(grepl("asthma", disease_status_orig_column), 1, class)
# Treating subjects with rhinitis as healthy controls
class <- ifelse(grepl("rhinitis", disease_status_orig_column), 0, class)
class <- ifelse(grepl("healthy", disease_status_orig_column), 0, class)
class <- ifelse(grepl("nasal epithelial cell", tissue), 2, class)
GSE44037$pheno$class <- class
GSE44037$class <- class
condition <- rep(2, length(disease_status_orig_column))
condition <- ifelse(grepl("asthma", disease_status_orig_column), "asthma", condition)
condition <- ifelse(grepl("healthy", disease_status_orig_column), "healthy", condition)
condition <- ifelse(grepl("rhinitis", disease_status_orig_column), "allergic rhinitis", condition)
GSE44037$pheno$condition <- condition
severity <- ifelse(grepl("asthma", disease_status_orig_column), "Mild", "Healthy")
severity <- ifelse(grepl("rhinitis", disease_status_orig_column), "Allergic rhinitis", severity)
GSE44037$pheno$severity <- severity
GSE44037$pheno$allergic <- ifelse(grepl("asthma|rhinitis", disease_status_orig_column), "allergic", "no")
GSE44037$pheno$age <- "20-30"
GSE44037$pheno$sex <- imputeSex(GSE44037)
###############################################################


###############################################################
#GSE64913
# 40 subjects, 17 severe asthma, 23 healthy control; attempted to compare central vs peripheral airway gene expression. Not all subjects had both samples available in the study. 
# Each patient had 4 samples taken from each location but they are represented and presumably analyzed together
disease_status_orig_column <- as.character(GSE64913$pheno$description)
class <- rep(2, length(disease_status_orig_column))
class <- ifelse(grepl("healthy", disease_status_orig_column), 0, class)
class <- ifelse(grepl("severe", disease_status_orig_column), 1, class)
GSE64913$pheno$class <- class
GSE64913$class <- class
GSE64913$pheno$condition <- ifelse(grepl("healthy",disease_status_orig_column), "healthy", "asthma")
severity <- rep(2, length(disease_status_orig_column))
severity <- ifelse(grepl("healthy", disease_status_orig_column), "Healthy", severity)
severity <- ifelse(grepl("severe", disease_status_orig_column), "Severe", severity)
GSE64913$pheno$severity <- severity
GSE64913$pheno$tissue <- rep("bronchial epithelial cell", length(disease_status_orig_column))
GSE64913$pheno$age <- as.numeric(GSE64913$pheno$`age:ch1`) 
GSE64913$pheno$allergic <- "Unknown"
GSE64913$pheno$sex <- imputeSex(GSE64913)
###############################################################


###############################################################
# GSE89809
# 46 patients with asthma, 19 healthy controls
disease_status_orig_column <- as.character(GSE89809$pheno$source_name_ch1)
class <- rep(2, length(disease_status_orig_column))
class <- ifelse(grepl("Epithelial brushings gene expression data from healthy control", disease_status_orig_column), 0, class)
class <- ifelse(grepl("Epithelial brushings gene expression data from mild asthmatic|Epithelial brushings gene expression data from moderate asthmatic|Epithelial brushings gene expression data from severe asthmatic", disease_status_orig_column), 1, class)
GSE89809$pheno$class <- class
GSE89809$class <- class
GSE89809$pheno$condition <- ifelse(grepl("asthmatic",GSE89809$pheno$source_name_ch1), "asthma", "healthy")
severity <- rep(2,length(disease_status_orig_column))
severity <- ifelse(grepl("heatlhy|healthy", disease_status_orig_column), "Healthy", severity)
severity <- ifelse(grepl("mild", disease_status_orig_column), "Mild", severity)
severity <- ifelse(grepl("moderate", disease_status_orig_column), "Moderate", severity)
severity <- ifelse(grepl("severe", disease_status_orig_column), "Severe", severity)
GSE89809$pheno$severity <- severity
tissue <- rep(NA,length(disease_status_orig_column))
tissue <- ifelse(grepl("Epithelial", disease_status_orig_column), "bronchial epithelial cell", tissue)
tissue <- ifelse(grepl("Bronchoalveolar", disease_status_orig_column), "bronchoalveolar lavage CD3+", tissue)
GSE89809$pheno$tissue <- tissue
GSE89809$pheno$age <- as.numeric(GSE89809$pheno$`age (yrs):ch1`)
GSE89809$pheno$allergic <- "Unknown"
GSE89809$pheno$sex <- imputeSex(GSE89809)
###############################################################



###############################################################
#GSE104468
# 12 allergic asthma, 12 controls; all non-smoking white; all on ICS, no nasal steroids; PBMC, nasal, and bronchial samples
disease_status_orig_column <- GSE104468$pheno$`disease state:ch1`
tissue <- rep(2, length(disease_status_orig_column))
tissue <- ifelse(grepl("PBMC", GSE104468$pheno$source_name_ch1), "PBMC", tissue)
tissue <- ifelse(grepl("Nasal", GSE104468$pheno$source_name_ch1), "nasal epithelial cell", tissue)
tissue <- ifelse(grepl("Bronchial", GSE104468$pheno$source_name_ch1), "bronchial epithelial cell", tissue)
GSE104468$pheno$tissue <- tissue
class <- rep(2, length(disease_status_orig_column))
class <- ifelse(grepl("Normal", disease_status_orig_column), 0,class)
class <- ifelse(grepl("Asthma", disease_status_orig_column), 1,class)
class <- ifelse(grepl("nasal epithelial cell|PBMC", tissue), 2, class)
GSE104468$pheno$class <- class
GSE104468$class <- class
condition <- rep("healthy", length(disease_status_orig_column))
condition <- ifelse(grepl("Asthma", disease_status_orig_column), "asthma", condition)
GSE104468$pheno$condition <- condition
severity <- rep("Healthy", length(disease_status_orig_column))
severity <- ifelse(grepl("Asthma", disease_status_orig_column), "Unspecified", severity)
GSE104468$pheno$allergic <- "Unknown"
GSE104468$pheno$severity <- severity
GSE104468$pheno$age <- as.numeric(GSE104468$pheno$`age:ch1`)
GSE104468$pheno$sex <- imputeSex(GSE104468)
###############################################################


###############################################################
#GSE67472
# 105 subjects, 62 asthma, all are Mild-moderate and were off treatment at the time of bronchcoscopy; some of these are sarp; select for the GSE67472 with "MASTA" or "MASTB" in "description"
disease_status_orig_column <- GSE67472$pheno$`disease state:ch1`
GSE67472$pheno$tissue <- rep("bronchial epithelial cell", length(disease_status_orig_column))
class <- ifelse(grepl("asthma", disease_status_orig_column), 1, 0)
# Remove samples that are already measured in SARP 1&2 data by keeping only those with "MASTA" or "MASTB" in "description"
class <- ifelse(grepl("MAST", GSE67472$pheno$description), class, 2)
class[rownames(GSE67472$pheno) == "GSM1647713"] <- 2
GSE67472$pheno$class <- class
GSE67472$class <- class
names(GSE67472$class) <- colnames(GSE67472$expr)
severity <- ifelse(grepl("asthma", disease_status_orig_column), "Mild/moderate", "Healthy")
GSE67472$pheno$severity <- severity
GSE67472$pheno$condition <- ifelse(grepl("asthma", disease_status_orig_column), "asthma", "healthy")
GSE67472$pheno$allergic <- "Unknown"
GSE67472$pheno$age <- as.numeric(GSE67472$pheno$`age:ch1`)
GSE67472$pheno$sex <- imputeSex(GSE67472)
###############################################################


###############################################################
save(GSE4302, GSE18965, GSE41861, GSE44037, GSE64913, GSE89809, GSE104468, GSE67472, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/downloaded_geo_data_annotated.rds")
###############################################################
