###############################################################
library(MetaIntegrator)
# Load the datasets downloaded and annotated from GEO
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/downloaded_geo_data_annotated.rds")
load("/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/SARP_meta_objs.rds")
###############################################################


###############################################################
# Format datasets to be compatible for meta-analysis by meta-integrator
# Removes samples not from cases or controls and fixes negative/NA expressions
clean_up_datasets_obtain_meta_obj <- function(list_of_datasets)
{
  cleaned_datasets <- list()
  for(x in 1:length(list_of_datasets))
  {
    dataset <- list_of_datasets[[x]]
    names(dataset$class) <- rownames(dataset$pheno)
    # Check if there are samples to be removed
    remove_index <- which(dataset$pheno$class == 2)
    if(length(remove_index) > 0)
    {
      dataset$pheno <- dataset$pheno[-remove_index, ]
      dataset$expr <- dataset$expr[, -remove_index]
      dataset$class <- dataset$class[-remove_index]
    }
    # If minimum expression is negative, add that to make expression matrix positive
    # Also update NAs to be the minimum expression value
    dataset_min <- min(dataset$expr, na.rm = TRUE)
    dataset$expr[is.na(dataset$expr)] <- dataset_min
    if(dataset_min < 0)
    {
      dataset$expr <- dataset$expr - dataset_min + 1
    }
    cleaned_datasets[[x]] <- dataset
    print(table(dataset$class))
  }
  # Set names for easier readability
  dataset_names <- sapply(cleaned_datasets, function(x) x$formattedName)
  names(cleaned_datasets) <- dataset_names
  # Check meta-object for validity
  meta_obj <- list()
  meta_obj$originalData <- cleaned_datasets
  meta_analysis_obj_check <- checkDataObject(meta_obj, "Meta", "Pre-Analysis")
  if(meta_analysis_obj_check)
  {
    meta_analysis_results <- runMetaAnalysis(meta_obj, maxCores = 11)
    return(meta_analysis_results)
  }
  else
  {
    return(meta_analysis_obj_check)
  }
}
###############################################################


###############################################################
discovery_datasets <- list(GSE4302, GSE18965, GSE41861, GSE44037, GSE64913, GSE89809, GSE104468)
discovery_meta_obj <- clean_up_datasets_obtain_meta_obj(discovery_datasets)

validation_datasets <- list(GSE67472, SARP_1_2, SARP_3)
validation_meta_obj <- clean_up_datasets_obtain_meta_obj(validation_datasets)

save(discovery_meta_obj, validation_meta_obj, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/discovery_validation_meta_objs.rds")
###############################################################