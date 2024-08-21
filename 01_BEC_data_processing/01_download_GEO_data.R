###############################################################
# Download and save the following BEC datasets from NCBI-GEO to file:
# GSE4302, GSE18965, GSE44037, GSE41861, GSE64913, GSE89809, GSE104468, GSE67472.
# SARP 1&2 and SARP 3 data were provided by csv and have been processesed in a separate file in this folder
###############################################################


###############################################################
library(MetaIntegrator)
###############################################################


###############################################################
# Download data from GEO
GEO_data <- getGEOData("GSE4302")
GSE4302 <- GEO_data$originalData$GSE4302
GEO_data <- getGEOData("GSE18965")
GSE18965 <- GEO_data$originalData$GSE18965
GEO_data <- getGEOData("GSE41861")
GSE41861 <- GEO_data$originalData$GSE41861
GEO_data <- getGEOData("GSE44037")
GSE44037 <- GEO_data$originalData$GSE44037
GEO_data <- getGEOData("GSE64913")
GSE64913 <- GEO_data$originalData$GSE64913
GEO_data <- getGEOData("GSE89809")
GSE89809 <- GEO_data$originalData$GSE89809
GEO_data <- getGEOData("GSE104468")
GSE104468 <- GEO_data$originalData$GSE104468
GEO_data <- getGEOData("GSE67472")
GSE67472 <- GEO_data$originalData$GSE67472
###############################################################


###############################################################
save(GSE4302, GSE18965, GSE41861, GSE44037, GSE64913, GSE89809, GSE104468, GSE67472, file = "/labs/khatrilab/ananthg/asthma_sarp/sarp_paper/data/downloaded_geo_data_unannotated.rds")
###############################################################
