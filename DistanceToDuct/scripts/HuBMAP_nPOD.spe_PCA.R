
library(standR)
library(SpatialExperiment)
library(readxl)
library(dplyr)
library(ggplot2)
data_file="data/HuBMAP_nPOD_Counts.xlsx"

# Get the sheet names
sheet_names <- excel_sheets(data_file)

# Read each sheet and assign it to a data frame named after the sheet's name
for (sheet in sheet_names) {
  assign(sheet, read_excel(data_file, sheet = sheet))
}

# Create the count data frame. Samples in columns and features/genes in rows. The first column is the gene names/ids
select <- dplyr::select
counts <- BioProbeCountMatrix %>% select(TargetName,all_of(SegmentProperties$SegmentDisplayName))

# Create the sample annotation data frame.
sampleAnno <- SegmentProperties

# Create the feature annotation data frame. 
featureAnno <- BioProbeCountMatrix %>% select(!all_of(SegmentProperties$SegmentDisplayName))
source("scripts/helper_functions/geomx_import_fun_HRK.R")
spe <- geomx_import_fun_HRK(
  counts,
  sampleAnno,
  featureAnno,
  rmNegProbe = TRUE,
  NegProbeName = "NegProbe-WTX",
  colnames.as.rownames = c("TargetName", "SegmentDisplayName", "TargetName"),
  coord.colnames = c("ROICoordinateX", "ROICoordinateY")
)
spe@colData$AOI_target <- gsub("endothelial_cells","endothelial",spe@colData$AOI_target)
spe <- spe[,colData(spe)$organ=="Pancreas"]
#########################################
#####Batch Correction####################
#########################################

drawPCA(spe, assay = 2, color = Case)
spe <- findNCGs(spe, batch_name = "Case", top_n = 300)

for(i in c(2:5)){
  spe_ruv <- geomxBatchCorrection(spe, factors = "AOI_target", 
                                  NCGs = metadata(spe)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 2, color = Case, title = paste0("k = ", i)))
  
}
spe_ruv <- geomxBatchCorrection(spe, factors = "AOI_target", 
                                NCGs = metadata(spe)$NCGs, k = 6)
spe_ruv <- scater::runPCA(spe_ruv)

pca_results_ruv <- reducedDim(spe_ruv, "PCA")
spe <- scater::runPCA(spe)
pca_results <- reducedDim(spe, "PCA")

p1 <- plotPairPCA(spe, precomputed = pca_results, color = AOI_target, title = "PCA without batch correction", n_dimension = 2)
p2 <- plotPairPCA(spe, precomputed = pca_results, color = Case, title = "PCA without batch correction", n_dimension = 2)
p3 <- plotPairPCA(spe_ruv, precomputed = pca_results_ruv, color = AOI_target, title = "PCA with batch correction (RUV4, k = 2)", n_dimension = 2)
p4 <- plotPairPCA(spe_ruv, precomputed = pca_results_ruv, color = Case, title = "PCA with batch correction (RUV4, k = 2)", n_dimension = 2)

library(patchwork)

# Combine the four ggplot objects into a 2x2 grid layout
combined_plot <-  (p1 | p2) / (p3 | p4)

# Save the combined plot to a file (e.g., PNG)
ggsave("results/plots/HuBMAP_nPOD_beforeandafterbatchCorrect.png", combined_plot, width = 14, height = 10, dpi = 300)
