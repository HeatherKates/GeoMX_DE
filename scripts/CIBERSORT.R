library(standR)
data_path="/home/hkates/blue_garrett/Campbell-Thompson/P1-P4_DE_Analysis/data"
data_file="P1-P4.QC.v3.xlsx"

library(readxl)
# Get the sheet names
sheet_names <- excel_sheets(paste(data_path,data_file,sep="/"))

# Read each sheet and assign it to a data frame named after the sheet's name
for (sheet in sheet_names) {
  assign(sheet, read_excel(paste(data_path,data_file,sep="/"), sheet = sheet))
}

# Create the count data frame. Samples in columns and features/genes in rows. The first column is the gene names/ids
library(dplyr)
select <- dplyr::select
counts <- BioProbeCountMatrix %>% select(TargetName,all_of(SegmentProperties$SegmentDisplayName))

# Create the sample annotation data frame.
sampleAnno <- SegmentProperties

# Create the feature annotation data frame. 
featureAnno <- BioProbeCountMatrix %>% select(!all_of(SegmentProperties$SegmentDisplayName))
source("/home/hkates/blue_garrett/Campbell-Thompson/P1-P4_DE_Analysis/scripts/geomx_import_fun_HRK.R")
spe <- geomx_import_fun_HRK(
  counts,
  sampleAnno,
  featureAnno,
  rmNegProbe = TRUE,
  NegProbeName = "NegProbe-WTX",
  colnames.as.rownames = c("TargetName", "SegmentDisplayName", "TargetName"),
  coord.colnames = c("ROICoordinateX", "ROICoordinateY")
)

#Variables
Subgroup_var="AOI_target"

#compare to previously used spe for debugging
#spe_compare <- readRDS("/blue/timgarrett/hkates/Campbell-Thompson/GeoMX_DSP_atUF/R/PostQC_P1-P3_spe.Rds")

#Metadata
library(tidyverse)
DistanceToMainPancreaticDuct <- read_csv("/home/hkates/blue_garrett/Campbell-Thompson/P1-P4_DE_Analysis/data/DistanceToMainPancreaticDuct.csv")
DistanceToMainPancreaticDuct$ROILabel <- sprintf("%03d", as.numeric(DistanceToMainPancreaticDuct$ROILabel))
DistanceToMainPancreaticDuct$merge <- paste(DistanceToMainPancreaticDuct$PancreasSection,DistanceToMainPancreaticDuct$ROILabel,sep="_")
# Use gsub to extract the first part before the first space
spe@colData@listData$Slide <- gsub(" .*", "", spe@colData@listData$SlideName)
spe@colData@listData$merge <- paste(spe@colData@listData$Slide,spe@colData@listData$ROILabel,sep="_")

# Merge the data frames
merged_data <- merge(
  as.data.frame(spe@colData@listData),
  DistanceToMainPancreaticDuct,
  by="merge",
  all.x = TRUE
)

# Remove the temporary segment_mapped column
merged_data$merge<- NULL

#convert literal "N/A" to NA
merged_data[merged_data == "N/A"] <- NA

# Convert the new columns (last 8 items) to numeric
cols_to_convert <- tail(names(merged_data), 8)
merged_data[cols_to_convert] <- lapply(merged_data[cols_to_convert], as.numeric)

# Convert merged_data back to a list
merged_data_list <- as.list(merged_data)

# Update spe@colData@listData with the merged data
spe@colData@listData <- merged_data_list
names(spe@colData@listData)[names(spe@colData@listData) == "Distance from Main Duct (µm)"] <- "Distance_from_Main_Duct"

library(SummarizedExperiment)
library(CIBERSORT)
source("cibersort_mod.R") #modify not to quit if a sample fails
#preprocessCore MUST be installed manually to disable threading or cibersort will fail. Use steps below
#git clone https://github.com/bmbolstad/preprocessCore.git
#cd preprocessCore
#R CMD INSTALL --configure-args="--disable-threading" -l /home/hkates/R/x86_64-pc-linux-gnu-library/4.4/ .
library(preprocessCore)
library(pheatmap)

# Define the function
cibersort_func <- function(subset) {
  
  # Subset the data based on the AOI_target
  subset_ROIs <- colnames(spe)[colData(spe)[["AOI_target"]] == subset]
  spe_subset <- spe[, subset_ROIs]
  
  # Extract and process the logcounts data
  logcounts <- as.data.frame(spe_subset@assays@data@listData[["logcounts"]])
  logcounts$`Gene symbol` <- rownames(logcounts)
  logcounts <- logcounts %>% dplyr::relocate("Gene symbol")
  write_tsv(logcounts, file = paste0("/blue/timgarrett/hkates/Campbell-Thompson/P1-P4_DE_Analysis/scripts/data/",subset, "ROI.txt"))
  
  # Run CIBERSORT and clean the results
  results <- cibersort(sig_matrix, paste0("/blue/timgarrett/hkates/Campbell-Thompson/P1-P4_DE_Analysis/scripts/data/",subset, "ROI.txt"))
  results_clean <- results[, 1:(ncol(results) - 3)]
  
  # Replace NA values with a placeholder value, e.g., 0
  df_clean <- results_clean
  df_clean[is.na(df_clean)] <- 0
  
  # Create a data frame for the row annotations
  annotation_row <- data.frame(ROI = spe_subset@colData$ROI_type)
  rownames(annotation_row) <- rownames(df_clean)  # Ensure rownames match
  
  # Create the heatmap with row annotations
  heatmap_plot <- pheatmap(
    as.matrix(df_clean),
    treeheight_row = 0,
    treeheight_col = 0,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    annotation_row = annotation_row,  # Add the row annotation
    main = paste("Heatmap of CIBERSORT results:", subset, "AOIs")
  )
  
  # Return a list with the cleaned results and the heatmap plot
  return(list(results_clean = df_clean, heatmap_plot = heatmap_plot))
}

# Example of running the function for each subset
results_acinar <- cibersort_func("acinar_and_other")
results_beta <- cibersort_func("beta_cells")
results_duct <- cibersort_func("duct_cells")
results_endothelial <- cibersort_func("endothelial_cells")

#Print output files
# List of results for easier iteration
results_list <- list(
  acinar_and_other = results_acinar,
  beta_cells = results_beta,
  duct_cells = results_duct,
  endothelial_cells = results_endothelial
)

# Save each heatmap to a PNG file
for (subset_name in names(results_list)) {
  heatmap_plot <- results_list[[subset_name]]$heatmap_plot
  
  # Set the filename using the subset name
  filename <- paste0( "/blue/timgarrett/hkates/Campbell-Thompson/P1-P4_DE_Analysis/results/CIBERSORT/","Heatmap_of_CIBERSORT_results_", subset_name, "AOIs.png")
  
  # Save the plot to a PNG file
  png(filename, width = 3000, height = 2500, res = 300)
  print(heatmap_plot)  # Print the plot to the device
  dev.off()  # Close the device
}
