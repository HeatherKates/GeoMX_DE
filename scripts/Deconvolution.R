data_path="data/"
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

# Calculate the total negative probe counts per sample
neg_probe_sums <- colSums(spe@metadata[["NegProbes"]])

# View the summary of the negative probe sums
summary(neg_probe_sums)
# Define a threshold (e.g., 90th percentile)
threshold <- quantile(neg_probe_sums, 0.90)

# Identify samples with low and high negative probe counts
low_neg_samples <- names(neg_probe_sums)[neg_probe_sums <= threshold]
high_neg_samples <- names(neg_probe_sums)[neg_probe_sums > threshold]

# Create the filtered matrix with low negative probe count samples
filtered_matrix <- spe[,colnames(spe) %in% low_neg_samples]
filtered_matrix <- filtered_matrix@assays@data$counts

# Create the raw matrix with all samples
raw_matrix <- spe@assays@data$counts
rownames(filtered_matrix) <- rownames(raw_matrix)
# Save the matrices
write.csv(as.data.frame(filtered_matrix), "scripts/data/filtered_matrix.csv")
write.csv(as.data.frame(raw_matrix), "scripts/data/raw_matrix.csv")

library(SoupX)
# Convert the matrices to a format compatible with SoupX (Matrix format)
library(Matrix)
raw_matrix <- as(as.matrix(raw_matrix), "dgCMatrix")
filtered_matrix <- as(as.matrix(filtered_matrix), "dgCMatrix")
library(Matrix)
rowSums = Matrix::rowSums
colSums = Matrix::colSums
# Create a SoupChannel object
soup_channel <- SoupChannel(raw_matrix, filtered_matrix)
soup_channel$tod <- filtered_matrix
spe@colData$AOI_target <- gsub("endothelial_cells","endothelial",spe@colData$AOI_target)
clusters <- as.factor(spe@colData$AOI_target)
clusters <- gsub("endothelial",0,clusters)
clusters <- gsub("duct_cells",1,clusters)
clusters <- gsub("acinar_and_other",2,clusters)
clusters <- gsub("beta_cells",3,clusters)
clusters <- gsub("immune_other",4,clusters)
clusters <- as.numeric(clusters)

# Set these clusters in the SoupChannel object
soup_channel = setClusters(soup_channel, setNames(clusters, rownames(spe@colData)))
soup_channel <- autoEstCont(soup_channel)

# Correct the raw expression data
corrected_matrix <- adjustCounts(soup_channel)

