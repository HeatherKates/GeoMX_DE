# Load necessary libraries
library(readxl)
library(Biobase)

# Step 1: Read the sheets from the Excel file
file_path <- "../data/HuBMAP_nPOD_ProteomeAtlas_NegativeControlNorm_BackgroundCorrected.xlsx"

segment_properties <- read_excel(file_path, sheet = "SegmentProperties")
target_properties <- read_excel(file_path, sheet = "TargetProperties")
target_count_matrix <- read_excel(file_path, sheet = "TargetCountMatrix")
bio_probe_properties <- read_excel(file_path, sheet = "BioProbeProperties")

#Step 1b: Remove special chars
# Sanitize the SegmentDisplayName by replacing special characters
segment_properties$SegmentDisplayName <- gsub("/", "_", segment_properties$SegmentDisplayName)
segment_properties$SegmentDisplayName <- gsub("\\|", "_", segment_properties$SegmentDisplayName)
segment_properties$SegmentDisplayName <- gsub(" ", "_", segment_properties$SegmentDisplayName)

# Set the sanitized SegmentDisplayName as the rownames for sample annotations
rownames(segment_properties) <- segment_properties$SegmentDisplayName

# Step 2: Process the data

# Extract the expression matrix
exprs_data <- as.matrix(target_count_matrix[, -1])  # Remove first column (TargetName)
# Update the column names in the expression matrix (exprs_data) accordingly
colnames(exprs_data) <- segment_properties$SegmentDisplayName
rownames(exprs_data) <- target_count_matrix$TargetName  # Set rownames as TargetName

# Step 3: Create the phenoData (sample annotations)
sample_annotations <- segment_properties
rownames(sample_annotations) <- sample_annotations$SegmentDisplayName  # Set rownames for sample annotations

# Step 4: Create the featureData (target/probe annotations)
feature_annotations <- target_properties
rownames(feature_annotations) <- feature_annotations$TargetName  # Set rownames for feature annotations

# Step 5: Create the ExpressionSet object
exprs_set <- ExpressionSet(
  assayData = exprs_data,
  phenoData = AnnotatedDataFrame(data = sample_annotations),
  featureData = AnnotatedDataFrame(data = feature_annotations)
)

# Step 6: Verify the ExpressionSet object
exprs_set

# Do a lil test with it

# Load required packages
library(lme4)
library(pheatmap)
library(Biobase)
library(dplyr)

# Step 1: Subset ExpressionSet for segments with Segment_name == "beta_cells"
endothelial_cells_exprs <- exprs_set[, pData(exprs_set)$Segment_name == "endothelial_cells"]

# Extract expression data and phenoData (sample annotations)
exprs_data_sub <- exprs(endothelial_cells_exprs)
sample_data_sub <- pData(endothelial_cells_exprs)

# Step 2: Test for DE between ROI_type "islet_present" and "islet_absent"
# Initialize the results dataframe to store p-values
results <- data.frame(protein = rownames(exprs_data_sub), p_value = NA)

# Loop through each protein/probe and fit a linear mixed model
for (i in 1:nrow(exprs_data_sub)) {
  # Fit the linear mixed model using lmerTest for p-values
  model <- lmer(exprs_data[i, ] ~ ROI_type + (1 | Donor) + (1 | Anatomical_region), data = sample_data_sub)
  
  # Extract the summary from the model
  model_summary <- summary(model)
  
  # Store the p-value for ROI_type (2nd row of fixed effects, under ROI_type)
  results$p_value[i] <- model_summary$coefficients["ROI_typeislet_present", "Pr(>|t|)"]
}

# Adjust p-values using the Benjamini-Hochberg method (FDR)
results$adj.p.val <- p.adjust(results$p_value, method = "BH")

# Step 3: Print the top 20 DE proteins based on p-value
top_20_proteins <- results %>%
  arrange(adj.p.val) %>%
  head(20)

print(top_20_proteins)

# Step 4: Plot a heatmap of the top 20 DE proteins

# Z-score normalization of expression data (standardize rows)
top_20_exprs_z <- t(apply(top_20_exprs, 1, scale))  # Z-score normalize each protein

# Ensure row names and column names are preserved after scaling
rownames(top_20_exprs_z) <- rownames(top_20_exprs)
colnames(top_20_exprs_z) <- colnames(top_20_exprs)

# Create an annotation for ROI_type to display at the top of the heatmap
annotation <- data.frame(ROI_type = sample_data$ROI_type, Region= sample_data$Anatomical_region,
                         Scan = sample_data$ScanLabel, Donor = sample_data$Donor)
rownames(annotation) <- colnames(top_20_exprs_z)

# Define custom colors for annotations
annotation_colors <- list(
  ROI_type = c(islet_present = "darkgreen", islet_absent = "orange")  # Red for present, blue for absent
)

# Plot the heatmap with custom annotation colors
pheatmap(top_20_exprs_z, 
         annotation_col = annotation, 
         annotation_colors = annotation_colors,  # Add custom colors here
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),  # Blue-white-red color scale
         scale = "none")  # Scaling already done manually


#PCA

# Load ggplot2
library(ggplot2)

# Step 1: Standardize the expression data (rows are proteins, columns are samples)
exprs_data <- exprs(exprs_set)
# Z-score standardization by protein (row-wise scaling)
exprs_data_z <- t(apply(exprs_data, 1, scale))

# Extract expression data and phenoData (sample annotations)
sample_data <- pData(exprs_set)

# Step 2: Perform PCA on the standardized data
pca_result <- prcomp(t(exprs_data_z))  # Transpose to get samples in rows, proteins in columns

# Step 3: Create a data frame with PCA results and segment annotations
pca_df <- as.data.frame(pca_result$x)  # Extract the PCA scores
pca_df$Segment_name <- sample_data$Segment_name  # Add segment name information
pca_df$Donor <- sample_data$Donor

# Step 4: Plot the first two principal components with ggplot2 and color by segment name
ggplot(pca_df, aes(x = PC1, y = PC2, color = Segment_name)) +
  geom_point(size = 3) +
  labs(title = "PCA of Expression Data",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal() +
  scale_color_discrete(name = "Segment Name")


# Step 4: Plot the first two principal components with ggplot2 and color by Donor
ggplot(pca_df, aes(x = PC1, y = PC2, color = Donor)) +
  geom_point(size = 3) +
  labs(title = "PCA of Expression Data",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal() +
  scale_color_discrete(name = "Donor")


