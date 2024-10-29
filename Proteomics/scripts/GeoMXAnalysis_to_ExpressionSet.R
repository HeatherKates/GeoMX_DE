# Load necessary libraries
library(readxl)
library(Biobase)

# Step 1: Read the sheets from the Excel file
file_path <- "../data/HuBMAP_nPOD_ProteomeAtlas_NegativeControlNorm_BackgroundCorrected.xlsx"

# Each donor has four blocks, one per region, one slide per donor with all four regions together. Two donors
# Had eight slices on one slide/scan. There is only one slice per block. 
# Control for slide and block (random) and donor (fixed) and ROI type (fixed).

# Subset P1-P4 only

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

# Visualize relationships between variables
# Load necessary libraries
library(dplyr)
library(ggplot2)

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

# Subset P1-P4
exprs_set <- exprs_set[,pData(exprs_set)$Donor %in% c("P1","P2","P3","P4")]
# Step 6: Convert it to a summarized experiment to do the normalization
# Step 1: Normalize and correct the expression data
norm_exprs <- exprs(exprs_set) * 
  pData(exprs_set)$NormalizationFactor * 
  pData(exprs_set)$`BackgroundCorrectionFactor (Protein)`

library(SummarizedExperiment)
library(S4Vectors)  # for DataFrame class

# Convert pData and featureData to DataFrame-compatible formats
col_data <- as.data.frame(pData(exprs_set))
col_data <- DataFrame(col_data)

row_data <- as.data.frame(featureData(exprs_set)@data)  # Extract data slot
row_data <- DataFrame(row_data)

# Create a list of assays
assays <- list(
  exprs = exprs(exprs_set),
  norm = norm_exprs
)

# Construct the SummarizedExperiment object
se <- SummarizedExperiment(
  assays = assays,
  colData = col_data,                  # Sample metadata
  rowData = row_data,                  # Feature metadata
  metadata = list(experimentData = experimentData(exprs_set))
)

# Differential Expression of anatomical region (pairwise comparisons) using mixed effects models and linear models (model testing)

# Load required packages
library(lme4)
library(lmerTest)
library(dplyr)
library(doParallel)
library(foreach)

# Set up parallel backend to use 6 CPUs
num_cores <- 6
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define model formulas
model_list <- list(
  ~ Anatomical_region + (1 | Donor),
  ~ Anatomical_region,
  ~ Anatomical_region + (1 | Donor) + (1 | Block_ID)
)

# Prepare list for pairwise comparisons
anatomical_region_pairs <- combn(unique(colData(se)$Anatomical_region), 2, simplify = FALSE)

# Initialize a list to store results for each model
final_results_list <- list()

# Loop over each model in model_list
for (model_formula in model_list) {
  model_name <- paste0("Model_", which(model_list == model_formula))
  model_results <- list()  # Store results for this model
  
  # Outer loop for segment types and anatomical region pairs
  for (segment_type in unique(colData(se)$Segment_name)) {
    segment_se <- se[, colData(se)$Segment_name == segment_type]
    
    for (pair in anatomical_region_pairs) {
      region_se <- segment_se[, colData(segment_se)$Anatomical_region %in% pair]
      
      # Define data matrices and sample data outside of parallel loop
      norm_exprs_data <- assay(region_se, "norm")
      original_exprs_data <- assay(region_se, "exprs")
      sample_data <- as.data.frame(colData(region_se))
      
      # Parallel loop for each protein/gene
      comparison_results <- foreach(i = 1:nrow(norm_exprs_data), .combine = rbind,
                                    .packages = c("lme4", "lmerTest")) %dopar% {
                                      
                                      # Construct the complete formula dynamically by adding the response
                                      full_formula <- as.formula(paste("norm_exprs_data[i,] ~", deparse(model_formula[[2]])))
                                      
                                      # Check if there is a random effect in the model formula
                                      has_random_effect <- grepl("\\|", deparse(model_formula[[2]]))
                                      
                                      # Use lmerTest::lmer if random effects are present, otherwise use lm
                                      if (has_random_effect) {
                                        model <- lmerTest::lmer(full_formula, data = sample_data)
                                        p_value <- summary(model)$coefficients[2, "Pr(>|t|)"]  # Extract p-value from lmerTest
                                        singular_message <- if (isSingular(model)) "Model fit is singular; results may be unreliable" else "Model fit is not singular"
                                      } else {
                                        model <- lm(full_formula, data = sample_data)
                                        p_value <- summary(model)$coefficients[2, "Pr(>|t|)"]  # Extract p-value from lm
                                        singular_message <- "No random effects"
                                      }
                                      
                                      # Calculate fold change and log2 fold change from original expression values
                                      mean_exprs <- tapply(original_exprs_data[i,], sample_data$Anatomical_region, mean)
                                      fold_change <- mean_exprs[pair[1]] / mean_exprs[pair[2]]
                                      log2_fold_change <- log2(fold_change)
                                      
                                      # Extract model fit statistics
                                      aic_value <- AIC(model)
                                      formula_str <- paste(deparse(formula(model)), collapse = " ")  # Updated to avoid line breaks
                                      
                                      
                                      # Store results for the current protein
                                      data.frame(
                                        protein = rownames(norm_exprs_data)[i],
                                        segment_type = segment_type,
                                        comparison = paste(pair, collapse = " vs "),
                                        p_value = p_value,
                                        fold_change = fold_change,
                                        log2_fold_change = log2_fold_change,
                                        singular_fit = singular_message,
                                        AIC = aic_value,
                                        model_formula = formula_str
                                      )
                                    }
      
      # Adjust p-values using Benjamini-Hochberg correction after loop
      comparison_results$adj_p_value <- p.adjust(comparison_results$p_value, method = "BH")
      
      # Add to results list with informative name
      result_name <- paste(segment_type, paste(pair, collapse = "_"), sep = ".")
      model_results[[result_name]] <- comparison_results
    }
  }
  
  # Combine results for this model
  final_results_list[[model_name]] <- bind_rows(model_results)
}

# Stop the cluster
stopCluster(cl)

# Combine all models' results for comparison
combined_results <- bind_rows(final_results_list, .id = "model")

# Find the best model for each row based on lowest AIC
model_testing_df <- combined_results %>%
  group_by(protein, segment_type, comparison) %>%
  slice_min(AIC, with_ties = FALSE) %>%
  select(protein, segment_type, comparison, model, AIC, p_value, adj_p_value, fold_change, log2_fold_change, model_formula, singular_fit)

# Calculate the range of adj_p_values and determine the best model based on AIC
model_testing_df <- combined_results %>%
  group_by(protein, segment_type, comparison) %>%
  mutate(
    adj_p_value_diff = paste0(min(adj_p_value), "-", max(adj_p_value))  # Calculate range of adj_p_values
  ) %>%
  slice_min(AIC, with_ties = FALSE) %>%
  select(protein, segment_type, comparison, model, AIC, p_value, adj_p_value, fold_change, log2_fold_change, model_formula, singular_fit, adj_p_value_diff)

# Add the additional columns based on the existing adj_p_value_diff column
model_testing_df <- model_testing_df %>%
  # Separate the min and max values from adj_p_value_diff
  mutate(
    min_adj_p_value = as.numeric(sub("-.*", "", adj_p_value_diff)),  # Extract the min p-value
    max_adj_p_value = as.numeric(sub(".*-", "", adj_p_value_diff)),  # Extract the max p-value
    # Calculate the difference between max and min adjusted p-values
    adj_p_value_range_diff = max_adj_p_value - min_adj_p_value,
    # Check if 0.05 falls within the range of adjusted p-values
    is_significant_difference = min_adj_p_value < 0.05 & max_adj_p_value > 0.05
  ) %>%
  # Drop the intermediate columns for clarity if not needed
  select(-min_adj_p_value, -max_adj_p_value)

# Histogram of the difference in adjusted p-value range
ggplot(model_testing_df, aes(x = adj_p_value_range_diff, fill = is_significant_difference)) +
  geom_histogram(binwidth = 0.01, position = "stack", color = "black", alpha = 0.6) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray"), name = "Crosses 0.05") +
  labs(
    title = "Histogram of Adjusted p-Value Differences Across Models",
    x = "Difference in Adjusted p-Value Range (max - min)",
    y = "Count"
  ) +
  theme_minimal()


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


