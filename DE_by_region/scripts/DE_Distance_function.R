library(standR)
data_path="../data/"
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
source("helpers.R")
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
#Temporary test remove P4 and see if P1-P3 results are reproducible
#spe <- spe[,spe@colData$Case %in% (c("P1","P2","P3"))]

#Variables
Subgroup_var="AOI_target"

#Metadata
library(tidyverse)
# Use gsub to extract the first part before the first space
spe@colData@listData$Slide <- gsub(" .*", "", spe@colData@listData$SlideName)
spe@colData@listData$merge <- paste(spe@colData@listData$Slide,spe@colData@listData$ROILabel,sep="_")
spe@colData@listData$sample_id <- rownames(spe@colData)

library(SummarizedExperiment)

#DE analysis function
library(limma)
library(edgeR)

# Define the subgroups
subgroups <- c("endothelial", "beta_cells", "duct_cells", "acinar_and_other")

# Initialize lists to store results and combined plots
results <- list()
combined_plots <- list()

library(lme4)
library(lmerTest)  # For p-values in mixed models
library(edgeR)

perform_de_analysis_edgeR <- function(subgroup) {
  samples_subset <- colnames(spe)[colData(spe)[[Subgroup_var]] == subgroup]
  spe_sub <- spe[, samples_subset]
  spe_sub <- spe_sub[, !is.na(colData(spe_sub)$organ_region)]
  spe_sub <- spe_sub[, !colData(spe_sub)$organ_region=="NA"]
  
  # Create DGEList object from the counts and sample metadata
  dge <- SE2DGEList(spe_sub)
  
  # Calculate normalization factors
  dge <- calcNormFactors(dge)
  
  # Create the design matrix with organ_region and Case as fixed effects
  colData(spe_sub)$organ_region <- as.factor(colData(spe_sub)$organ_region)
  colData(spe_sub)$Case <- as.factor(colData(spe_sub)$Case)  # Ensure Case is a factor
  design <- model.matrix(~ 0 + organ_region + Case, data = colData(spe_sub))
  
  # Estimate dispersion
  dge <- estimateDisp(dge, design)
  
  # Fit the GLM model
  fit <- glmFit(dge, design)
  
  # Get levels of organ_region
  organ_region_levels <- levels(colData(spe_sub)$organ_region)
  
  # Dynamically create all pairwise contrasts
  contrast_pairs <- combn(organ_region_levels, 2, simplify = FALSE)
  
  results <- list()
  
  # Iterate through each contrast and perform DE analysis
  for (pair in contrast_pairs) {
    level1 <- pair[1]
    level2 <- pair[2]
    contrast_name <- paste(level1, "vs", level2)
    
    # Create a contrast vector for this comparison
    contrast_vector <- rep(0, ncol(design))
    contrast_vector[which(colnames(design) == paste0("organ_region", level1))] <- 1
    contrast_vector[which(colnames(design) == paste0("organ_region", level2))] <- -1
    
    # Perform the contrast
    lrt <- glmLRT(fit, contrast = contrast_vector)
    
    # Extract top differentially expressed genes
    de_results <- topTags(lrt, n = Inf)$table
    de_results$contrast <- contrast_name
    
    # Visualization
    de_results$color_factor <- factor(sign(de_results$logFC) * (de_results$FDR < 0.05),
                                      levels = c(-1, 0, 1),
                                      labels = c("DOWN", "NOT DE", "UP"))
    
    DE_plot <- ggplot(de_results, aes(x = logCPM, y = logFC, color = color_factor)) +
      geom_point(shape = 16, size = 1) +
      geom_text_repel(data = subset(de_results, FDR < 0.05),
                      aes(label = ifelse(FDR < 0.05, rownames(de_results), ""))) +
      scale_color_manual(values = c("DOWN" = "blue", "NOT DE" = "grey50", "UP" = "red")) +
      labs(title = contrast_name, x = "Log Counts per Million (logCPM)", y = "Log Fold Change (logFC)") +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Store the results and the plot
    results[[contrast_name]] <- list(de_results = de_results, DE_plot = DE_plot)
  }
  
  return(results)
}


# Perform DE analysis for each subgroup and save the results
library(ggplot2)
library(ggrepel)
library(cowplot)
results <- lapply(subgroups, perform_de_analysis_edgeR)

# Optionally, you can name the elements of the results list
names(results) <- subgroups

library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Iterate over each cell type in the results list
for (cell_type in names(results)) {
  # Iterate over each comparison within the current cell type
  for (comparison in names(results[[cell_type]])) {
    # Get the de_results data frame for the current comparison
    de_results_df <- results[[cell_type]][[comparison]]$de_results
    
    # Check if rownames are present, and add them as a new column "Gene" if needed
    if (!"Gene" %in% colnames(de_results_df)) {
      de_results_df$Gene <- rownames(de_results_df)
    }
    
    # Add the "isSignificant" column: TRUE if FDR < 0.05, otherwise FALSE
    de_results_df$isSignificant <- de_results_df$FDR < 0.05
    
    # Rename "color_factor" to "change_direction" if it exists
    if ("color_factor" %in% colnames(de_results_df)) {
      colnames(de_results_df)[which(names(de_results_df) == "color_factor")] <- "change_direction"
    }
    
    # Rename "logFC" to "coefficient" if it exists
    if ("logFC" %in% colnames(de_results_df)) {
      colnames(de_results_df)[which(names(de_results_df) == "logFC")] <- "coefficient"
    }
    
    # Add "model" and "subset" columns
    de_results_df$model <- "Gene Expression ~ organ region + Case"
    de_results_df$subset <- cell_type
    
    # Create a new sheet for each comparison within the current cell type
    sheet_name <- paste(cell_type, comparison, sep = "_")
    sheet_name <- gsub("[^A-Za-z0-9_]", "_", sheet_name)  # Clean up the sheet name if necessary
    
    # Add a new sheet with the cleaned-up name
    addWorksheet(wb, sheet_name)
    
    # Write the entire data frame to the new sheet
    writeData(wb, sheet = sheet_name, de_results_df)
  }
}

# Save the workbook to a file
output_file <- paste0("../results/By_Region_DE_Results_", Sys.Date(), ".xlsx")
saveWorkbook(wb, file = output_file, overwrite = TRUE)

cat("Workbook saved to", output_file, "\n")


library(gridExtra)

# Define the path for saving plots
plot_save_path <- "../results/plots/"

# Loop through each cell type in the results list
for (cell_type in names(results)) {
  # Initialize an empty list to store the plots for the current cell type
  plot_list <- list()
  
  # Loop through each comparison within the current cell type
  for (comparison in names(results[[cell_type]])) {
    # Get the DE_plot object for the current comparison
    de_plot <- results[[cell_type]][[comparison]]$DE_plot
    
    # Add a title to the plot based on the cell type and comparison
    de_plot <- de_plot + ggtitle(paste("Differential Expression by organ region in", cell_type, "-", comparison))
    
    # Store the plot in the plot_list
    plot_list[[comparison]] <- de_plot
  }
  
  # Arrange the plots in a 2x3 grid and save them as a single PNG file
  combined_plot <- marrangeGrob(grobs = plot_list, nrow = 2, ncol = 3)
  plot_file_name <- paste0(plot_save_path, cell_type, "_combined_DE_plots", ".png")
  
  # Save the combined plot
  ggsave(filename = plot_file_name, combined_plot, width = 15, height = 10)
}


