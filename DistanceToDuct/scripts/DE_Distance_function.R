library(standR)
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
spe@colData$AOI_target <- gsub("endothelial_cells","endothelial",spe@colData$AOI_target)
#Temporary test remove P4 and see if P1-P3 results are reproducible
spe <- spe[,spe@colData$Case %in% (c("P1","P2","P3"))]

#Variables
Subgroup_var="AOI_target"

#Metadata
library(tidyverse)
DistanceToMainPancreaticDuct <- read_csv("data/DistanceToMainPancreaticDuct.csv")
DistanceToMainPancreaticDuct$ROILabel <- sprintf("%03d", as.numeric(DistanceToMainPancreaticDuct$ROILabel))
DistanceToMainPancreaticDuct$merge <- paste(DistanceToMainPancreaticDuct$PancreasSection,DistanceToMainPancreaticDuct$ROILabel,sep="_")
# Use gsub to extract the first part before the first space
spe@colData@listData$Slide <- gsub(" .*", "", spe@colData@listData$SlideName)
spe@colData@listData$merge <- paste(spe@colData@listData$Slide,spe@colData@listData$ROILabel,sep="_")
spe@colData@listData$sample_id <- rownames(spe@colData)

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
rownames(spe@colData) <- merged_data$sample_id

library(SummarizedExperiment)
islet_ROIs <- colnames(spe)[colData(spe)[["ROI_type"]] == "islet_present"]
spe <- spe[, islet_ROIs]

#DE analysis function
library(limma)
library(edgeR)

# Define the subgroups
subgroups <- c("endothelial_cells", "beta_cells", "duct_cells", "acinar_and_other")

# Initialize lists to store results and combined plots
results <- list()
combined_plots <- list()

# Function to perform DE analysis and create plots
perform_de_analysis <- function(subgroup) {
  samples_subset <- colnames(spe)[colData(spe)[[Subgroup_var]] == subgroup]
  spe_sub <- spe[, samples_subset]
  spe_sub <- spe_sub[, !is.na(colData(spe_sub)$Distance_from_Main_Duct)]
  
  # Create the design matrix with the continuous variable and Case
  design <- model.matrix(~ Distance_from_Main_Duct + Case, data = colData(spe_sub))
  
  # Create DGEList object
  dge <- SE2DGEList(spe_sub)
  
  # Calculate normalization factors
  dge <- calcNormFactors(dge)
  
  # Apply voom transformation
  v <- voom(dge, design)
  
  # Fit linear model
  fit <- lmFit(v, design)
  efit <- eBayes(fit, robust = TRUE)
  
  # Extract differential expression results
  de_results <- topTable(efit, coef = 2, sort.by = "P", n = Inf)
  de_results$contrast <- "Distance from Main Duct"
  
  # Visualization
  de_results$color_factor <- factor(sign(de_results$logFC) * (de_results$P.Value < 0.05),
                                    levels = c(-1, 0, 1),
                                    labels = c("DOWN", "NOT DE", "UP"))
  
  DE_plot <- ggplot(de_results, aes(x = AveExpr, y = logFC, color = color_factor)) +
    geom_point(shape = 16, size = 1) +
    geom_text_repel(data = subset(de_results, P.Value < 0.05),
                    aes(label = ifelse(P.Value < 0.05, rownames(de_results), ""))) +
    scale_color_manual(values = c("DOWN" = "blue", "NOT DE" = "grey50", "UP" = "red")) +
    labs(title = subgroup, x = "Average Expression", y = "Regression Coefficient") +
    theme_minimal() +
    theme(legend.position = "none")
  
  results[[subgroup]] <- list(de_results = de_results, DE_plot = DE_plot, spe_sub = spe_sub)
}

# Perform DE analysis for each subgroup and save the results
library(ggplot2)
library(ggrepel)
library(cowplot)
results <- lapply(subgroups, perform_de_analysis)

# Optionally, you can name the elements of the results list
names(results) <- subgroups

#Write the DE results to file
library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# Iterate over each element in the results list
for (cell_type in names(results)) {
  # Get the de_results data frame
  de_results_df <- results[[cell_type]]$de_results
  
  # Add rownames as a new column named "Gene"
  de_results_df$Gene <- rownames(de_results_df)
  
  # Select and reorder columns
  de_results_df <- de_results_df[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "contrast", "color_factor")]
  
  # Rename "color_factor" to "change_direction"
  colnames(de_results_df)[which(names(de_results_df) == "color_factor")] <- "change_direction"
  colnames(de_results_df)[which(names(de_results_df) == "logFC")] <- "coefficient"
  # Add "model" and "subset" columns
  de_results_df$model <- "Gene Expression ~ Distance_from_Main_Duct + Case"
  de_results_df$subset <- "islet present ROIs"
  
  # Add a new sheet with the name of the cell type
  addWorksheet(wb, cell_type)
  
  # Write the data frame to the new sheet
  writeData(wb, sheet = cell_type, de_results_df)
}

# Save the workbook to a file
saveWorkbook(wb, file = paste0("results/DE_Results/Distance_from_duct_DE_Results.",Sys.Date(),".xlsx"), overwrite = TRUE)

# Define the path for saving plots
plot_save_path <- "results/plots/"

# Loop through each element in the results list
for (cell_type in names(results)) {
  # Get the DE_plot object
  de_plot <- results[[cell_type]]$DE_plot
  de_plot <- de_plot + ggtitle(paste("Differential Expression by Distance from Main Duct (µm) in", cell_type))
  # Save the modified plot as a PNG file
   ggsave(filename = paste0(plot_save_path, cell_type, "_DE_plot",Sys.Date(),".png"), plot = de_plot, width = 10, height = 8)
}

# Function to plot gene expression
plot_gene_expression <- function(gene_id, spe_sub) {
  plot_data <- data.frame(
    Distance_from_Main_Duct = colData(spe_sub)$Distance_from_Main_Duct,
    Case = colData(spe_sub)$Case,
    Gene_Expression = assay(spe_sub, "logcounts")[gene_id, ]
  )
  
  lm_model <- lm(Gene_Expression ~ Distance_from_Main_Duct + Case, data = plot_data)
  summary_lm <- summary(lm_model)
  slope <- coef(lm_model)[["Distance_from_Main_Duct"]]
  p_value <- summary_lm$coefficients[2, "Pr(>|t|)"]
  r_squared <- summary_lm$r.squared
  
  plot <- ggplot(plot_data, aes(x = Distance_from_Main_Duct, y = Gene_Expression, color = Case)) +
    geom_point(size = 1.5) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(
      title = gene_id,
      subtitle = paste("R^2 =", round(r_squared, 2), 
                       ", Slope =", round(slope, 5), 
                       ", p-value =", round(p_value, 4)),
      x = "Distance from Main Duct (µm)",
      y = "Gene Expression"
    ) +
    theme_minimal() +
    theme(
      plot.subtitle = element_text(size = 8),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8)
    )
  
  return(plot)
}

# Combine plots for each subgroup
for (subgroup in subgroups) {
  de_results <- results[[subgroup]]$de_results
  spe_sub <- results[[subgroup]]$spe_sub
  top_genes <- rownames(de_results[order(de_results$P.Value), ])[1:20]
  plots <- lapply(top_genes, plot_gene_expression, spe_sub = spe_sub)
  
  combined_plot <- plot_grid(plotlist = plots, ncol = 3)
  
  combined_plot <- plot_grid(
    ggdraw() + 
      draw_label(paste("Gene Expression vs Distance from Main Duct in", unique(colData(spe_sub)$AOI_target), "enriched AOI (islet present ROI)\n (top 20 DE genes)"), size = 14),
    combined_plot,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )
  
  combined_plots[[subgroup]] <- combined_plot
}

for (subgroup in subgroups) {
  # Get the DE_plot object
  plot <- combined_plots[[subgroup]]
  #plot <- plot + ggtitle(paste("Differential Expression by Distance from Main Duct (µm) in", cell_type))
  # Save the modified plot as a PNG file
  ggsave(filename = paste0(plot_save_path, subgroup, "_top_DE_genes_regression_plot",Sys.Date(),".png"), plot = plot, width = 10, height = 15)
}
