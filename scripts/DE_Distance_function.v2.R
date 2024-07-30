library(limma)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)
library(tidyverse)


library(standR)
data_path="/home/hkates/blue_garrett/Campbell-Thompson/P1-P4_DE_Analysis/data"
data_file="P1-P4.QC.v2.xlsx"

library(readxl)
# Get the sheet names
sheet_names <- excel_sheets(paste(data_path,data_file,sep="/"))

# Read each sheet and assign it to a data frame named after the sheet's name
for (sheet in sheet_names) {
  assign(sheet, read_excel(paste(data_path,data_file,sep="/"), sheet = sheet))
}

# Create the count data frame. Samples in columns and features/genes in rows. The first column is the gene names/ids
library(dplyr)
counts <- BioProbeCountMatrix %>% select(TargetName,all_of(SegmentProperties$SegmentDisplayName))

# Create the sample annotation data frame.
sampleAnno <- SegmentProperties

# Create the feature annotation data frame. 
featureAnno <- BioProbeCountMatrix %>% select(!all_of(SegmentProperties$SegmentDisplayName))

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
DistanceToMainPancreaticDuct <- read_csv("~/blue_garrett/Campbell-Thompson/GeoMX_DSP_atUF/DistanceToMainPancreaticDuct.csv")
DistanceToMainPancreaticDuct$merge <- paste(DistanceToMainPancreaticDuct$`Pancreas Section`,DistanceToMainPancreaticDuct$ROILabel,sep="_")
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

# Convert the new columns (last 8 items) to numeric
cols_to_convert <- tail(names(merged_data), 8)
merged_data[cols_to_convert] <- lapply(merged_data[cols_to_convert], as.numeric)

# Convert merged_data back to a list
merged_data_list <- as.list(merged_data)

# Update spe@colData@listData with the merged data
spe@colData@listData <- merged_data_list
names(spe@colData@listData)[names(spe@colData@listData) == "Distance from Main Duct (µm)"] <- "Distance_from_Main_Duct"
islet_ROIs <- colnames(spe)[colData(spe)[["ROI_type"]] == "islet_present"]
spe <- spe[, islet_ROIs]
#DE analysis function
library(limma)
library(edgeR)
library(voom)

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
results <- lapply(subgroups, perform_de_analysis)

# Optionally, you can name the elements of the results list
names(results) <- subgroups

# Combine individual DE plots into a single plot using cowplot
combined_DE_plot <- plot_grid(plotlist = lapply(results, function(res) res$DE_plot), ncol = 2)

# Add a title to the combined plot
combined_DE_plot <- plot_grid(
  ggdraw() + 
    draw_label("Differential Expression by Distance from Main Duct (µm)", size = 14),
  combined_DE_plot,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Print the combined DE plot
print(combined_DE_plot)

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
  
  combined_plot <- plot_grid(plotlist = plots, ncol = 4)
  
  combined_plot <- plot_grid(
    ggdraw() + 
      draw_label(paste("Gene Expression vs Distance from Main Duct in", unique(colData(spe_sub)$AOI_target), "enriched AOI (islet present ROI)\n (top 20 DE genes)"), size = 14),
    combined_plot,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )
  
  combined_plots[[subgroup]] <- combined_plot
}

# Print all combined plots
for (subgroup in subgroups) {
  print(combined_plots[[subgroup]])
}
# Combine plots for each subgroup
for (subgroup in subgroups) {
  de_results <- results[[subgroup]]$de_results
  spe_sub <- results[[subgroup]]$spe_sub
  top_genes <- rownames(de_results[order(de_results$P.Value), ])[1:20]
  plots <- lapply(top_genes, plot_gene_expression, spe_sub = spe_sub)
  
  combined_plot <- plot_grid(plotlist = plots, ncol = 4)
  
  combined_plot <- plot_grid(
    ggdraw() + 
      draw_label(paste("Gene Expression vs Distance from Main Duct in", unique(colData(spe_sub)$AOI_target), "enriched AOI (islet present ROI)\n (top 20 DE genes)"), size = 14),
    combined_plot,
    ncol = 1,
    rel_heights = c(0.1, 1)
  )
  
  combined_plots[[subgroup]] <- combined_plot
}

# Print all combined plots
for (subgroup in subgroups) {
  print(combined_plots[[subgroup]])
}