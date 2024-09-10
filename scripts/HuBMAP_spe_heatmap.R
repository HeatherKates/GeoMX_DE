#########################################
#####Read in data from GeoMX DSP#########
#########################################

library(standR)
data_file="data/HuBMAP_nPOD_Counts.xlsx"

library(readxl)
# Get the sheet names
sheet_names <- excel_sheets(data_file)

# Read each sheet and assign it to a data frame named after the sheet's name
for (sheet in sheet_names) {
  assign(sheet, read_excel(data_file, sheet = sheet))
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

#########################################
#####Normalize###########################
#########################################
spe_tmm <- geomxNorm(spe, method = "TMM")


#########################################
#####Heatmap############################
#########################################

# Load required libraries
library(pheatmap)
library(SpatialExperiment)
library(dplyr)
library(ggplot2)

#####################################
###Gene selection and annotation#####
#####################################

#Marker genes
Markers <- read.table("scripts/data/marker_genes.csv")
colnames(Markers) <- c("CellTypes","Gene")
#Markers <- read.csv("scripts/data/Established_markers.txt",header=TRUE)
#colnames(Markers) <- c("Gene","CellTypes")
Markers <- Markers %>% filter(CellTypes %in% c("Acinar","Alpha","Beta","Ductal","Endothelial"))
Markers <- Markers %>%
  group_by(Gene) %>%
  filter(n() == 1) %>%
  ungroup()
spe <- spe_tmm
# Subset the genes in Markers that are present in spe
valid_genes <- Markers$Gene[Markers$Gene %in% rownames(spe)]
Markers_filtered <- Markers[Markers$Gene %in% valid_genes, ]

#####################################
###AOI selection and annotation######
#####################################

Case_Subset=c("HuBMAP_P1", "HuBMAP_P2", "HuBMAP_P3", "HuBMAP_P4" ,"nPOD_6318" ,"nPOD_6488", "nPOD_6584")
Cases="AllCases"

# Subset the SpatialExperiment object for the valid genes
expr_matrix <- assay(spe, "logcounts")[valid_genes, ]

# Get the AOI_target annotation and Case annotation
col_annotation1 <- data.frame(AOI_target = colData(spe)$AOI_target)
col_annotation2 <- data.frame(Region=colData(spe)$organ_region)
col_annotation3 <- data.frame(Case = colData(spe)$Case)


# Combine the two annotations into one dataframe
col_annotation <- cbind(col_annotation1, col_annotation3,col_annotation2)
rownames(col_annotation) <- colnames(spe)

# Filter out any rows where AOI_target contains "immune"
col_annotation_filtered <- col_annotation[!grepl("NA", col_annotation$Region), , drop = FALSE]
col_annotation_filtered <- col_annotation_filtered[!grepl("immune", col_annotation_filtered$AOI_target), , drop = FALSE]
col_annotation_filtered <- col_annotation_filtered[col_annotation_filtered$Case %in% Case_Subset, , drop = FALSE]

# Order by AOI_target first, then by Case
col_order <- order(col_annotation_filtered[[1]], col_annotation_filtered[[2]],col_annotation_filtered[[3]])

# Apply the ordering to both the annotations and the expression matrix
col_annotation_filtered <- col_annotation_filtered[col_order, , drop = FALSE]

# Subset the expression matrix to remove the corresponding columns
expr_matrix_filtered <- expr_matrix[, rownames(col_annotation_filtered)]
#expr_matrix_filtered <- expr_matrix_filtered[, col_order, drop = FALSE]

# Ensure rownames of expr_matrix match with Markers_filtered$Gene
row_annotation <- data.frame(CellType = Markers_filtered$CellTypes)
rownames(row_annotation) <- Markers_filtered$Gene

# Sort the row annotation (markers) and expression matrix by CellType
#row_order <- order(row_annotation$CellType)
#row_annotation <- row_annotation[row_order, , drop = FALSE]
#expr_matrix_filtered <- expr_matrix_filtered[row_order, , drop = FALSE]


##########################
#####Plotting#############
##########################

# Define a custom color palette (blue to white to red)
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)
# Load RColorBrewer for diverging color palettes
library(RColorBrewer)

# Define a color palette with enough distinct colors for the Donor variable
donor_colors <- brewer.pal(7, "Set1")

# Assign colors to Donor levels
names(donor_colors) <- levels(as.factor(colData(spe)$Case))
# Define a custom color palette for CellType annotation
annotation_colors <- list(
  CellType = c(
    Acinar = "darkorange",
    Alpha = "lightgreen",
    Beta = "darkgreen",
    Ductal = "red",
    Endothelial = "blue"
  ),
  AOI_target = c(
    acinar_and_other = "darkorange",
    beta_cells = "darkgreen",
    duct_cells = "red",
    endothelial = "blue",
    `CD3+ T cells`="pink",
    `CD4+ T cells`   ="magenta",
    Other_Non_beta_Tcells="chartreuse"
    
  ),
  Region=c(
    Body="black",
    Head="gray",
    Neck="purple",
    Tail="tan"
  ),
  Case=donor_colors
)

# Plot heatmap with distinct row annotation colors
p <- pheatmap(expr_matrix_filtered,
         cluster_rows = FALSE,       # No clustering for rows
         cluster_cols = FALSE,       # No clustering for columns
         annotation_row = row_annotation,      # Add row annotations (cell types)
         annotation_col = col_annotation_filtered,  # Add column annotations (AOI_target)
         annotation_colors = annotation_colors,  # Custom colors for annotations
         show_rownames = FALSE,     # Option to hide gene names
         show_colnames = FALSE,     # Option to hide sample names
         scale = "row",             # Scale expression values per gene
         color = my_palette,
         fontsize = 10,
         main = paste("Heatmap of Elgamal et al. 2023 (HPAP) Cluster Marker Gene Expression\nacross",ncol(expr_matrix_filtered), "AOIs by AOI Target, HuBMAP and nPOD Case(",Cases,")", "and organ region"))        

ggsave(p, file=paste0("results/plots/HuBMAP_and_nPOD_CellTypeMarker_Heatmap.TmmNorm.",Cases,".png"), width = 10, height = 8, dpi = 300)

