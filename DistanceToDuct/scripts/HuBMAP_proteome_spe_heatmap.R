#########################################
#####Read in data from GeoMX DSP#########
#########################################

library(standR)
data_file="data/HuBMAP_proteome_atlas_QCFiltered_Scaled_NCNorm_BGCorrected.xlsx"

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
  rmNegProbe = FALSE,
  colnames.as.rownames = c("TargetName", "SegmentDisplayName", "TargetName"),
  coord.colnames = c("ROICoordinateX", "ROICoordinateY")
)
#spe@colData$AOI_target <- gsub("endothelial_cells","endothelial",spe@colData$AOI_target)

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
Markers=NULL
#Marker genes
#Markers <- read.table("scripts/data/marker_genes.csv")
if(!is.null(Markers)){
colnames(Markers) <- c("CellTypes","Gene")
#Markers <- read.csv("scripts/data/Established_markers.txt",header=TRUE)
#colnames(Markers) <- c("Gene","CellTypes")
Markers <- Markers %>% filter(CellTypes %in% c("Acinar","Alpha","Beta","Ductal","Endothelial"))
Markers <- Markers %>%
  group_by(Gene) %>%
  filter(n() == 1) %>%
  ungroup()

# Subset the genes in Markers that are present in spe
valid_genes <- Markers$Gene[Markers$Gene %in% rownames(spe)]
Markers_filtered <- Markers[Markers$Gene %in% valid_genes, ]

# Subset the SpatialExperiment object for the valid genes
expr_matrix <- assay(spe, "logcounts")[valid_genes, ]
} else {
expr_matrix <- assay(spe,"logcounts")
}

#####################################
###AOI selection and annotation######
#####################################

Case_Subset=levels(as.factor(colData(spe)$Donor))
Cases="AllCases"

# Get the AOI_target annotation and Case annotation
col_annotation1 <- data.frame(Segment_type= colData(spe)$Segment_type)
col_annotation2 <- data.frame(Anatomical_Region=colData(spe)$Anatomical_region)
#col_annotation3 <- data.frame(Donor = colData(spe)$Donor)


# Combine the two annotations into one dataframe
col_annotation <- cbind(col_annotation1, col_annotation2)#,col_annotation3)
rownames(col_annotation) <- colnames(spe)

# Filtering
#col_annotation_filtered <- col_annotation[col_annotation$Donor %in% Case_Subset, , drop = FALSE]
col_annotation_filtered <- col_annotation
# Order by AOI_target first, then by Case
col_order <- order(col_annotation_filtered[[1]], col_annotation_filtered[[2]])#,col_annotation_filtered[[3]])

# Apply the ordering to both the annotations and the expression matrix
col_annotation_filtered <- col_annotation_filtered[col_order, , drop = FALSE]

# Subset the expression matrix to remove the corresponding columns
expr_matrix_filtered <- expr_matrix[, rownames(col_annotation_filtered)]
#expr_matrix_filtered <- expr_matrix_filtered[, col_order, drop = FALSE]

# Ensure rownames of expr_matrix match with Markers_filtered$Gene
if(!is.null(Markers)){
row_annotation <- data.frame(CellType = Markers_filtered$CellTypes)
rownames(row_annotation) <- Markers_filtered$Gene
}else{
  
}
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
#donor_colors <- c(brewer.pal(9, "Set1"), brewer.pal(3, "Set3"))  # Combining two color sets for 12 donors

# Assign colors to Donor levels
#names(donor_colors) <- levels(as.factor(colData(spe)$Donor))

# Define a custom color palette for CellType annotation
annotation_colors <- list(
  Segment_type = c(
    endocrine_cells     = "darkorange",
    exocrine_and_other   = "darkgreen",
    pancreatic_ducts = "red",
    vasculature  = "blue"
  ),
  Anatomical_Region=c(
    Body="black",
    Head="gray",
    Neck="purple",
    Tail="tan"
  )#,
  #Donor = donor_colors 
)
# Define your desired min and max for the color scale
color_min <- -7  # Adjust to desired minimum limit for color scale
color_max <- 7  # Adjust to desired maximum limit for color scale

# Create breaks for the color scale
breaks <- seq(color_min, color_max, length.out = length(my_palette) + 1)

# Plot heatmap with distinct row annotation colors
p <- pheatmap(expr_matrix_filtered,
         cluster_rows = TRUE,       # No clustering for rows
         cluster_cols = FALSE,       # No clustering for columns
         annotation_col = col_annotation_filtered,  # Add column annotations (AOI_target)
         annotation_colors = annotation_colors,  # Custom colors for annotations
         show_rownames = FALSE,     # Option to hide gene names
         show_colnames = FALSE,     # Option to hide sample names
         scale = "row",             # Scale expression values per gene
         color = my_palette,
         fontsize = 10,
         treeheight_row = 0,
         breaks=breaks,
         main = paste("Heatmap of GeoMx IO Proteome Atlas scaled probe counts across",ncol(expr_matrix_filtered), "AOIs\n from 12 donors by Segment Type and organ region"))        

ggsave(p, file=paste0("results/Proteome_Heatmap_noDonor.",Cases,".png"), width = 10, height = 8, dpi = 300)

