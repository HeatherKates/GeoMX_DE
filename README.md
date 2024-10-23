# Spatial RNAseq and Proteomics of HuBMAP Pancreas

**Overview**
This repository contains scripts and files related to the differential expression (DE) analysis of **RNAseq** and **protein expression data** generated from GeoMX DSP libraries from four donor pancreases.

## Repository Structure

```plaintext  
.
├── DE_by_region
│   ├── data
│   │   └── P1-P4.QC.v3.xlsx
│   ├── results
│   │   ├── By_Region_DE_Results_2024-10-23.xlsx
│   │   └── plots
│   │       ├── endothelial_Body_vs_Head_DE_plot_2024-10-23.png
│   │       ├── endothelial_Body_vs_Neck_DE_plot_2024-10-23.png
│   └── scripts
│       └── DE_Distance_function.R
├── DistanceToDuct
│   ├── data
│   │   └── DistanceToMainPancreaticDuct.csv
│   ├── results
│   │   └── CIBERSORT
│   └── scripts
│       └── CIBERSORT.R
├── HuBMAP_nPOD
│   └── data
├── Proteomics
│   └── data
└── README.md
```

### Directory Breakdown:

- **`DE_by_region/`**  
  Contains analyses of differential expression by pancreatic region (e.g., Body vs Head, Neck vs Tail).

  - **data/**: Input data used for the analyses.
    - **`P1-P4.QC.v3.xlsx`**: QC and data for four donors.
  
  - **results/**: Contains results and plots of the DE analysis.
    - **`By_Region_DE_Results_2024-10-23.xlsx`**: Output summarizing DE analysis.
    - **plots/**: Contains PNG plots for each analysis.

- **`DistanceToDuct/`**  
  Contains analyses related to **Does gene expression level change with increasing (or decreasing) distance fro
m main pancreatic duct?**

 Data was subset to regions of inerest that included an islet and each enriched-cell-type target (AOI) was analyzed separately (endothelial cells, acinar and other cells, beta cells, and duct cells)
  
  - **data/**: Input data for distance-based analyses.
    - **`DistanceToMainPancreaticDuct.csv`**: Data for pancreatic duct distances. Created by Sam Ewing.
  
  - **results/**: DE results based on distance.
    - **CIBERSORT/**: Contains CIBERSORT results and heatmaps.

- **`Proteomics/`**  
  Contains analyses of proteomics data.

  - **data/**: Contains raw and processed data.
  
- **`HuBMAP_nPOD/`**  
  Additional HuBMAP and nPOD analyses, including raw counts.

## Reproducibility

To reproduce the analysis:

1. **Download Original Data**: Contact for access to the original `.xlsx` files.
2. **Run Scripts**: Use the provided R scripts to process data and generate results.

## Contact

For questions, please contact:

Heather Kates  
Email: hkates@ufl.edu

## Acknowledgements

This research was supported by NIH NIDDK grant 1U54DK127823-01 Multi-omic 3D tissue maps for a Human BioMolecular Atlas.

