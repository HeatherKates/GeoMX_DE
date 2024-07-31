# Spatial RNAseq of HuBMAP pancreas

This repository contains scripts and files related to the DE analysis of RNAseq data generated from GeoMX DSP libraries. 

## Overview

The code provided here is fully reproducible. Except for the original dataset .xlsx files (generated on the GeoMX DSP analysis suite based on proprietary counts and QC DCC files generated from .fastq.gx files by Nanostring NGS pipeline on BaseSpace by ICBR staff), all other files and results can be recreated using the scripts included in this repository. Due to size constraints and data ownership, large files and non public files are omitted. For access to these files, please see the `.gitignore` and email [hkates@ufl.edu](mailto:hkates@ufl.edu).

## Directory Structure

- **scripts/**: Contains all the scripts necessary to reproduce the analysis.
- **results/**: Directory for storing analysis results.

## Analysis

** Does gene expression level change with increasing (or decreasing) distance from main pancreatic duct? Data was subset to regions of inerest that included an islet and each enriched-cell-type target (AOI) was analyzed separately (endothelial cells, acinar and other cells, beta cells, and duct cells)


## Reproducibility

To ensure full reproducibility, follow the steps below:

1. **Download Original Data**: Obtain the original .xlsx files. Authorized collaborators email: [hkates@ufl.edu](mailto:hkates@ufl.edu) or access file of same name on dropbox.
2. **Run Scripts**: Use the scripts provided in the `scripts` directory to process the data and generate results.

## Contact

For access to data by authorized project collaborators or any other inquiries, please contact:

Heather Kates  
Email: [hkates@ufl.edu](mailto:hkates@ufl.edu)

## License

None

## Acknowledgements

This research was supported by NIH NIDDK grant 1U54DK127823-01 Multi-omic 3D tissue maps for a Human BioMolecular Atlas
