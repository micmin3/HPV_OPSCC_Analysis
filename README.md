# HPV_OPSCC_Analysis
The scripts in this repository will allow you to perform the computational analyses from the paper "Cellular states are coupled to genomic and viral heterogeneity in HPV-related oropharyngeal carcinoma". Be aware that external files not included in the repository are needed for all analyses. How to access the relevant files is described at the top of each script. Description of the scripts below:

functions.R - Functions needed for all scripts, source this file before running anything.
preproc_celltype.R - Script for QC and assignment of cell type
cna.R - Defining malignant and nonmalignant epithelial cells and division of malignant cells into subclones
epithelial_analysis.R - Analysis focusing on TCGA subtypes, intertumour heterogeneity and nonmalignant epithelial cells
nmf.R - Defining metaprograms through NMF
hpv_analysis.R - Specific analyses to define the role of HPV in intratumour heterogeneity
ext_validation.R - Validation in external datasets

Tested in R 4.1.0
