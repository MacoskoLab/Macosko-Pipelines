## [Sample data](https://drive.google.com/drive/folders/1BvupJwPw2le1KyIL0-4qyzA6yWOwXCjz?usp=drive_link)

**RNA_data**:
* filtered_feature_bc_matrix.h5 or filtered_feature_bc_matrix/ (required)
* molecule_info.h5 (optional, used for QC plots)
* metrics_summary.csv (optional, used for QC plots)

**spatial_fastqs**:
* FASTQs with spatial barcode data (required)

**pucks**:
* Puck.csv files with coordinates for each spatial barcode (required)

**example_output**:  
* SBcounts.h5: spatial barcode count matrix produced from the spatial FASTQs
* cb_whitelist.txt: list of cell barcodes extracted from the RNA matrix
* matrix.csv.gz: spatial barcode count matrix filtered+matched to cb_whitelist.txt
* spatial_metadata.json: metadata recorded during conversion from FASTQs into filtered matrix
* coords.csv: spatial coordinates for each cell in cb_whitelist.txt
* seurat.qs: Seurat object with both RNA data and spatial coordinates
* summary.pdf: a visual summary of the RNA and spatial data

## Pipeline commands

**commands**:
1. julia [spatial-count.jl](https://github.com/MacoskoLab/Macosko-Pipelines/tree/main/spatial-count) spatial_fastqs/ pucks/
2. Rscript [run-positioning.R](https://github.com/MacoskoLab/Macosko-Pipelines/tree/main/positioning) RNA_data/ SBcounts.h5
