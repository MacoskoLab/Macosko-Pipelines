## Sample data

**RNA_data**:
* filtered_feature_bc_matrix.h5 or filtered_feature_bc_matrix/ (required)
* molecule_info.h5 (optional, used for QC plots)
* metrics_summary.csv (optional, used for QC plots)

**spatial_data**:
* FASTQs with spatial barcode data (required)
* Puck.csv file with coordinates for each spatial barcode (required)

## Sample outputs

**example_output**:  
* SBcounts.h5: spatial barcode count matrix produced from the spatial FASTQs
* cb_whitelist.txt: list of cell barcodes extracted from the RNA matrix
* matrix.csv.gz: spatial barcode count matrix filtered+matched to cb_whitelist.txt
* spatial_metadata.json: metadata recorded when converting spatial FASTQs into filtered matrix
* coords.csv: spatial coordinates for each cell in cb_whitelist.txt
* seurat.qs: Seurat object with both RNA data and spatial coordinates
* summary.pdf: a visual summary of the RNA and spatial data

## Pipeline commands

**commands**:
1. julia spatial-count.jl spatial_data spaial_data/Puck.csv
2. Rscript run-positioning.R 

Scripts and .wdl tasks can be found in Macosko-Pipelines/spatial-count and Macosko-Pipelines/positioning for steps 1 and 2, respectively