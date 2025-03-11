library(glue) ; g=glue ; len=length
library(magrittr)
library(purrr)
library(tidyr)
library(dplyr)
library(data.table)
library(Seurat)
library(stringr)
library(ggplot2)
library(cowplot)
library(viridisLite)
library(viridis)
library(qpdf)
library(qs)

# Helper files required for loading spatial data and positioning cells
stopifnot(file.exists("dbscan.R", "helpers.R"))
source("helpers.R")
# setDTthreads(parallelly::availableCores())

# Load arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  RNApath <- args[[1]]
  SBpath <- args[[2]]
  out_path <- "."
} else if (length(args) == 3) {
  RNApath <- args[[1]]
  SBpath <- args[[2]]
  out_path <- args[[3]]
} else {
  stop("Usage: Rscript run-positioning.R RNApath SBpath [output_path]", call. = FALSE)
}

# Check input arguments
stopifnot("RNA dir not found" = dir.exists(RNApath))
stopifnot("RNA dir empty" = len(list.files(RNApath)) > 0)
if (dir.exists(SBpath)) { SBpath = file.path(SBpath, "SBcounts.h5") }
stopifnot("SB path is not SBcounts.h5" = str_sub(SBpath, -11, -1) == "SBcounts.h5")
stopifnot("SBcounts.h5 not found" = file.exists(SBpath))
if (!dir.exists(out_path)) { dir.create(out_path, recursive = T) }
stopifnot("Could not create output path" = dir.exists(out_path))

# Print arguments
print(g("RNA dir: {normalizePath(RNApath)}"))
print(g("RNA files: {list.files(RNApath) %>% paste0(collapse=', ')}"))
print(g("SB file: {normalizePath(SBpath)}"))
print(g("Output dir: {normalizePath(out_path)}"))

# Determine the RNA method
if (dir.exists(file.path(RNApath, "filtered_feature_bc_matrix"))) {
  rna_source = "10X"
  matrix_path = file.path(RNApath, "filtered_feature_bc_matrix")
  intronic_path = file.path(RNApath, "molecule_info.h5")
  metrics_path = file.path(RNApath, "metrics_summary.csv")
} else if (file.exists(file.path(RNApath, "filtered_feature_bc_matrix.h5"))) {
  rna_source = "10X_h5"
  matrix_path = file.path(RNApath, "filtered_feature_bc_matrix.h5")
  intronic_path = file.path(RNApath, "molecule_info.h5")
  metrics_path = file.path(RNApath, "metrics_summary.csv")
} else if (any(str_ends(list.files(RNApath), "_filtered_mtx_files.tar"))) {
  rna_source = "Optimus"
  matrix_path <-  list.files(RNApath, full.names=T) %>% keep(~str_sub(.,-5) == ".h5ad") %T>% {stopifnot(len(.)==1)}
  intronic_path <- list.files(RNApath, full.names=T) %>% keep(~str_sub(.,-5) == ".h5ad") %T>% {stopifnot(len(.)<=1)}
  metrics_path <- list.files(RNApath, full.names=T) %>% keep(~str_sub(.,-20) == "_library_metrics.csv") %T>% {stopifnot(len(.)<=1)}
} else {
  stop("Unknown RNA technique", call. = FALSE)
}
print(g("RNA source: {rna_source}\n"))

### Load the RNA ###############################################################

# Load the matrix
if (rna_source == "10X") {
  mat <- Read10X(matrix_path)
} else if (rna_source == "10X_h5") {
  mat <- Read10X_h5(matrix_path)
} else if (rna_source == "Optimus") {
  mat <- ReadOptimus(matrix_path)
}
obj <- CreateSeuratObject(mat, project = "Slide-tags")
rm(mat) ; gc()

# Add metadata
obj[["logumi"]] <- log10(obj$nCount_RNA+1)
obj[["cb"]] <- map_chr(colnames(obj), ~sub("-[0-9]*$", "", .)) %T>% {stopifnot(!duplicated(.))}
obj[["cb_index"]] <- 1:ncol(obj)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(MT-|mt-)")

# Compute %intronic
obj$pct.intronic = ReadIntronic(intronic_path, unname(obj$cb))

# Normalize, HVG, Scale, PCA, Neighbors, Clusters, UMAP
obj %<>%
  Seurat::NormalizeData(verbose=F) %>%
  Seurat::FindVariableFeatures(verbose=F) %>%
  Seurat::ScaleData(verbose=F) %>%
  Seurat::RunPCA(npcs=50, verbose=F) %>%
  Seurat::FindNeighbors(dims=1:30, verbose=F) %>%
  Seurat::FindClusters(resolution=0.8, verbose=F) %>%
  Seurat::RunUMAP(dims=1:30, umap.method="uwot", metric="cosine", verbose=F)

# Load the RNA metadata
metrics_df <- ReadLibraryMetrics(metrics_path)
Misc(obj, "RNA_metadata") <- setNames(metrics_df[,2], metrics_df[,1])

print(g("Loaded {ncol(obj)} cells"))

# Page 1
plot_metrics_csv(metrics_df) %>% make.pdf(file.path(out_path, "RNAmetrics.pdf"), 7, 8)

# Page 2
plot_cellcalling(intronic_path, obj$cb) %>% make.pdf(file.path(out_path,"RNAcalls.pdf"), 7, 8)

# Page 3
plot_umaps(obj) %>% make.pdf(file.path(out_path,"RNAumap.pdf"), 7, 8)

### Load the Spatial ###########################################################

# Load the spatial barcode counts matrix and fuzzy match to our called-cells whitelist
print("Generating the matrix...")
cb_whitelist = unname(obj$cb)
writeLines(cb_whitelist, file.path(out_path, "cb_whitelist.txt"))

system(g("Rscript load_matrix.R {SBpath} {file.path(out_path, 'cb_whitelist.txt')} {out_path}"))
stopifnot(file.exists(file.path(out_path, "matrix.csv.gz"), file.path(out_path, "spatial_metadata.json")))
Misc(obj, "spatial_metadata") <- jsonlite::fromJSON(file.path(out_path, "spatial_metadata.json"))

### DBSCAN #####################################################################

# Assign a position to each whitelist cell
print("Positioning cells...")
system(g("Rscript positioning.R {file.path(out_path, 'matrix.csv.gz')} {out_path}"))
stopifnot(file.exists(file.path(out_path,"coords.csv")))
coords <- read.table(file.path(out_path,"coords.csv"), header=T, sep=",")
coords %<>% right_join(data.frame(cb_index=1:len(cb_whitelist)), by = "cb_index") %>% arrange(cb_index) # fill in cells with no spatial data
Misc(obj, "coords") <- coords

# Create spatial reduction
stopifnot(nrow(coords) == ncol(obj), coords$cb_index == obj$cb_index)
obj$x_um <- coords$x_um
obj$y_um <- coords$y_um
# Add DBSCAN coords
emb = coords %>% select(x_um_dbscan, y_um_dbscan)
colnames(emb) = c("d_1","d_2") ; rownames(emb) = rownames(obj@meta.data)
obj[["dbscan"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "d_", assay = "RNA")
# Add KDE coords
emb = coords %>% mutate(across(everything(), ~ifelse(ratio > 1/3, NA, .))) %>% select(x_um_kde, y_um_kde)
colnames(emb) = c("k_1","k_2") ; rownames(emb) = rownames(obj@meta.data)
obj[["kde"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "k_", assay = "RNA")
# Add KDE-filtered DBSCAN coords
emb = obj@meta.data[,c("x_um","y_um")] ; colnames(emb) = c("s_1","s_2")
obj[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_", assay = "RNA")

plot <- plot_clusters(obj, reduction="dbscan")
make.pdf(plot, file.path(out_path,"DimPlotDBSCAN.pdf"), 7, 8)
plot <- plot_clusters(obj, reduction="kde")
make.pdf(plot, file.path(out_path,"DimPlotKDE.pdf"), 7, 8)
plot <- plot_clusters(obj, reduction="spatial")
make.pdf(plot, file.path(out_path,"DimPlot.pdf"), 7, 8)

plot <- plot_RNAvsSB(obj)
make.pdf(plot, file.path(out_path, "RNAvsSB.pdf"), 7, 8)

# Merge the PDF files
plotlist <- c(c("SB.pdf","beadplot.pdf","SBmetrics.pdf"),
              c("DBSCAN.pdf","KDE.pdf","DBSCANvsKDE.pdf","beadplots.pdf"),
              c("RNAmetrics.pdf","RNA.pdf","UMAP.pdf","DimPlot.pdf","DimPlotDBSCAN.pdf","DimPlotKDE.pdf","RNAvsSB.pdf"))
plotorder <- c(8, 9, 10, 1, 2, 4, 5, 6, 11, 12, 13, 14, 3, 7)
pdfs <- file.path(out_path, plotlist[plotorder])
pdfs %<>% keep(file.exists)
qpdf::pdf_combine(input=pdfs, output=file.path(out_path,"summary.pdf"))
file.remove(pdfs)

# Save the seurat object
qsave(obj, file.path(out_path,"seurat.qs"))
