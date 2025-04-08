library(Seurat)
library(qs)

# Helper files required for loading spatial data and positioning cells
stopifnot(file.exists("positioning.R", "helpers.R"))
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/usr/local/bin", sep = ":"))
stopifnot(Sys.which("Rscript") != "")
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
} else if (any(str_ends(list.files(RNApath), ".h5ad"))) {
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

### DBSCAN #####################################################################

cb_whitelist = unname(obj$cb)
writeLines(cb_whitelist, file.path(out_path, "cb_whitelist.txt"))

# Assign a position to each whitelist cell
print("Positioning cells...")
system(g("Rscript --vanilla positioning.R {SBpath} {file.path(out_path, 'cb_whitelist.txt')} {out_path}"))

stopifnot(file.exists(file.path(out_path, "matrix.csv.gz"),
                      file.path(out_path, "spatial_metadata.json"),
                      file.path(out_path, "coords_global.csv"),
                      file.path(out_path, "coords_dynamic.csv")))
Misc(obj, "SB_metadata") <- jsonlite::fromJSON(file.path(out_path, "spatial_metadata.json"))

coords_global <- fread(file.path(out_path, "coords_global.csv"), header=TRUE, sep=",")
coords_dynamic <- fread(file.path(out_path, "coords_dynamic.csv"), header=TRUE, sep=",")
Misc(obj, "coords_global") <- coords_global
Misc(obj, "coords_dynamic") <- coords_dynamic

# Create spatial reduction
coords <- obj@meta.data %>% select(cb) %>% left_join(coords_global, by="cb") %>% select(x,y,clusters)
obj %<>% AddMetaData(coords)

emb <- coords %>% select(x,y)
colnames(emb) <- c("d_1","d_2") ; rownames(emb) = rownames(obj@meta.data)
obj[["dbscan"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "d_", assay = "RNA")
colnames(emb) <- c("s_1","s_2") ; rownames(emb) = rownames(obj@meta.data)
obj[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_", assay = "RNA")

plot_clusters(obj, reduction="dbscan") %>% make.pdf(file.path(out_path,"DimPlot.pdf"), 7, 8)
plot_RNAvsSB(obj) %>% make.pdf(file.path(out_path, "RNAvsSB.pdf"), 7, 8)

# Merge the PDF files
plotlist <- c(c("RNAmetrics.pdf","RNAcalls.pdf","RNAumap.pdf"),
              c("SBlibrary.pdf","SBmetrics.pdf","SBplot.pdf"),
              c("GDBSCANopt.pdf", "GDBSCAN1.pdf", "DDBSCANxy.pdf"),
              c("DimPlot.pdf", "RNAvsSB.pdf", "DBSCAN.pdf"))
pdfs <- file.path(out_path, plotlist)
pdfs %<>% keep(file.exists)
qpdf::pdf_combine(input=pdfs, output=file.path(out_path,"summary.pdf"))
file.remove(pdfs)

# Save the seurat object
qsave(obj, file.path(out_path, "seurat.qs"))
print("Done!")
