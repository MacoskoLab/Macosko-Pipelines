# Helper files required for loading spatial data and positioning cells
stopifnot(file.exists("positioning.R", "helpers.R"))
suppressMessages(source("helpers.R"))
suppressMessages(library(Seurat))
suppressMessages(library(qs))
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/usr/local/bin", sep = ":"))
if (Sys.which("Rscript") == "") {
  stop("Rscript not found")
}

# setwd("/broad/macosko/mshabet/testing/positioning")
# rna_path = "RNA"
# sb_path = "RNA"
# out_path = "output"
# positioning_args = ""

# rna_path="GG/slidetag_gg/outs"
# sb_path="GG" 
# out_path="GGo"

# Load arguments
library(optparse)
arguments <- OptionParser(
  usage = "Usage: Rscript run-positioning.R rna_path sb_path output_path [options]",
  option_list = list(
    make_option("--knn", type="integer", help = "Passed to positioning.R"),
    make_option("--cmes", type="double", help = "Passed to positioning.R"),
    make_option("--prob", type="double", help = "Passed to positioning.R")
  )
) %>% parse_args(positional_arguments=3)

rna_path <- arguments$args[[1]]
sb_path <- arguments$args[[2]]
out_path <- arguments$args[[3]]
positioning_args <- paste0(arguments$options$knn %>% {ifelse(is.null(.), "", g(" --knn={.}"))},
                           arguments$options$mask %>% {ifelse(is.null(.), "", g(" --cmes={.}"))},
                           arguments$options$prob %>% {ifelse(is.null(.), "", g(" --prob={.}"))}) %>% trimws
rm(arguments)

# Check arguments
print(g("RNA dir: {normalizePath(rna_path)}"))
stopifnot("RNA dir not found" = dir.exists(rna_path))
stopifnot("RNA dir empty" = len(list.files(rna_path)) > 0)

print(g("SB file: {normalizePath(sb_path)}"))
if (dir.exists(sb_path)) { sb_path = file.path(sb_path, "SBcounts.h5") }
stopifnot("SB path is not SBcounts.h5" = str_sub(sb_path, -11, -1) == "SBcounts.h5")
stopifnot("SBcounts.h5 not found" = file.exists(sb_path))

if (!dir.exists(out_path)) { dir.create(out_path, recursive = TRUE) }
print(g("Output dir: {normalizePath(out_path)}"))
stopifnot("Could not create output path" = dir.exists(out_path))

if (nchar(positioning_args) > 0) {
  print(g("Positioning arguments override: {positioning_args}"))
}

cat("\n")

# Determine the RNA method
if (dir.exists(file.path(rna_path, "filtered_feature_bc_matrix"))) {
  rna_source = "10X"
  matrix_path = file.path(rna_path, "filtered_feature_bc_matrix")
  intronic_path = file.path(rna_path, "molecule_info.h5")
  metrics_path = file.path(rna_path, "metrics_summary.csv")
} else if (file.exists(file.path(rna_path, "filtered_feature_bc_matrix.h5"))) {
  rna_source = "10X_h5"
  matrix_path = file.path(rna_path, "filtered_feature_bc_matrix.h5")
  intronic_path = file.path(rna_path, "molecule_info.h5")
  metrics_path = file.path(rna_path, "metrics_summary.csv")
} else if (any(str_ends(list.files(rna_path), ".h5ad"))) {
  rna_source = "Optimus"
  matrix_path <-  list.files(rna_path, full.names=T) %>% keep(~str_sub(.,-5) == ".h5ad") %T>% {stopifnot(len(.)==1)}
  intronic_path <- list.files(rna_path, full.names=T) %>% keep(~str_sub(.,-5) == ".h5ad") %T>% {stopifnot(len(.)<=1)}
  metrics_path <- list.files(rna_path, full.names=T) %>% keep(~str_sub(.,-20) == "_library_metrics.csv") %T>% {stopifnot(len(.)<=1)}
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
rm(mat) ; invisible(gc())

# Add metadata
obj[["logumi"]] <- log10(obj$nCount_RNA+1)
obj[["cb"]] <- map_chr(colnames(obj), ~sub("-[0-9]*$", "", .)) %T>% {stopifnot(!duplicated(.))}
obj[["cb_index"]] <- 1:ncol(obj)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(MT-|mt-)")

# Compute %intronic
obj$pct.intronic <- ReadIntronic(intronic_path, unname(obj$cb))

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
plot_cellcalling(intronic_path, obj$cb) %>% make.pdf(file.path(out_path,"RNAlibrary.pdf"), 7, 8)

# Page 3
plot_umaps(obj) %>% make.pdf(file.path(out_path,"RNAplot.pdf"), 7, 8)

### DBSCAN #####################################################################

cb_whitelist = unname(obj$cb)
writeLines(cb_whitelist, file.path(out_path, "cb_whitelist.txt"))

# Assign a position to each whitelist cell
print(g("\nRunning positioning.R"))
system(g("Rscript --vanilla positioning.R {sb_path} {file.path(out_path, 'cb_whitelist.txt')} {out_path}"))

stopifnot(file.exists(file.path(out_path, "matrix.csv.gz"),
                      file.path(out_path, "spatial_metadata.json"),
                      file.path(out_path, "coords.csv"),
                      file.path(out_path, "coords2.csv")))

### Add to Seurat object #######################################################

Misc(obj, "SB_metadata") <- jsonlite::fromJSON(file.path(out_path, "spatial_metadata.json"))

coords_global <- fread(file.path(out_path, "coords.csv"), header=TRUE, sep=",")
#coords_dynamic <- fread(file.path(out_path, "coords2.csv"), header=TRUE, sep=",")
Misc(obj, "coords_global") <- coords_global
#Misc(obj, "coords_dynamic") <- coords_dynamic

# Create spatial reduction
coords <- obj@meta.data %>% select(cb) %>% left_join(coords_global, by="cb") %>% select(x,y,clusters)
obj %<>% AddMetaData(coords)

emb <- coords %>% select(x,y)
colnames(emb) <- c("d_1","d_2") ; rownames(emb) = rownames(obj@meta.data)
obj[["dbscan"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "d_", assay = "RNA")
#colnames(emb) <- c("s_1","s_2") ; rownames(emb) = rownames(obj@meta.data)
#obj[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_", assay = "RNA")

plot_clusters(obj, reduction="dbscan") %>% make.pdf(file.path(out_path,"DimPlot.pdf"), 7, 8)
plot_RNAvsSB(obj) %>% make.pdf(file.path(out_path, "RNAvsSB.pdf"), 7, 8)

# Merge the PDF files
plotlist <- c("RNAmetrics.pdf","RNAlibrary.pdf","RNAplot.pdf",
              "SBlibrary.pdf","SBplot.pdf","SBmetrics.pdf",
              "GDBSCANopt.pdf", "GDBSCAN1.pdf", "GDBSCAN2.pdf",
              "DimPlot.pdf", "RNAvsSB.pdf", "GDBSCANs.pdf")
pdfs <- file.path(out_path, plotlist)
pdfs %<>% keep(file.exists)
qpdf::pdf_combine(input=pdfs, output=file.path(out_path, "summary.pdf"))
#file.remove(pdfs)

# Save the seurat object
qsave(obj, file.path(out_path, "seurat.qs"))
print("Done!")
