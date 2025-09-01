### Wrapper script #############################################################
### Input: RNA directory
###        - GEX matrix (e.g. filtered_feature_bc_matrix.h5 or .h5ad)
###        - (optional) molecule_info.h5
###        - (optional) metrics_summary.csv or library_metrics.csv
### Input: path to SBcounts.h5 (output of spatial-count.jl)
### Output: obj.qs, summary.pdf (in addition to positioning.R outputs)
################################################################################

stopifnot(file.exists("positioning.R", "helpers.R", "plots.R"))
suppressMessages(source("helpers.R"))
suppressMessages(source("plots.R"))
suppressMessages(library(Seurat))
suppressMessages(library(qs))
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/usr/local/bin", sep = ":"))
if (Sys.which("Rscript") == "") {
  stop("Rscript not found")
}

# Load arguments
library(optparse)
arguments <- OptionParser(
  usage = "Usage: Rscript run-positioning.R rna_path sb_path out_path [options]",
  option_list = list(
    make_option("--cells", type="character", help = "Path to barcodes file"),
    make_option("--dropsift", action="store_true", help = "Add is_cell to obj"),
    make_option("--cores", type="integer", default=-1L, help = "The number of parallel processes to use [default: -1]"),
    make_option("--args",   type="character", default = "", help = "Passed to positioning.R")
  )
) %>% parse_args(positional_arguments=3)

# Load required arguments
rna_path <- arguments$args[[1]]
sb_path <- arguments$args[[2]]
out_path <- arguments$args[[3]]

print(g("RNA dir: {normalizePath(rna_path)}"))
stopifnot("RNA dir not found" = dir.exists(rna_path))
stopifnot("RNA dir empty" = len(list.files(rna_path)) > 0)

print(g("SB path: {normalizePath(sb_path)}"))
if (dir.exists(sb_path)) { sb_path %<>% file.path("SBcounts.h5") }
stopifnot("SB path is not an .h5" = str_sub(sb_path, -3, -1) == ".h5")
stopifnot("SBcounts.h5 not found" = file.exists(sb_path))

if (!dir.exists(out_path)) { dir.create(out_path, recursive = TRUE) }
print(g("Output dir: {normalizePath(out_path)}"))
stopifnot("Could not create output path" = dir.exists(out_path))

# Load optional arguments
cells <- arguments$options$cells
dropsift <- arguments$options$dropsift %>% {ifelse(is.null(.), FALSE, .)}
positioning_args <- arguments$options$args %>% trimws
cores <- arguments$options$cores %>% ifelse(.<1, parallelly::availableCores(), .) ; print(g("cores: {cores}"))
setDTthreads(cores)

print(g("Cell barcode whitelist: {cells}"))
print(g("DropSift: {dropsift}"))
if (nchar(positioning_args) > 0) {
  print(g("Positioning arguments override: {positioning_args}"))
}

cat("\n")

### Load the RNA ###############################################################

# Load the matrix
if (!is.null(cells) || dropsift) {
  cells <- file.path(rna_path, cells)
  stopifnot(file.exists(cells))
  
  # Load the raw matrix
  if (file.exists(file.path(rna_path, "raw_feature_bc_matrix.h5"))) {
    print("Loading raw_feature_bc_matrix.h5")
    mat <- Read10X_h5(file.path(rna_path, "raw_feature_bc_matrix.h5"))
  } else if (dir.exists(file.path(rna_path, "raw_feature_bc_matrix"))) {
    print("Loading raw_feature_bc_matrix")
    mat <- Read10X(file.path(rna_path, "raw_feature_bc_matrix"))
  } else if (any(str_ends(list.files(rna_path), ".h5ad"))) {
    print("Loading .h5ad")
    matrix_path <- list.files(rna_path, full.names=T) %>% keep(~str_sub(.,-5) == ".h5ad") %T>% {stopifnot(len(.)==1)}
    mat <- ReadAnnDataX(matrix_path, calledonly=FALSE)
  } else {
    rna_files <- list.files(rna_path) %>% paste0(collapse="\n")
    stop(g("Unable to find raw data:\n{rna_files}"), call. = FALSE)
  }
  
  # Filter out noise
  mat <- mat[, colSums(mat) >= 10]
  
} else {
  # Load the filtered matrix
  if (file.exists(file.path(rna_path, "filtered_feature_bc_matrix.h5"))) {
    print("Loading filtered_feature_bc_matrix.h5")
    mat <- Read10X_h5(file.path(rna_path, "filtered_feature_bc_matrix.h5"))
  } else if (dir.exists(file.path(rna_path, "filtered_feature_bc_matrix"))) {
    print("Loading filtered_feature_bc_matrix")
    mat <- Read10X(file.path(rna_path, "filtered_feature_bc_matrix"))
  } else if (file.exists(file.path(rna_path, "sample_filtered_feature_bc_matrix.h5"))) {
    print("Loading sample_filtered_feature_bc_matrix.h5")
    mat <- Read10X_h5(file.path(rna_path, "sample_filtered_feature_bc_matrix.h5"))
  } else if (dir.exists(file.path(rna_path, "sample_filtered_feature_bc_matrix"))) {
    print("Loading sample_filtered_feature_bc_matrix")
    mat <- Read10X(file.path(rna_path, "sample_filtered_feature_bc_matrix"))
  } else if (any(str_ends(list.files(rna_path), ".h5ad"))) {
    print("Loading .h5ad")
    matrix_path <-  list.files(rna_path, full.names=TRUE) %>% keep(~str_sub(.,-5) == ".h5ad") %T>% {stopifnot(len(.)==1)}
    mat <- ReadAnnDataX(matrix_path, calledonly=TRUE)
  } else {
    rna_files <- list.files(rna_path) %>% paste0(collapse="\n")
    stop(g("Unable to find filtered data:\n{rna_files}"), call. = FALSE)
  }
}

# Load the %intronic
if (file.exists(file.path(rna_path, "molecule_info.h5"))) {
  print("Loading molecule_info.h5")
  pct_intronic <- ReadIntronic(file.path(rna_path, "molecule_info.h5"), colnames(mat))
} else if (exists("matrix_path") && str_ends(matrix_path, ".h5ad")) {
  print("Loading .h5ad")
  pct_intronic <- ReadIntronic(matrix_path, colnames(mat))
} else {
  print("No Intronic file found")
  pct_intronic <- rep(NA, ncol(mat))
}
stopifnot(ncol(mat) == len(pct_intronic))

# Load the %mt
pct_mt <- colSums(mat[grepl("^MT-", rownames(mat), ignore.case=TRUE) |
                      grepl("_MT-", rownames(mat), ignore.case=TRUE),]) / colSums(mat)
pct_mt %<>% unname
stopifnot(ncol(mat) == len(pct_mt))

# DropSift
if (dropsift) {
  if (!file.exists(file.path(rna_path, "dropsift.csv"))) {
    stopifnot(!is.na(pct_intronic), !is.na(pct_mt))
    stopifnot(any(pct_intronic > 0), any(pct_mt > 0))
    
    print("Running DropSift")
    library(DropSift)
    cellFeatures <- data.frame(cell_barcode = colnames(mat),
                               pct_intronic = pct_intronic,
                               pct_mt = pct_mt)
    svmNucleusCaller <- SvmNucleusCaller(cellFeatures=cellFeatures,
                                         dgeMatrix=mat,
                                         useCBRBFeatures = FALSE)
    is_cell <- svmNucleusCaller$cell_features$is_cell %>% as.logical
    is_cell_prob <- svmNucleusCaller$cell_features$is_cell_prob
    dropsift <- data.frame(cell_barcode = colnames(mat)[is_cell==TRUE],
                           is_cell_prob = is_cell_prob[is_cell==TRUE])
    write.table(dropsift, file.path(rna_path, "dropsift.csv"),
                sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
  
  print("Loading dropsift.csv")
  dropsift <- read.csv(file.path(rna_path, "dropsift.csv"), header=FALSE)
  dropsift %<>% setNames(c("cell_barcode", "is_cell_prob"))
  stopifnot(dropsift$cell_barcode %in% colnames(mat))
  
  # Plot dropsift selection
  data.table(umi=colSums(mat),
             pct_intronic=pct_intronic,
             called=colnames(mat) %in% dropsift$cell_barcode) %>% 
    plot_cellcalling("DropSift") %>% make.pdf(file.path(out_path, "RNAlibrary2.pdf"), 7, 8)
}

# Load cell barcode whitelist
if (!is.null(cells)) {
  print(g("Loading {basename(cells)}"))
  cb_whitelist <- readLines(cells)
  stopifnot(cb_whitelist %in% colnames(mat))
  
  # Plot cell selection
  data.table(umi=colSums(mat),
             pct_intronic=pct_intronic,
             called=colnames(mat) %in% cb_whitelist) %>% 
    plot_cellcalling(basename(cells)) %>% make.pdf(file.path(out_path, "RNAlibrary.pdf"), 7, 8)
  
  if (class(dropsift) == "data.frame") {cb_whitelist %<>% union(dropsift$cell_barcode) %>% sort}
  
} else if (class(dropsift) == "data.frame") {
  cb_whitelist <- dropsift$cell_barcode
} else {
  cb_whitelist <- colnames(mat)
}

# Subset to called cells
stopifnot(cb_whitelist %in% colnames(mat))
mask <- match(cb_whitelist, colnames(mat))
mat <- mat[,mask]
pct_intronic <- pct_intronic[mask]
pct_mt <- pct_mt[mask]
rm(mask)

# Create Seurat object
obj <- CreateSeuratObject(mat, project="Slide-tags")
rm(mat) ; invisible(gc())

# Add metadata
if (class(dropsift) == "data.frame") {
  obj$is_cell <- colnames(obj) %in% dropsift$cell_barcode
  #obj$is_cell_prob <- dropsift$is_cell_prob[match(colnames(obj), dropsift$cell_barcode)]
}
obj$pct_intronic <- pct_intronic
obj$pct_mt <- pct_mt

# Add .h5ad annotations
if (exists("matrix_path") && str_ends(matrix_path, ".h5ad")) {
  Misc(obj, "obs") <- rhdf5::h5read(matrix_path, "/obs")
  Misc(obj, "var") <- rhdf5::h5read(matrix_path, "/var")
}

# Plot library metrics + add to object
if (file.exists(file.path(rna_path, "metrics_summary.csv"))) { # 10X
  print("Loading metrics_summary.csv")
  df <- file.path(rna_path, "metrics_summary.csv") %>% 
    read.table(header=FALSE, sep=",", comment.char="") %>% t %>% as.data.frame
} else if (any(str_ends(list.files(rna_path), "_library_metrics.csv"))) { # Optimus
  print("Loading library_metrics.csv")
  df <- list.files(rna_path, full.names=TRUE) %>% 
    keep(~str_sub(., -20) == "_library_metrics.csv") %T>% 
    {stopifnot(len(.) == 1)} %>% 
    read.table(header=FALSE, sep=",", comment.char="")
} else {
  df <- data.frame()
}
if (ncol(df) > 0) {
  Misc(obj, "gex_metrics") <- setNames(df[,2], df[,1])
  plot_metrics_csv(df) %>% make.pdf(file.path(out_path, "RNAmetrics.pdf"), 7, 8)
}

# Normalize, HVG, Scale, PCA, Neighbors, Clusters, UMAP
obj %<>%
  Seurat::NormalizeData(verbose=FALSE) %>%
  Seurat::FindVariableFeatures(verbose=FALSE) %>%
  Seurat::ScaleData(verbose=FALSE) %>%
  Seurat::RunPCA(npcs=50, verbose=FALSE) %>%
  Seurat::FindNeighbors(dims=1:30, verbose=FALSE) %>%
  Seurat::FindClusters(resolution=0.8, verbose=FALSE) %>%
  Seurat::RunUMAP(dims=1:30, umap.method="uwot", metric="cosine", verbose=FALSE)

# Plot UMAPs
plot_umaps(obj) %>% make.pdf(file.path(out_path, "RNAplot.pdf"), 7, 8)

print(g("Loaded {ncol(obj)} cells"))

### DBSCAN #####################################################################

# Write cell barcode sequences
colnames(obj) %>% trim_10X_CB %>% writeLines(file.path(out_path, "cb_whitelist.txt"))

# Assign a position to each whitelist cell
print(g("\nRunning positioning.R"))
system(g("Rscript --vanilla positioning.R {sb_path} {file.path(out_path, 'cb_whitelist.txt')} {out_path} --cores={cores} {positioning_args}"))

stopifnot(file.exists(file.path(out_path, "matrix.csv.gz"),
                      file.path(out_path, "spatial_metadata.json"),
                      file.path(out_path, "SBsummary.pdf"),
                      file.path(out_path, "DBSCANs.pdf"),
                      file.path(out_path, "coords.csv"),
                      file.path(out_path, "coords2.csv")))

### Add to Seurat object #######################################################

Misc(obj, "sb_metrics") <- jsonlite::fromJSON(file.path(out_path, "spatial_metadata.json"))

coords <- fread(file.path(out_path, "coords.csv"), header=TRUE, sep=",")
coords2 <- fread(file.path(out_path, "coords2.csv"), header=TRUE, sep=",")
stopifnot(nrow(coords) == ncol(obj), coords$cb == trim_10X_CB(colnames(obj)))
stopifnot(nrow(coords2) == ncol(obj), coords2$cb == trim_10X_CB(colnames(obj)))
Misc(obj, "coords") <- coords
Misc(obj, "coords2") <- coords2

# Add metadata
obj %<>% AddMetaData(cbind(coords2[,.(x=x,y=y)],
                           coords[,.(dbscan_clusters=clusters)],
                           coords2[,.(dbscan_score=score)]))

# Create spatial reduction
emb <- obj@meta.data %>% transmute(x = if_else(dbscan_clusters == 1, x, NA_real_),
                                   y = if_else(dbscan_clusters == 1, y, NA_real_))
colnames(emb) <- c("spatial_1", "spatial_2") ; rownames(emb) = rownames(obj@meta.data)
obj[["spatial"]] <- CreateDimReducObject(embeddings=as.matrix(emb), key="spatial_", assay="RNA")

plot_clusters(obj, reduction="spatial") %>% make.pdf(file.path(out_path, "DimPlot.pdf"), 7, 8)
plot_RNAvsSB(obj) %>% make.pdf(file.path(out_path, "RNAvsSB.pdf"), 7, 8)

# Merge the PDF files
plotlist <- c("RNAmetrics.pdf","RNAlibrary.pdf","RNAlibrary2.pdf","RNAplot.pdf",
              "SBsummary.pdf",
              "DimPlot.pdf", "RNAvsSB.pdf", "DBSCANs.pdf")
pdfs <- file.path(out_path, plotlist)
pdfs %<>% keep(file.exists)
qpdf::pdf_combine(input=pdfs, output=file.path(out_path, "summary.pdf"))
file.remove(pdfs)

# Save the seurat object
qsave(obj, file.path(out_path, "seurat.qs"))
print("Done!")
