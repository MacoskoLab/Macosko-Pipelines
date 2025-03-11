################################################################################
### Input: SBcounts.h5 (output of spatial-count.jl)
### Input: cb_whitelist.txt (newline-delimited list of cell barcodes to position)
### Output: coords.csv (+metadata)
################################################################################
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

source("helpers.R")

# Load arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  sb_path <- args[[1]]
  cb_path <- args[[2]]
  out_path <- "."
} else if (length(args) == 3) {
  sb_path <- args[[1]]
  cb_path <- args[[2]]
  out_path <- args[[3]]
} else {
  stop("Usage: Rscript load_matrix.R SBcounts_path cb_whitelist_path [output_path]", call. = FALSE)
}

# Check input arguments
if (!dir.exists(out_path)) {dir.create(out_path, recursive=T)}
stopifnot(filter(rhdf5::h5ls(sb_path), group=="/matrix")$name == c("cb_index", "reads", "sb_index", "umi"))

# Load the spatial barcode count matrix
f <- function(p){return(rhdf5::h5read(sb_path, p))}
dt <- ReadSpatialMatrix(f) # stopifnot(nrow(dt) == nrow(unique(dt[,.(cb, umi, sb)])))
metadata <- ReadSpatialMetadata(f)
print(g("{add.commas(sum(dt$reads))} spatial barcode reads loaded"))

# load the CB whitelist
cb_whitelist = readLines(cb_path)

remap <- determine_remap_10X_CB(cb_whitelist, dt)
if (remap) { cb_whitelist %<>% remap_10X_CB }
metadata$SB_info$remap_10X_CB = remap

stopifnot(class(cb_whitelist) == "character")
stopifnot(!duplicated(cb_whitelist))
stopifnot(uniqueN(nchar(cb_whitelist)) == 1)
stopifnot(map_lgl(strsplit(cb_whitelist, ""), ~all(. %in% c("A","C","G","T"))))
print(g("{len(cb_whitelist)} cell barcodes loaded"))

### Fuzzy matching #############################################################

print("Performing fuzzy cell-barcode matching")
setnames(dt, "cb", "cb_fuzzy")
df <- data.table(cb=cb_whitelist)[, .(cb_fuzzy=listHD1neighbors(cb), match="fuzzy"), cb
                                  ][!cb_fuzzy %in% cb_whitelist
                                    ][cb_fuzzy %in% dt$cb_fuzzy]
df[, `:=`(cb = ifelse(.N > 1, NA, cb),
       match = ifelse(.N > 1, "ambig", match)), cb_fuzzy]
df %<>% unique
df <- rbindlist(list(df, data.table(cb=cb_whitelist,
                              cb_fuzzy=cb_whitelist,
                                 match="exact")[cb_fuzzy %in% dt$cb_fuzzy]))
stopifnot(df$cb_fuzzy %in% dt$cb_fuzzy)
df[, cb_fuzzy := factor(cb_fuzzy, levels = levels(dt$cb_fuzzy))]
stopifnot(sort(unique(df$cb), na.last=NA) == sort(cb_whitelist))
df[, cb := factor(cb, levels = cb_whitelist)]
stopifnot(df$match %in% c("exact", "fuzzy", "ambig"), !is.na(df$match))
df[, match := factor(match, levels = c("exact", "fuzzy", "ambig"))]
stopifnot(!any(duplicated(df$cb_fuzzy)), !is.na(df$cb_fuzzy))
stopifnot(!xor(df$match=="ambig", is.na(df$cb)))
stopifnot(!xor(df$match=="exact", df$cb_fuzzy %in% cb_whitelist))
stopifnot(levels(df$cb_fuzzy) == levels(dt$cb_fuzzy))
stopifnot(is.null(key(dt)), is.null(key(df)), is.null(indices(dt)), is.null(indices(df)))
# Perform the match
dt %<>% merge(df, by = "cb_fuzzy", all.x = TRUE, sort = FALSE)
stopifnot(is.null(key(dt)), is.null(indices(dt)))

# Record metadata
metadata$CB_matching <- dt[, .(reads=sum(reads)), match][order(-reads)] %>% 
            {setNames(.$reads, .$match %>% as.character %>% replace_na("none"))}
stopifnot(sum(dt$reads) == sum(metadata$CB_matching))

# Clean up
dt[, match := NULL]
setnames(dt, "cb_fuzzy", "cb_raw")
setcolorder(dt, c("cb_raw", "cb", "umi", "sb", "reads"))
rm(df) ; invisible(gc())

# Collapse post-matched barcodes
dt <- rbindlist(list(dt[is.na(cb)],
                     dt[!is.na(cb), .(reads=sum(reads)), .(cb,umi,sb)]),
                fill = TRUE)

### Load puck ##################################################################

# Load the puck
print("Loading the puck")
res <- ReadPuck(f)
puckdf <- res[[1]]
metadata$puck_info %<>% c(res[[2]])
epsilon <- res[[3]]
rm(res) ; invisible(gc())

# Filter reads with a low-quality spatial barcode
print("Removing low-quality spatial barcodes")
dt[, m := sb %in% puckdf$sb]
metadata$SB_filtering %<>% c(reads_lqsb=dt[m == FALSE, sum(reads)])
dt <- dt[m == TRUE]
dt[, m := NULL]

# Page 4
plot_SBlibrary(dt, f) %>% make.pdf(file.path(out_path, "SBlibrary.pdf"), 7, 8)

# Page 5
plot_SBplot(dt, puckdf) %>% make.pdf(file.path(out_path, "SBplot.pdf"), 7, 8)

# Delete cell barcodes for cells that were not called
print("Removing non-whitelist cells")
metadata$SB_filtering %<>% c(reads_nocb = dt[is.na(cb), sum(reads)])
metadata$SB_info$UMI_pct_in_called_cells <- round(nrow(dt[!is.na(cb)])/nrow(dt)*100, digits=2) %>% paste0("%")
metadata$SB_info$sequencing_saturation <- round((1 - nrow(dt) / sum(dt$reads)) * 100, 2) %>% paste0("%")
dt <- dt[!is.na(cb)]
dt[, cb_raw := NULL]
invisible(gc())
stopifnot(nrow(dt) == nrow(unique(dt[,.(cb,umi,sb)])))

# Remove chimeric UMIs
print("Removing chimeras")
dt[, m := reads == max(reads) & sum(reads == max(reads)) == 1, .(cb, umi)]
metadata$SB_filtering %<>% c(reads_chimeric = dt[m==FALSE, sum(reads)])
dt <- dt[m == TRUE]
dt[, m := NULL]

# Aggregate reads into UMIs
print("Counting UMIs")
metadata$SB_filtering %<>% c(reads_final = dt[,sum(reads)])
dt <- dt[, .(umi=.N), .(cb, sb)]
metadata$SB_filtering %<>% c(UMIs_final = dt[, sum(umi)])
invisible(gc())

# Merge spatial barcode count matrix with puck coordinates
print("Joining puck coordinates")
stopifnot(dt$sb %in% puckdf$sb)
stopifnot(is.factor(dt$sb), is.factor(puckdf$sb))
dt <- merge(dt, puckdf, by = "sb", all = FALSE, sort = FALSE)[order(cb, -umi)]
setcolorder(dt, c("cb","umi","sb","x","y"))
metadata$puck_info$umi_final <- map2_int(metadata$puck_info$puck_boundaries %>% head(-1),
                                         metadata$puck_info$puck_boundaries %>% tail(-1),
                                         ~dt[x >= .x & x <= .y, sum(umi)])
rm(puckdf) ; invisible(gc())

# Page 6
plot_SBmetrics(metadata) %>% make.pdf(file.path(out_path, "SBmetrics.pdf"), 7, 8)

# Write intermediate output
print("Writing intermediate matrix")
fwrite(dt, file.path(out_path, "matrix.csv.gz"), quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE, compress="gzip")
metadata %>% map(as.list) %>% jsonlite::toJSON(pretty=T) %>% writeLines(file.path(out_path, "spatial_metadata.json"))
dt[, sb := NULL]

stopifnot(dt$cb %in% cb_whitelist)
stopifnot(dt$umi > 0)
stopifnot(!any(is.na(dt)))
stopifnot(sort(names(dt)) == c("cb", "umi", "x", "y"))

### DBSCAN #####################################################################

xlims = range(dt$x) ; xrange = max(dt$x) - min(dt$x)
ylims = range(dt$y) ; yrange = max(dt$y) - min(dt$y)
data.list <- split(dt, by="cb", drop=TRUE, keep.by=FALSE)
for (dl in data.list) {stopifnot(nrow(dl) > 0)}
rm(dt) ; invisible(gc())
print(g("Running positioning on {len(data.list)} cells"))

# Prepare for positioning
library(future)
library(furrr)
myoptions <- furrr::furrr_options(packages=c("data.table"), seed=TRUE, scheduling=1)
future::plan(future::multisession, workers=parallelly::availableCores())
mydbscan <- function(dl, eps, minPts) {dbscan::dbscan(x = dl[,.(x,y)],
                                                      eps = eps,
                                                      minPts = minPts,
                                                      weights = dl$umi,
                                                      borderPoints = FALSE)$cluster}

# Dynamic DBSCAN
ms <- 1 # minPts step size (must be integer for IRanges)
eps = epsilon * 15
res <- furrr::future_map(data.list, function(dl){
  # base case
  minpts <- ms
  a1 <- mydbscan(dl, eps, minpts)
  a2 <- a1
  if (max(a2) <= 1) {vec <- a1 ; minpts2 <- minpts}
  
  # loop
  while (max(a2) > 0) {
    minpts <- minpts + ms
    a1 <- a2
    a2 <- mydbscan(dl, eps, minpts)
    if (max(a1) > 1 & max(a2) <= 1) {vec <- a1 ; minpts2 <- minpts}
  }
  
  # [[1]] are the DBSCAN clusters at highest minPts that produces 2+ clusters
  # [[2]] are the DBSCAN clusters at highest minPts that produces 1 cluster
  # [[3]][[1]] is the highest minPts that produces DBSCAN=2, plus 1 step
  # [[3]][[2]] is the lowest minPts that produces DBSCAN=0, minus 1 step
  return(list(vec, a1, c(minpts2, minpts-ms)))
}, .options = myoptions)
invisible(map2(data.list, res, ~.x[, c2 := .y[[1]]]))
invisible(map2(data.list, res, ~.x[, c1 := .y[[2]]]))
mranges <- IRanges::IRanges(start=map_int(res, ~.[[3]][[1]]),
                            end=map_int(res, ~.[[3]][[2]])) %>% 
           IRanges::coverage() %>% as.integer
rm(res)

# DBSCAN
minPts <- which.max(mranges) * ms
coords_dbscan <- map(data.list, function(dl){
  d <- mydbscan(dl, eps, minPts)
  if (max(d) == 1) {
    return(dl[d==1, .(x=weighted.mean(x, umi),
                      y=weighted.mean(y, umi),
                      sumi=sum(umi),
                      umi=sum(dl$umi),
                      SNR=sum(umi)/sum(dl$umi))])
  } else {
    return(dl[,.(x=NA,y=NA,sumi=NA,umi=sum(umi),SNR=NA)])
  }
}) %>% rbindlist

fwrite(coords_dbscan, file.path(out_path, "coords.csv"))
