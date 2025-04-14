################################################################################
### Input: SBcounts.h5 (output of spatial-count.jl)
### Input: cb_whitelist.txt (newline-delimited list of cell barcodes to position)
### Output: coords.csv (placements+metadata)
################################################################################

suppressMessages(source("helpers.R"))

# Load arguments
library(optparse)
arguments <- OptionParser(
  usage = "Usage: Rscript positioning.R SBcounts_path CBwhitelist_path output_path [options]",
  option_list = list(
    make_option(c("-k", "--knn"), type="integer", default=51, help = "Number of bead neighbors used to compute eps [default: %default]", metavar="K")
  )
) %>% parse_args(positional_arguments=3)

sb_path <- arguments$args[[1]]
cb_path <- arguments$args[[2]]
out_path <- arguments$args[[3]]
k <- arguments$options$knn
rm(arguments)

# Print arguments
print(g("sb_path: {sb_path}"))
print(g("cb_path: {cb_path}"))
print(g("out_path: {out_path}"))
print(g("k: {k}"))

# Check input arguments
if (!dir.exists(out_path)) {dir.create(out_path, recursive=TRUE)}
stopifnot(filter(rhdf5::h5ls(sb_path), group=="/matrix")$name == c("cb_index", "reads", "sb_index", "umi"))

# Load the spatial barcode count matrix
f <- function(p){return(rhdf5::h5read(sb_path, p))}
dt <- ReadSpatialMatrix(f)
metadata <- ReadSpatialMetadata(f)
print(g("{add.commas(sum(dt$reads))} spatial barcode reads loaded"))

# load the CB whitelist
cb_whitelist <- readLines(cb_path)

remap <- determine_remap_10X_CB(cb_whitelist, dt)
if (remap) { cb_whitelist %<>% remap_10X_CB }
metadata$SB_info$remap_10X_CB <- remap

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
stopifnot(na.omit(unique(df$cb)) %in% cb_whitelist)
df[, cb := factor(cb, levels = cb_whitelist)]
stopifnot(df$match %in% c("exact", "fuzzy", "ambig"), !is.na(df$match))
df[, match := factor(match, levels = c("exact", "fuzzy", "ambig"))]
stopifnot(!any(duplicated(df$cb_fuzzy)), !is.na(df$cb_fuzzy))
stopifnot(!xor(df$match=="ambig", is.na(df$cb)))
stopifnot(!xor(df$match=="exact", df$cb_fuzzy %in% cb_whitelist))
stopifnot(levels(df$cb_fuzzy) == levels(dt$cb_fuzzy))
stopifnot(is.null(key(dt)), is.null(key(df)), is.null(indices(dt)), is.null(indices(df)))
# Perform the match
dt <- merge(dt, df, by = "cb_fuzzy", all.x = TRUE, sort = FALSE)
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
rm(res) ; invisible(gc())

# Compute eps scale using the kth neighbor
eps <- RANN::nn2(data = puckdf[, .(x,y)],
                 query = puckdf[sample(.N, 10000), .(x,y)],
                 k = k)$nn.dists[,k] %>% median
print(g("eps: {eps}"))

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
invisible(gc())

# Page 6
plot_SBmetrics(metadata) %>% make.pdf(file.path(out_path, "SBmetrics.pdf"), 7, 8)

# Write intermediate output
print("Writing intermediate matrix")
fwrite(dt, file.path(out_path, "matrix.csv.gz"), quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE, compress="gzip")
metadata %>% map(as.list) %>% jsonlite::toJSON(pretty=TRUE) %>% writeLines(file.path(out_path, "spatial_metadata.json"))
dt[, sb := NULL]

stopifnot(dt$cb %in% cb_whitelist)
stopifnot(dt$umi > 0)
stopifnot(!any(is.na(dt)))
stopifnot(sort(names(dt)) == c("cb", "umi", "x", "y"))

### DBSCAN #####################################################################

data.list <- split(dt, by="cb", drop=TRUE, keep.by=FALSE)
for (dl in data.list) {stopifnot(nrow(dl) > 0, dl$umi > 0)}
rm(dt) ; invisible(gc())
print(g("Running positioning on {len(data.list)} cells"))

# prepare for positioning
library(future)
library(furrr)
options(future.globals.maxSize = 1024 * 1024 * 1024)
myoptions <- furrr::furrr_options(packages=c("data.table"), seed=TRUE, scheduling=1)
future::plan(future::multisession, workers=parallelly::availableCores())
mydbscan <- function(dl, eps, minPts) {dbscan::dbscan(x = dl[,.(x,y)],
                                                      eps = eps,
                                                      minPts = minPts,
                                                      weights = dl$umi,
                                                      borderPoints = TRUE)$cluster}

### Run DBSCAN optimization ###
ms <- 1L # minPts step size
mranges <- furrr::future_map(data.list, function(dl) {
  # base case
  v2 <- mydbscan(dl, eps, 1L * ms)
  if (max(v2) == 0) {return(data.table(i2=0L, i1=0L))}
  if (max(v2) == 1) {i2=0L}
  
  # loop
  i <- 1L
  while (max(v2) > 0) {
    i <- i + 1L
    v1 <- v2
    v2 <- mydbscan(dl, eps, i * ms)
    if (max(v1) > 1 & max(v2) <= 1) {i2=i-1L}
  }
  
  # i2 is the index at the highest minPts that produces DBSCAN=2+
  # i1 is the index at the highest minPts that produces DBSCAN=1
  return(data.table(i2=i2, i1=i-1L))
}, .options = myoptions) %>% rbindlist
stopifnot(nrow(mranges) == len(data.list))

# TODO: Walk it back by filtering nearby DBSCAN=2

### Global DBSCAN ###

# Compute minPts that places the greatest number of cells
i_opt <- IRanges::IRanges(start=mranges$i2+1L, end=mranges$i1) %>% 
  IRanges::coverage() %>% as.integer %>% 
  {which(.==max(.))} %>% tail(1)
minPts <- i_opt * ms
print(g("Global optimum minPts: {minPts}"))

# Rerun DBSCAN with the optimal minPts
for (i in seq_along(data.list)) {
  data.list[[i]][, cluster := mydbscan(data.list[[i]], eps, minPts)]
}

# Assign the centroid via a weighted mean
coords_global <- imap(data.list, function(dl, cb){
  ret <- dl[,.(cb=cb,
               umi=sum(umi),
               beads=.N,
               max=max(umi),
               clusters=max(cluster))]
  
  # Add DBSCAN=0 data
  ret0 <- dl[cluster==0, .(umi0=sum(umi),
                           beads0=.N,
                           max0=max(umi, default=0))]
  ret[, names(ret0) := ret0]
  
  # Add DBSCAN=1 data
  ret1 <- dl[cluster==1, .(x1=weighted.mean(x, umi),
                           y1=weighted.mean(y, umi),
                           rmsd1=sqrt(weighted.mean((x-weighted.mean(x,umi))^2+(y-weighted.mean(y,umi))^2, umi)),
                           umi1=sum(umi),
                           beads1=.N,
                           max1=max(umi, default=0),
                           h1=h_index(umi))]
  ret[, names(ret1) := ret1]
  
  # Add DBSCAN=2 data
  ret2 <- dl[cluster==2, .(x2=weighted.mean(x, umi),
                           y2=weighted.mean(y, umi),
                           rmsd2=sqrt(weighted.mean((x-weighted.mean(x,umi))^2+(y-weighted.mean(y,umi))^2, umi)),
                           umi2=sum(umi),
                           beads2=.N,
                           max2=max(umi, default=0),
                           h2=h_index(umi))]
  ret[, names(ret2) := ret2]
  
  return(ret)
  
}) %>% rbindlist(use.names=TRUE, fill=TRUE) # ignore.attr=TRUE
coords_global[, `:=`(eps=eps, minPts=minPts)]
coords_global[clusters==1, `:=`(x=x1, y=y1)]
print(g("Placed: {round(coords_global[,sum(clusters==1)/.N]*100, 2)}%"))

# Plots
plot_gdbscan_opt(coords_global, mranges) %>% make.pdf(file.path(out_path, "GDBSCANopt.pdf"), 7, 8)
plot_gdbscan_1(coords_global) %>% make.pdf(file.path(out_path, "GDBSCAN1.pdf"), 7, 8)
plot_gdbscan_2(coords_global) %>% make.pdf(file.path(out_path, "GDBSCAN2.pdf"), 7, 8)

# Save coords
fwrite(coords_global, file.path(out_path, "coords.csv"))


### Dynamic DBSCAN ###

# Rerun DBSCAN at the computed minPts values
for (i in seq_along(data.list)) {
  data.list[[i]][, cluster2 := mydbscan(data.list[[i]], eps, mranges[i,i2])]
  data.list[[i]][, cluster1 := mydbscan(data.list[[i]], eps, mranges[i,i1])]
}

# Assign the centroid via a weighted mean
coords_dynamic <- imap(data.list, function(dl, cb) {
  ret <- dl[,.(cb=cb, umi=sum(umi), beads=.N)]
  
  # Position the cell, using the highest minPts that produces DBSCAN=1
  if (max(dl$cluster1) == 1) {
    s <- dl[cluster1 == 1, .(x=weighted.mean(x, umi),
                             y=weighted.mean(y, umi),
                             umi1s=sum(umi),
                             beads1s=.N,
                             h1s=h_index(umi))]
    
  } else {
    s <- data.table(x=NA_real_, y=NA_real_, umi1s=NA_real_, beads1s=NA_integer_, h1s=NA_integer_)
  }
  ret[, names(s) := s]
  
  # Compute the score, using the highest minPts that produces DBSCAN=2
  if (max(dl$cluster2) >= 2) {
    s <- dl[, .(x=weighted.mean(x, umi),
                y=weighted.mean(y, umi),
                umi=sum(umi),
                beads=.N,
                h=h_index(umi)
               ), cluster2][cluster2 > 0][order(-umi, cluster2)]
    setcolorder(s, c("x", "y", "umi", "beads", "h", "cluster2"))
    ret[, c("x1","y1","umi1","beads1","h1","cluster1") := s[1]] # Highest UMI DBSCAN cluster
    ret[, c("x2","y2","umi2","beads2","h2","cluster2") := s[2]] # Second-highest UMI DBSCAN cluster
  } else if (max(dl$cluster2) == 1) {
      s <- dl[cluster2 == 1, .(x1=weighted.mean(x, umi),
                               y1=weighted.mean(y, umi),
                               umi1=sum(umi),
                               beads1=.N,
                               h1=h_index(umi),
                               cluster1=1L)]
      ret[, names(s) := s]
      s <- data.table(x2=NA_real_, y2=NA_real_, umi2=NA_real_, beads2=NA_integer_, h2=NA_integer_, cluster2=NA_integer_)
      ret[, names(s) := s]
  } else {
    s <- data.table(x1=NA_real_, y1=NA_real_, umi1=NA_real_, beads1=NA_integer_, h1=NA_integer_, cluster1=NA_integer_)
    ret[, names(s) := s]
    s <- data.table(x2=NA_real_, y2=NA_real_, umi2=NA_real_, beads2=NA_integer_, h2=NA_integer_, cluster2=NA_integer_)
    ret[, names(s) := s]
  }
  
}) %>% rbindlist

coords_dynamic[, c("eps", "minPts2", "minPts1") := data.table(eps=eps,
                                                              minPts2=mranges$i2*ms,
                                                              minPts1=mranges$i1*ms)]

# Plots
plot_ddbscan_xy(coords_dynamic) %>% make.pdf(file.path(out_path, "DDBSCANxy.pdf"), 7, 8)

# Save coords
fwrite(coords_dynamic, file.path(out_path, "coords2.csv"))


### Cell plots ###
plot_global_cellplots(data.list) %>% make.pdf(file.path(out_path, "DBSCAN.pdf"), 7, 8)
#dynamic_plots <- plot_dynamic_cellplots(data.list, coords_dynamic)

### Final check ###
stopifnot(coords_global$cb_index == coords_dynamic$cb_index)
stopifnot(file.exists(file.path(out_path, c("coords.csv", "coords2.csv"))))
