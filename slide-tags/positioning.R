################################################################################
### Input: SBcounts.h5 (output of spatial-count.jl)
### Input: cb_whitelist.txt (newline-delimited list of cell barcodes to position)
### Output: coords.csv (placements+metadata)
################################################################################

suppressMessages(source("helpers.R"))
suppressMessages(source("plots.R"))

# Load arguments
library(optparse)
arguments <- OptionParser(
  usage = "Usage: Rscript positioning.R SBcounts_path CBwhitelist_path out_path [options]",
  option_list = list(
    make_option("--minPts", type="integer", help = "dbscan optimization override"),
    make_option("--eps", type="integer", help = "dbscan optimization override"),
    make_option("--knn", type="integer", default=51L, help = "Number of bead neighbors used to compute eps [default: %default]"),
    make_option("--cmes", type="double", default=0.0, help = "Reconstruction parameter"),
    make_option("--prob", type="double", default=1.0, help = "Proportion of reads to retain [default: 1.0]"),
    make_option("--cores", type="integer", default=-1L, help = "The number of parallel processes to use [default: -1]")
  )
) %>% parse_args(positional_arguments=3)

sb_path <- arguments$args[[1]]  ; print(g("sb_path: {sb_path}"))
cb_path <- arguments$args[[2]]  ; print(g("cb_path: {cb_path}"))
out_path <- arguments$args[[3]] ; print(g("out_path: {out_path}"))

minPts <- arguments$options$minPts ; print(g("minPts: {minPts}"))
eps <- arguments$options$eps       ; print(g("eps: {eps}"))
knn <- arguments$options$knn       ; print(g("knn: {knn}"))
cmes <- arguments$options$cmes     ; print(g("cmes: {cmes}"))
prob <- arguments$options$prob     ; print(g("prob: {prob}"))
cores <- arguments$options$cores %>% ifelse(.<1, parallelly::availableCores(), .) ; print(g("cores: {cores}"))
setDTthreads(cores)

rm(arguments)

# Check input arguments
stopifnot(filter(rhdf5::h5ls(sb_path), group=="/matrix")$name == c("cb_index", "reads", "sb_index", "umi"))
stopifnot(file.exists(cb_path))
if (!dir.exists(out_path)) {dir.create(out_path, recursive=TRUE)}
stopifnot(dir.exists(out_path))

# Load the spatial barcode count matrix
f <- function(p){return(rhdf5::h5read(sb_path, p))}
dt <- ReadSpatialMatrix(f)
metadata <- ReadSpatialMetadata(f)
print(g("{add.commas(sum(dt$reads))} spatial barcode reads loaded"))
stopifnot(names(dt) == c("cb","umi","sb","reads"))

if (prob < 1) {
  dt[, reads := rbinom(.N, reads, prob)]
  dt <- dt[reads > 0]
  print(g("{add.commas(sum(dt$reads))} downsampled spatial barcode reads"))
}

# load the CB whitelist
cb_whitelist <- readLines(cb_path)

# determine CB remap
remap <- determine_remap_10X_CB(cb_whitelist, dt)
if (remap) { cb_whitelist %<>% remap_10X_CB }
metadata$SB_info$remap_10X_CB <- remap

# validate CB whitelist
stopifnot(class(cb_whitelist) == "character")
stopifnot(!duplicated(cb_whitelist))
stopifnot(uniqueN(nchar(cb_whitelist)) == 1)
stopifnot(map_lgl(strsplit(cb_whitelist, ""), ~all(. %in% c("A","C","G","T"))))
print(g("{len(cb_whitelist)} cell barcodes loaded"))
invisible(gc())

### Fuzzy matching #############################################################

print("Performing fuzzy cell-barcode matching")
setnames(dt, "cb", "cb_fuzzy")
# Add HD1 matches
df <- data.table(cb=cb_whitelist)
df <- df[, .(cb_fuzzy=listHD1neighbors(cb), match="fuzzy"), cb
         ][!cb_fuzzy %in% cb_whitelist
           ][cb_fuzzy %in% dt$cb_fuzzy]
# Label ambiguous matches
df[, `:=`(cb = ifelse(.N > 1, NA, cb),
          match = ifelse(.N > 1, "ambig", match)), cb_fuzzy]
df %<>% unique
# Add exact matches
df <- rbindlist(list(df, data.table(cb=cb_whitelist,
                                    cb_fuzzy=cb_whitelist,
                                    match="exact")[cb_fuzzy %in% dt$cb_fuzzy]))
# Check + factorize
stopifnot(!any(duplicated(df[,.(cb_fuzzy, cb)])))
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
dt <- merge(dt, df, by = "cb_fuzzy", all.x = TRUE, all.y = FALSE, sort = FALSE)
stopifnot(is.null(key(dt)), is.null(indices(dt)))

# Record metadata
metadata$CB_matching <- dt[, .(reads=sum(reads)), match][order(-reads)] %>% 
  {setNames(.$reads, .$match %>% as.character %>% replace_na("none"))}
stopifnot(sum(dt$reads) == sum(metadata$CB_matching))

# Clean up
dt[, match := NULL]
setnames(dt, "cb_fuzzy", "cr")
setcolorder(dt, c("cr", "cb", "umi", "sb", "reads"))
rm(df) ; invisible(gc())

# Save matching result
# fwrite(dt, file.path(out_path, "matrix0.csv.gz"), quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE, compress="gzip")

print(g("{dt[!is.na(cb), .N] %>% add.commas} matched UMIs"))

### Load puck ##################################################################

# Load the puck
print("Loading the puck")
res <- ReadPuck(f)
puckdf <- res[[1]]
metadata$puck_info %<>% c(res[[2]])
rm(res) ; invisible(gc())

# Compute eps scale using the kth neighbor
if (!exists("eps") || !is.numeric(eps) || !(eps > 0)) {
  eps <- RANN::nn2(data = puckdf[, .(x,y)],
                   query = puckdf[sample(.N, 10000), .(x,y)],
                   k = knn)$nn.dists[,knn] %>% median
  print(g("eps: {eps}"))
} else {
  print(g("Override eps: {eps}"))
}

# Filter reads with a low-quality spatial barcode
print("Removing low-quality spatial barcodes")
dt[, m := sb %in% puckdf$sb]
metadata$SB_filtering %<>% c(reads_lqsb=dt[m == FALSE, sum(reads)])
dt <- dt[m == TRUE]
dt[, m := NULL]
invisible(gc())

# Page 4
plot_SBlibrary(dt, f) %>% make.pdf(file.path(out_path, "SBlibrary.pdf"), 7, 8)

# Page 5
plot_SBplot(dt, puckdf) %>% make.pdf(file.path(out_path, "SBplot.pdf"), 7, 8)

# Delete cell barcodes for cells that were not called
print("Removing non-whitelist cells")
metadata$SB_filtering %<>% c(reads_nocb = dt[is.na(cb), sum(reads)])
metadata$SB_info$UMI_pct_in_called_cells <- round(dt[,sum(!is.na(cb))/.N]*100, digits=2) %>% paste0("%")
metadata$SB_info$sequencing_saturation <- round((1 - nrow(dt) / sum(dt$reads)) * 100, 2) %>% paste0("%")
dt[, cr := NULL]
dt <- dt[!is.na(cb)]
dt <- dt[, .(reads=sum(reads)), .(cb,umi,sb)] # collapse post-matched barcodes
stopifnot(nrow(dt) == nrow(unique(dt[,.(cb,umi,sb)])))
invisible(gc())

# Remove chimeric UMIs
print("Removing chimeras")
dt[, m := reads == max(reads) & sum(reads == max(reads)) == 1, .(cb, umi)]
metadata$SB_filtering %<>% c(reads_chimeric = dt[m==FALSE, sum(reads)])
dt <- dt[m == TRUE]
dt[, m := NULL]
invisible(gc())

# Compute sequencing saturation per cell
metadata$sequencing_saturation <- dt[, .(ss=1-.N/sum(reads)), cb] %>% {setNames(.[,ss], .[,cb])}

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

# Write intermediate results
print("Writing intermediates")
fwrite(dt, file.path(out_path, "matrix.csv.gz"), quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE, compress="gzip")
metadata %>% map(as.list) %>% jsonlite::toJSON(pretty=TRUE) %>% writeLines(file.path(out_path, "spatial_metadata.json"))
dt[, sb := NULL]

# Intermediate checks
stopifnot(dt$cb %in% cb_whitelist)
stopifnot(levels(dt$cb) == cb_whitelist)
stopifnot(dt$umi > 0)
stopifnot(!any(is.na(dt)))
stopifnot(sort(names(dt)) == c("cb", "umi", "x", "y"))

### DBSCAN #####################################################################

data.list <- split(dt, by="cb", drop=FALSE, keep.by=FALSE)
stopifnot(sort(names(data.list)) == sort(cb_whitelist))
data.list <- data.list[cb_whitelist]
stopifnot(names(data.list) == cb_whitelist)
rm(dt) ; invisible(gc())
print(g("Running positioning on {len(data.list)} cells"))

# Prepare for positioning
library(future)
library(furrr)
options(future.globals.maxSize = 1024 * 1024 * 1024)
myoptions <- furrr::furrr_options(packages=c("data.table"), seed=TRUE, scheduling=1)
future::plan(future::multisession, workers=cores)
mydbscan <- function(dl, eps, minPts) {
  if (nrow(dl) == 0) {return(numeric(0))}
  dbscan::dbscan(x = dl[,.(x,y)],
                 eps = eps,
                 minPts = minPts,
                 weights = dl$umi,
                 borderPoints = TRUE)$cluster
}
ipercluster <- function(dl, v) {
  dls <- sort(unique(v[v>0])) %>% map(~dl[v==.])
  map_int(dls, function(sdl){
    dbscan::frNN(sdl[,.(x,y)], eps, sort=FALSE)$id %>% imap_dbl(~sdl[.x,sum(umi)]+sdl[.y,umi]) %>% max %>% divide_by_int(ms)
  })
}
iestimate <- function(dl) {
  dbscan::frNN(dl[,.(x,y)], eps, query=dl[umi==max(umi),.(x,y)], sort=FALSE)$id %>% 
    map_dbl(~dl[.,sum(umi)]) %>% max %>% divide_by_int(ms)
}
centroid_dists <- function(dl, v) {
  unique(v[v>0]) %>% {stopifnot(sort(.) == seq_along(.))}
  if (max(v) < 2) {return(0)}
  centroids <- dl[v>0, .(x=weighted.mean(x,umi),
                         y=weighted.mean(y,umi)), v[v>0]][order(v)]
  cdist(centroids[,.(x,y)], centroids[v==1,.(x,y)]) %>% as.numeric
}

### Run DBSCAN optimization ###
ms <- 1L # minPts step size

# Find DBSCAN 1->0 boundary
i1 <- future_map_int(data.list, function(dl) {
  if (nrow(dl) == 0) {return(0L)}
  if (nrow(dl) == 1) {return(sum(dl[,umi]))}
  
  # minPts estimation
  i <- iestimate(dl)
  v <- mydbscan(dl, eps, i * ms) %T>% {stopifnot(max(.)>=1)}
  
  # minPts estimation improvement
  i <- ipercluster(dl, v) %>% max %T>% {stopifnot(.>=i)}
  v <- mydbscan(dl, eps, i * ms) %T>% {stopifnot(max(.)>=1)}
  
  # Increase minPts until DBSCAN=0
  while (max(v) > 0) {
    i <- i + 1L
    v <- mydbscan(dl, eps, i * ms)
  }
  return(i-1L)
}, .options = myoptions)

# Find DBSCAN 2->1 boundary
i2 <- future_map_int(data.list, function(dl) {
  if (nrow(dl) == 0) {return(0L)}
  if (nrow(dl) == 1) {return(0L)}
  
  # minPts estimation
  i <- iestimate(dl)
  v <- mydbscan(dl, eps, i * ms) %T>% {stopifnot(max(.)>=1)}
  
  ii <- integer(0)
  repeat {
    ii %<>% c(i)
    
    if (any(v!=1)) {
      i <- iestimate(dl[v!=1]) # TODO
      v <- mydbscan(dl, eps, i * ms)
    } else {
      return(0L)
    }
    
    if (max(v) >= 2) {
      i <- ipercluster(dl, v) %>% sort %>% pluck(-2)
      v <- mydbscan(dl, eps, i * ms)
    }
    
    if (max(v) >= 2 || i %in% ii) {break}
  }
  
  # Decrease minPts until DBSCAN=2
  if (max(v) < 2) {
    # TODO: stopifnot(FALSE)
    while (max(v) < 2 && i >= 1L) {
      i <- i - 1L
      v <- mydbscan(dl, eps, i * ms)
    }
    return(i)
  }
  
  # Increase minPts until DBSCAN=1/0
  while (max(v) >= 2) {
    i <- i + 1L
    v <- mydbscan(dl, eps, i * ms)
  }
  return(i-1L)
}, .options = myoptions)

stopifnot(names(i1)==cb_whitelist)
stopifnot(names(i2)==cb_whitelist)
mranges <- data.table(i2=i2, i1=i1)
stopifnot(mranges$i2 <= mranges$i1)
rm(i1, i2) ; invisible(gc())

# Decrease 2->1 boundary if all centroids are within [eps * cmes] of DBSCAN=1
if (cmes > 0) {
  mranges$i2 <- future_map2_int(data.list, mranges$i2, function(dl, i2) {
    v <- mydbscan(dl, eps, i2 * ms)
    while (i2 >= 1L && max(centroid_dists(dl,v)) < eps * cmes) {
      i2 <- i2 - 1
      v <- mydbscan(dl, eps, i2 * ms)
    }
    return(i2)
  }, .options = myoptions)
}

### Global DBSCAN ###

# Compute minPts that places the greatest number of cells
if (exists("minPts") && is.numeric(minPts) && minPts >= 0) {
  print(g("Override minPts: {minPts}"))
} else {
  i_opt <- IRanges::IRanges(start=mranges$i2+1L, end=mranges$i1) %>% 
    IRanges::coverage() %>% as.integer %>% 
    {which(.==max(.))} %>% tail(1)
  minPts <- i_opt * ms
  print(g("Global optimum minPts: {minPts}"))
}

# Rerun DBSCAN with the optimal minPts
for (i in seq_along(data.list)) {
  data.list[[i]][, cluster := mydbscan(data.list[[i]], eps, minPts)]
}

# Merge clusters within [eps * cmes] of DBSCAN=1
if (cmes > 0) {
  for (i in seq_along(data.list)) {
    if (nrow(data.list[[i]]) > 0) {
      near_centroids <- which(centroid_dists(data.list[[i]], data.list[[i]][,cluster]) < eps * cmes)
      data.list[[i]][cluster %in% near_centroids, cluster := 1]
      data.list[[i]][cluster > 0, cluster := match(cluster, sort(unique(cluster)))]
    }
  }
}

# Assign the centroid via a weighted mean
coords <- imap(data.list, function(dl, cb){
  ret <- dl[,.(cb=cb,
               umi=sum(umi),
               beads=.N,
               max=ifelse(.N>0, max(umi), 0),
               clusters=ifelse(.N>0, max(cluster), 0))]
  
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
  
}) %>% rbindlist(use.names=TRUE, fill=TRUE)
coords[, `:=`(eps=eps, minPts=minPts)]
coords[clusters==1, `:=`(x=x1, y=y1)]
print(g("Placed: {round(coords[,sum(clusters==1)/.N]*100, 2)}%"))

# Final check
stopifnot(names(data.list) == cb_whitelist)
stopifnot(coords$cb == cb_whitelist)
stopifnot(nrow(mranges) == len(cb_whitelist))

# Write coords
setcolorder(coords, c("cb","x","y"))
if (remap) { coords$cb %<>% remap_10X_CB }
fwrite(coords, file.path(out_path, "coords.csv"))

# Plots
plot_dbscan_opt(coords, mranges) %>% make.pdf(file.path(out_path, "DBSCANopt.pdf"), 7, 8)
plot_dbscan_1(coords) %>% make.pdf(file.path(out_path, "DBSCAN1.pdf"), 7, 8)
plot_dbscan_2(coords, cmes) %>% make.pdf(file.path(out_path, "DBSCAN2.pdf"), 7, 8)

### Dynamic DBSCAN ###

# Rerun DBSCAN at the computed minPts values
for (i in seq_along(data.list)) {
  data.list[[i]][, cluster2 := mydbscan(data.list[[i]], eps, mranges[i,i2]*ms)]
  data.list[[i]][, cluster1 := mydbscan(data.list[[i]], eps, mranges[i,i1]*ms)]
}

# Merge clusters within [eps * cmes] of DBSCAN=1
if (cmes > 0) {
  for (i in seq_along(data.list)) {
    if (nrow(data.list[[i]]) > 0) {
      near_centroids <- which(centroid_dists(data.list[[i]], data.list[[i]][,cluster2]) < eps * cmes)
      data.list[[i]][cluster2 %in% near_centroids, cluster2 := 1]
      data.list[[i]][cluster2 > 0, cluster2 := match(cluster2, sort(unique(cluster2)))]
    }
  }
}

# Assign the centroid via a weighted mean
coords2 <- imap(data.list, function(dl, cb) {
  ret <- dl[,.(cb=cb, umi=sum(umi))]
  
  # Position the cell, using the highest minPts that produces DBSCAN=1
  if (nrow(dl) > 0 && dl[,max(cluster1)] == 1) {
    s <- dl[cluster1 == 1, .(x=weighted.mean(x, umi),
                             y=weighted.mean(y, umi))]
    ret[, names(s) := s]
  }
  
  # Compute statistics, using the highest minPts that produces DBSCAN=2
  s <- dl[cluster2 > 0, .(x=weighted.mean(x, umi),
                          y=weighted.mean(y, umi),
                          umi=sum(umi),
                          beads=.N,
                          h=h_index(umi)
  ), cluster2][order(-umi, cluster2)]
  setcolorder(s, c("x", "y", "umi", "beads", "h", "cluster2"))
  ret[, c("x1","y1","umi1","beads1","h1","cluster1") := s[1]] # Highest UMI DBSCAN cluster
  ret[, c("x2","y2","umi2","beads2","h2","cluster2") := s[2]] # Second-highest UMI DBSCAN cluster
  
  return(ret)
  
}) %>% rbindlist(use.names=TRUE, fill=TRUE)

# Compute the score
F1 <- ecdf(coords2$umi1)
coords2[, score := 1-F1(umi2)/F1(umi1)]

# Add DBSCAN parameters
coords2[, c("eps", "minPts2", "minPts1") := data.table(eps=eps,
                                                       minPts2=mranges$i2*ms,
                                                       minPts1=mranges$i1*ms)]

# Save results
plot_dbscan_score(coords2) %>% make.pdf(file.path(out_path, "DBSCANscore.pdf"), 7, 8)
fwrite(coords2, file.path(out_path, "coords2.csv"))

### Final check ################################################################

print("Writing output...")

# Combine into one PDF
plotlist <- c("SBmetrics.pdf", "SBlibrary.pdf", "SBplot.pdf",
              "DBSCANopt.pdf", "DBSCAN1.pdf", "DBSCAN2.pdf", "DBSCANscore.pdf")
pdfs <- file.path(out_path, plotlist)
pdfs %<>% keep(file.exists)
qpdf::pdf_combine(input=pdfs, output=file.path(out_path, "SBsummary.pdf"))
file.remove(pdfs)

# Plot individual cell plots
c(plot_dbscan_cellplots(data.list),
  plot_debug_cellplots(data.list, coords2)) %>%
  make.pdf(file.path(out_path, "DBSCANs.pdf"), 7, 8)

# rbindlist(data.list, use.names=TRUE) %>% fwrite(file.path(out_path, "matrixd.csv.gz"), quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE, compress="gzip")

stopifnot(coords$cb == coords2$cb)
