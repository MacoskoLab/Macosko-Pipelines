################################################################################
### Input: SBcounts.h5 (output of spatial-count.jl)
### Input: cb_whitelist.txt (newline-delimited list of cell barcodes to position)
### Output: coords.csv (+metadata)
################################################################################

source("plots.R")

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

# Print arguments
print("")

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
invisible(gc())

# Page 6
plot_SBmetrics(metadata) %>% make.pdf(file.path(out_path, "SBmetrics.pdf"), 7, 8)

# Write intermediate output
print("Writing intermediate matrix")
# order[]
fwrite(dt, file.path(out_path, "matrix.csv.gz"), quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE, compress="gzip")
metadata %>% map(as.list) %>% jsonlite::toJSON(pretty=T) %>% writeLines(file.path(out_path, "spatial_metadata.json"))
dt[, sb := NULL]

stopifnot(dt$cb %in% cb_whitelist)
stopifnot(dt$umi > 0)
stopifnot(!any(is.na(dt)))
stopifnot(sort(names(dt)) == c("cb", "umi", "x", "y"))

### DBSCAN #####################################################################

xlims = range(dt$x) ; ylims = range(dt$y)
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

# Run DBSCAN optimization
ms <- 1L # minPts step size
eps = epsilon * 15 # TODO: make this a vector, or pass-in-able
mranges <- furrr::future_map(data.list, function(dl){
  # base case
  v2 <- mydbscan(dl, eps, ms)
  if (max(v2) == 0) {return(list(i2=1L, c2=v2, i1=0L, c1=v2))}
  if (max(v2) == 1) {i2=1L ; c2=v2}
  
  # loop
  i <- 1L
  while (max(v2) > 0) {
    i <- i + 1L
    v1 <- v2
    v2 <- mydbscan(dl, eps, i*ms)
    if (max(v1) > 1 & max(v2) <= 1) {i2=i-1L ; c2=v1}
  }
  
  # i2/c2 is the index/clusters at the highest minPts that produces DBSCAN=2+ (if applicable)
  # i1/c1 is the index/clusters at the lowest minPts that produces DBSCAN=0, minus 1 step
  return(list(i2=i2, c2=c2, i1=i-1L, c1=v1))
}, .options = myoptions) %>% rbindlist

stopifnot(len(res)==len(data.list), names(res)==names(data.list))

# Consolidate results
mranges <- data.table(start = map_int(res, ~.$i2+1L),
                      end = map_int(res, ~.$i1))
for (i in seq_along(data.list)) {
  data.list[[i]][, clusters2 := mydbscan(data.list[[i]], eps, res[[i]]$i2) %T>% {stopifnot(.==res[[i]]$c2)}]
  data.list[[i]][, clusters1 := mydbscan(data.list[[i]], eps, res[[i]]$i1) %T>% {stopifnot(.==res[[i]]$c1)}]
}
#rm(res) ; invisible(gc())

### Global DBSCAN ###

# Compute minPts that places the greatest number of cells
i_opt <- IRanges::IRanges(start=mranges$start, end=mranges$end) %>% 
  IRanges::coverage() %>% as.integer %>% 
  {which(.==max(.))} %>% tail(1)
minPts <- i_opt*ms
print(g("Global optimum minPts: {minPts}"))

# Rerun DBSCAN with the optimal minPts
for (i in seq_along(data.list)) {
  data.list[[i]][, clusters := mydbscan(data.list[[i]], eps, minPts)]
}

# Assign the centroid via a weighted mean
coords_global <- imap(data.list, function(dl, cb){
  ret <- dl[,.(cb=cb,
               x=NA_real_, y=NA_real_,
               rmsd=NA_real_,
               clusters=max(clusters),
               umi=sum(umi), sumi=NA_real_,
               beads=.N, sbeads=NA_integer_,
               h=NA_integer_,
               eps=eps, minPts=minPts)]
  if (max(dl[,clusters]) == 1) {
    ret1 <- dl[clusters==1, .(x=weighted.mean(x, umi),
                              y=weighted.mean(y, umi),
                              rmsd=sqrt(weighted.mean((x-weighted.mean(x,umi))^2+(y-weighted.mean(y,umi))^2, umi)),
                              sumi=sum(umi),
                              sbeads=.N,
                              h=h_index(umi))]
    ret[, names(ret1) := ret1]
  }
  return(ret)
}) %>% rbindlist
print(g("Placed: {round(coords_global[,sum(clusters==1)/.N]*100, 2)}%"))

# Plots
plot_gdbscan_opt(coords_global, mranges) %>% make.pdf(file.path(out_path, "GDBSCANopt.pdf"), 7, 8)
plot_gdbscan_1(coords_global) %>% make.pdf(file.path(out_path, "GDBSCAN1.pdf"), 7, 8)

# Save coords
fwrite(coords_global, file.path(out_path, "coords_global.csv"))


### Dynamic DBSCAN ###
coords_dynamic <- pmap(list(data.list, names(data.list), res), function(dl, cb, r) {
  ret <- dl[,.(cb=cb,
               umi=sum(umi),
               beads=.N,
               x=NA_real_,  y=NA_real_,  sumi=NA_real_,  sbeads=NA_integer_, h=NA_integer_,
               x1=NA_real_, y1=NA_real_, sumi1=NA_real_, sbeads1=NA_integer_, h1=NA_integer_, cluster1=NA_integer_,
               x2=NA_real_, y2=NA_real_, sumi2=NA_real_, sbeads2=NA_integer_, h2=NA_integer_, cluster2=NA_integer_)]
  
  # Position the cell, using the highest minPts that produces DBSCAN=1
  if (max(r$c1) == 1) {
    s <- dl[r$c1 == 1, .(x=weighted.mean(x, umi),
                         y=weighted.mean(y, umi),
                         sumi=sum(umi),
                         sbeads=.N,
                         h=h_index(umi))]
    ret[, names(s) := s]
  }
  
  # Compute the score, using the highest minPts that produces DBSCAN=2
  if (max(r$c2) >= 2) {
    s <- dl[, .(x=weighted.mean(x, umi),
                y=weighted.mean(y, umi),
                sumi=sum(umi),
                sbeads=.N,
                h=h_index(umi)
               ), r$c2][r > 0][order(-sumi, r)]
    setcolorder(s, c("x", "y", "sumi", "sbeads", "h", "r"))
    ret[, c("x1","y1","sumi1","sbeads1","h1","cluster1") := s[1]] # Highest UMI DBSCAN cluster
    ret[, c("x2","y2","sumi2","sbeads2","h2","cluster2") := s[2]] # Second-highest UMI DBSCAN cluster
  } else if (max(r$c2) == 1) {
      s <- dl[r$c2 == 1, .(x1=weighted.mean(x, umi),
                           y1=weighted.mean(y, umi),
                           sumi1=sum(umi),
                           sbeads1=.N,
                           h1=h_index(umi),
                           cluster1=1L)]
      ret[, names(s) := s]
  }
  
  ret[, c("eps", "minPts2", "minPts1") := data.table(eps=eps,
                                                     minPts2=r$i2*ms,
                                                     minPts1=r$i1*ms)]
  
}) %>% rbindlist


# Plots
plot_ddbscan_xy(coords_dynamic) %>% make.pdf(file.path(out_path, "DDBSCANxy.pdf"), 7, 8)

# Save coords
fwrite(coords_dynamic, file.path(out_path, "coords_dynamic.csv"))


### Cell plots ###
plot_global_cellplots(data.list)

plot_dynamic_cellplots

dynamic_cellplot <- function(CB) {
  subdf <- data.list[[CB]][!is.na(x) & !is.na(y)][order(umi)]
  #subdf1 <- subdf[clusters==1]
  #subdf2 <- subdf[clusters==2]
  
}

# plots <- sample_bead_plots(data.list, coords)
# make.pdf(plots, file.path(out_path, "beadplots.pdf"), 7, 8)

cellbeadplot <- function(data.list, res, coords, CB) {
  stopifnot(CB %in% names(data.list))
  dl <- data.list[[CB]]
  
  stopifnot(CB %in% names(res))
  r <- res[[CB]]
  
  stopifnot(CB %in% coords$cb)
  row <- coords[cb == CB]
  
  ggplot() + geom_point(data=puckdf, mapping=aes(x=x,y=y), shape=1, col="grey") + 
             geom_point(data=dl[order(umi)], mapping=aes(x=x,y=y,col=umi)) + 
             geom_point(data=dl[clusters2==1], mapping=aes(x=x,y=y), col="red", shape=1, stroke=0.2, size=1.5) +
             geom_point(data=dl[clusters2==2], mapping=aes(x=x,y=y), col="green", shape=1, stroke=0.2, size=1.5) + 
             geom_point(data=row, mapping=aes(x=x1,y=y1),col="red",shape=0) + 
             geom_point(data=row, mapping=aes(x=x2,y=y2),col="green",shape=0) + 
             geom_point(data=row, mapping=aes(x=x,y=y),col="red",shape=5) +
    xlim(row[,c(min(x1,x2)-2*abs(x1-x2), max(x1,x2)+2*abs(x1-x2))]) + 
    ylim(row[,c(min(y1,y2)-2*abs(y1-y2), max(y1,y2)+2*abs(y1-y2))]) +
             coord_fixed(ratio=1)
  plot_grid(p, p)
  
}

### Final check ###
stopifnot(coords_global$cb_index == coords_dynamic$cb_index)
stopifnot(file.exists(file.path(out_path,c("coords_global.csv", "coords_dynamic.csv"))))
