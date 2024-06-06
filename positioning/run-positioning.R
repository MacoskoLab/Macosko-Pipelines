library(glue) ; g=glue ; len=length
library(matrixStats)
library(gridExtra)
library(magrittr)
library(ggplot2)
library(cowplot)
library(stringr)
library(Seurat)
library(dbscan)
library(future)
library(rlist)
library(dplyr)
library(purrr)
library(furrr)
library(rhdf5)
library(qpdf)
library(qs)

# setwd("~/spatial")

### Download files #############################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript spatial.R RNApath SBpath", call. = FALSE)
}

RNApath <- args[[1]] ; print(g("RNApath: {RNApath}"))
SBpath <- args[[2]] ; print(g("SBpath: {SBpath}"))

#RNApath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/03_COUNTS/240421_SL-EXB_0289_A22KHMWLT3/SI-TT-G8_04-09_Pu-GP_GEX_rxn16"
#SBpath = "gs://fc-fc3a2afa-6861-4da5-bb37-60ebb40f7e8d/04_SPATIAL/240421_SL-EXB_0289_A22KHMWLT3/SI-TT-G12"

stopifnot(!dir.exists("RNAcounts"))
stopifnot(!dir.exists("SBcounts"))

stopifnot(checkgsfile(RNApath))
stopifnot(checkgsfile(SBpath))

system("mkdir RNAcounts")
RNAtech <- download_RNA_data(RNApath)
stopifnot(length(list.files("RNAcounts")) >= 1)

system("mkdir SBcounts")
download_SB_data(SBpath)
stopifnot(length(list.files("SBcounts")) == 1)
f <- function(p){return(h5read("SBcounts/SBcounts.h5", p))}

system("mkdir plots")

### Load the seurat ############################################################
obj <- load_seurat("RNAcounts/filtered_feature_bc_matrix.h5", "RNAcounts/molecule_info.h5")

Misc(obj, "RNA_path") <- RNApath
Misc(obj, "SB_path") <- SBpath
Misc(obj, "method") <- RNAtech
rm(RNApath, SBpath, RNAtech)

cb_whitelist = unname(obj$cb)
stopifnot(!duplicated(cb_whitelist))
stopifnot(len(unique(nchar(cb_whitelist))) == 1)
stopifnot(map_lgl(strsplit(cb_whitelist, ""), ~all(. %in% c("A","C","G","T"))))
Misc(obj, "called_cells") <- len(cb_whitelist)
writeLines(cb_whitelist, "cb_whitelist.txt")
system("Rscript load_matrix.R SBcounts/SBcounts.h5 cb_whitelist.txt")

df = read.table("matrix.csv", header=T, sep=",")

### Positioning methods ########################################################

chunk_vector <- function(v, chunk_size) {return(split(v, ceiling(seq_along(v) / chunk_size)))}

# Do a grid search to find the ideal DBSCAN parameters
ncores = 20L ; plan(multisession, workers=ncores)
opt_dbscan <- function(data.list) {
  eps.vec = c(50) ; minPts.vec = c(3:42)
  res = data.frame() ; i = 0
  repeat{
    params = expand.grid(eps.vec,minPts.vec) %>% setNames(c("eps","minPts"))
    row_lists = chunk_vector(1:nrow(params), round(nrow(params)/ncores))
    
    params$pct = furrr::future_map(row_lists, function(v) {
      map_dbl(v, function(i) {
        m = map_lgl(data.list, ~max(dbscan::dbscan(.[c("x_um","y_um")], eps=params$eps[[i]], minPts=params$minPts[[i]], weights=.$umi)$cluster) == 1)
        return(sum(m)/length(m))
      })
    }, .options=furrr_options(seed=T)) %>% flatten_dbl
    
    res = rbind(res, params)
    
    if (which.max(res$pct)<0.9*nrow(res) || i >= 26) {
      break
    }
    
    minPts.vec = minPts.vec + 40
    i = i + 1
  }
  
  params = res ; rm(res)
  params$is.max = params$pct==max(params$pct)
  eps = params$eps[params$is.max][[1]] ; minPts = params$minPts[params$is.max][[1]]
  pct.placed = round(max(params$pct)*100,2)
  print(g("Optimal eps: {eps}    Optimal minPts: {minPts}    %placed: {pct.placed}"))
  
  return(c(eps,minPts,pct.placed))
}

# Add the DBSCAN clusters to the dataframes
run_dbscan <- function(data.list, eps, minPts) {
  lapply(data.list, function(df){
    df$cluster <- dbscan::dbscan(df[c("x_um","y_um")], eps=eps, minPts=minPts, weights=df$umi)$cluster
    return(df)
  })
}

# assign centroid and record metadata
create_coords <- function(data.list) {
  coords = lapply(data.list, function(df) {
    p = c(x_um=NA,
          y_um=NA,
          DBSCAN_clusters=max(df$cluster),
          SB_umi = sum(df$umi),
          SNR=NA,
          SB_bin = unique(df$bin) %>% {ifelse(is.null(.), NA, .)},
          minPts = unique(df$minPts) %>% {ifelse(is.null(.), NA, .)},
          eps = unique(df$eps) %>% {ifelse(is.null(.), NA, .)},
          pct.placed = unique(df$pct.placed) %>% {ifelse(is.null(.), NA, .)})
    if (max(df$cluster) == 1) {
      sdf = dplyr::filter(df, cluster==1)
      p[["x_um"]] = matrixStats::weightedMedian(sdf$x_um,w=sdf$umi)
      p[["y_um"]] = matrixStats::weightedMedian(sdf$y_um,w=sdf$umi)
      p[["SNR"]] = sum(sdf$umi)/sum(df$umi)
    }
    return(p)
  }) %>% bind_rows %>% as.data.frame %>% mutate(cb_index=as.numeric(names(data.list))) %>% select(cb_index, everything())
  return(coords)
}

# Run these methods on the entire data
normal_positioning <- function(df) {
  stopifnot("umi" %in% colnames(df)) # if this fails, you haven't grouped by cb,sb and counted umis
  data.list = split(df, df$cb_index)
  params = opt_dbscan(data.list)
  data.list %<>% run_dbscan(eps=params[[1]], minPts=params[[2]])
  data.list %<>% map(~mutate(.,eps=params[[1]],minPts=params[[2]],pct.placed=params[[3]]))
  coords <- create_coords(data.list)
  return(list(coords,data.list))
}

# Split the data into 10 SB UMI buckets and run on each
binned_positioning <- function(df) {
  stopifnot("umi" %in% colnames(df)) # if this fails, you haven't grouped by cb,sb and counted umis
  data.list = split(df, df$cb_index)
  
  # create the deciles
  quants = c(0,quantile(map_dbl(data.list,~sum(.$umi)), probs = seq(0.1, 1, by = .1))) %>% unname
  paste("Deciles: ", paste(round(quants),collapse=" "))
  umicounts = map(data.list,~sum(.$umi))
  data.lists=map(1:(len(quants)-1),~data.list[umicounts>quants[.] & umicounts<=quants[.+1]])
  
  # run positioning on each decile
  data.list = map2(data.lists,quants[-1],function(data.list,quant) {
    if (len(data.list) == 0) {print(g({"skipping quantile, empty"}));return(list())}
    params = opt_dbscan(data.list)
    data.list %<>% run_dbscan(eps=params[[1]], minPts=params[[2]])
    data.list %<>% map(~mutate(.,bin=quant,eps=params[[1]],minPts=params[[2]],pct.placed=params[[3]]))
    return(data.list)
  }) %>% list_flatten
  
  stopifnot(len(data.list)==len(unique(df$cb_index)))
  
  coords <- create_coords(data.list)
  return(list(coords,data.list))
}

### Positioning ################################################################

# add spatial positions from puck
df = merge(x=df, y=puckdf, all.x=T, by="sb_index")
df %<>% filter(!is.na(x_um),!is.na(y_um))

# run binned positioning (see method above)
res = normal_positioning(df)
coords = res[[1]] ; data.list = res[[2]]
rm(res) ; invisible(gc())

# merge with seurat object
coords %<>% mutate(cb = cb_whitelist[cb_index])
rownames(coords) = paste0(coords$cb,"-1")
obj = AddMetaData(obj,coords)
Misc(obj,"pct.placed") = round(sum(!is.na(obj$x_um))/ncol(obj)*100,2)
emb = obj@meta.data[,c("x_um","y_um")] ; colnames(emb) = c("s_1","s_2")
obj[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_", assay = "RNA")

coords %>% select(cb,everything()) %>% write.csv("coords.csv",quote=F,row.names=F)
qsave(obj, "seurat.qs")




################################################################################
### CREATE PDF #################################################################
################################################################################

### Page 5: DBSCAN #############################################################

plot_dbscan <- function(obj, coords) {
  # Panel 1: DBSCAN cluster distribution
  d = data.frame(x=coords$DBSCAN_clusters) %>% rowwise %>% mutate(x=min(x,10)) %>% ungroup
  p1 = ggplot(d,aes(x=x)) + geom_histogram(aes(y = after_stat(count)/sum(after_stat(count))*100), binwidth=.5) +
    geom_text(aes(label = sprintf("%1.0f%%", after_stat(count)/sum(after_stat(count))*100), y=after_stat(count)/sum(after_stat(count))*100), stat="bin", binwidth=1, vjust=-0.5)+
    theme_classic() + xlab("Num DBSCAN clusters") + ylab("Percent") +
    scale_y_continuous(limits=c(0,100)) +
    scale_x_continuous(breaks=min(d$x):max(d$x), labels=(min(d$x):max(d$x)) %>% {ifelse(.==10, "10+", .)}) +
    ggtitle("DBSCAN cluster distribution")
  
  # Panel 2: SNR density
  max_density_x = density(obj$SNR %>% na.omit) %>% {.$x[which.max(.$y)]}
  max_density_x = median(obj$SNR, na.rm = T)
  p2 = obj@meta.data %>% filter(!is.na(x_um)) %>% ggplot(aes(x = SNR)) +
    geom_density() + 
    theme_minimal() +
    labs(title = "SNR per cell (density)", x = "SNR", y = "Density") + 
    geom_vline(xintercept = max_density_x, color = "red", linetype = "dashed") +
    annotate(geom = 'text', label = round(max_density_x, 2), x = max_density_x+0.01, y = Inf, hjust = 0, vjust = 1, col="red")
  
  # Panel 3: RNA umi vs SB umi
  p3 = data.frame(x=obj$nCount_RNA,y=obj$SB_umi,placed=!is.na(obj$x_um)) %>% {ggplot(.,aes(x=log10(x),y=log10(y),col=placed))+geom_point(size=0.2)+theme_bw()+xlab("RNA UMI")+ylab("SB UMI")+ggtitle("SB UMI vs. RNA UMI")+theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"), legend.background = element_blank(), legend.title=element_text(size=10), legend.text=element_text(size=8), legend.margin=margin(0,0,0,0,"pt"), legend.box.margin=margin(0,0,0,0,"pt"), legend.spacing.y = unit(0.1,"lines"), legend.key.size = unit(0.5, "lines"))}
  
  # Panel 4: DBSCAN parameters
  df = coords %>% select(SB_bin,minPts,eps,pct.placed) %>% distinct %>% arrange(SB_bin) %>% mutate(SB_bin=round(SB_bin,2),pct.placed=round(pct.placed,2) %>% paste0("%"))
  rownames(df) <- NULL
  p4 = plot_grid(gdraw("DBSCAN parameters"),plot.tab(df),ncol=1,rel_heights=c(1,17))
  
  plot = plot_grid(p1,p2,p3,p4,ncol=2)
  return(plot)
}
plot <- plot_dbscan(obj, coords)
make.pdf(plot,"plots/5DBSCAN.pdf",7,8)

### Page 6: Spatial ############################################################

p1 = DimPlot(obj[,!is.na(obj$x_um)],reduction="spatial")+coord_fixed(ratio=1)+
  ggtitle(g("%placed: {round(sum(!is.na(coords$x_um))/nrow(coords)*100,2)} ({sum(!is.na(obj$x_um))}/{ncol(obj)})")) + 
  NoLegend() + xlab("x-position (\u00B5m)") + ylab("y-position (\u00B5m)")
p2 = DimPlot(obj[,!is.na(obj$x_um)], reduction="spatial",split.by="seurat_clusters",ncol=5) + theme_void() + coord_fixed(ratio=1) + NoLegend()
plot = plot_grid(p1, p2, ncol=1, rel_heights=c(1,1))
make.pdf(plot,"plots/6spatial.pdf",7,8)

### Sample bead plots ##########################################################

sample_bead_plots <- function(data.list, xlim=range(df$x_um),ylim=range(df$y_um)) {
  plot.sb <- function(subdf) {
    subdf %<>% arrange(lr)
    subdf1 <- filter(subdf,cluster==1)
    subdf2 <- filter(subdf,cluster==2)
    ggplot()+coord_fixed(ratio=1,xlim=xlim,ylim=ylim)+theme_void()+
      geom_point(data=subdf, mapping=aes(x=x_um,y=y_um,col=lr),size=2,shape=16)+
      geom_point(aes(x=matrixStats::weightedMedian(subdf1$x_um,w=subdf1$umi),
                     y=matrixStats::weightedMedian(subdf1$y_um,w=subdf1$umi)),
                 color="red",shape=0,size=3) + 
      geom_point(aes(x=matrixStats::weightedMedian(subdf2$x_um,w=subdf2$umi),
                     y=matrixStats::weightedMedian(subdf2$y_um,w=subdf2$umi)),
                 color="green",shape=0,size=3) + ggtitle(g("{unique(subdf$cb_index)}")) +
      theme(legend.key.width=unit(0.5,"lines"), legend.position="right", legend.key.height=unit(1,"lines"), legend.title=element_blank(), legend.spacing.y=unit(0.2,"lines"), legend.margin=margin(0,0,0,0,"lines"), legend.box.margin=margin(0,0,0,0,"pt"), legend.box.background=element_blank(), legend.background=element_blank(), legend.direction="vertical", legend.justification="left",legend.box.just="left",legend.box.spacing=unit(0,"cm"))
    
  }
  
  list0 = data.list %>% keep(~max(.$cluster)==0) %>% map(~arrange(.,umi) %>% filter(!is.na(x_um), !is.na(y_um)))
  list1 = data.list %>% keep(~max(.$cluster)==1) %>% map(~arrange(.,umi) %>% filter(!is.na(x_um), !is.na(y_um)))
  list2 = data.list %>% keep(~max(.$cluster)==2) %>% map(~arrange(.,umi) %>% filter(!is.na(x_um), !is.na(y_um)))
  
  list0 = list0[as.integer(map(list0,nrow))!=0]
  list1 = list1[as.integer(map(list1,nrow))!=0]
  list2 = list2[as.integer(map(list2,nrow))!=0]
  
  if(len(list0) > 0) {
    p1 = plot_grid(ggdraw()+draw_label("DBSCAN=0"), map(sample(list0,min(12,len(list0)),replace=F),plot.sb) %>% {plot_grid(plotlist=.,ncol=3)}, ncol=1,rel_heights=c(0.1,2))
  } else {p1 = gdraw("No DBSCAN = 0")}
  if(len(list1) > 0) {
    p2 = plot_grid(ggdraw()+draw_label("DBSCAN=1"), map(sample(list1,min(12,len(list1)),replace=F),plot.sb) %>% {plot_grid(plotlist=.,ncol=3)}, ncol=1,rel_heights=c(0.1,2))
  } else {p2 = gdraw("No DBSCAN = 1")}
  if(len(list2) > 0) {
    p3 = plot_grid(ggdraw()+draw_label("DBSCAN=2"), map(sample(list2,min(12,len(list2)),replace=F),plot.sb) %>% {plot_grid(plotlist=.,ncol=3)}, ncol=1,rel_heights=c(0.1,2))
  } else {p3 = gdraw("No DBSCAN = 2")}
  
  plots = list(p1,p2,p3)
  
  return(plots)
}
plots <- sample_bead_plots(data.list)
make.pdf(plots,"plots/SB.pdf",7,7)

### Save output ################################################################

pdfs = c("0cellranger.pdf","1cellcalling.pdf", "2umap.pdf", "3rawspatial.pdf", "4beadplot.pdf", "5DBSCAN.pdf","6spatial.pdf","7metrics.pdf", "SB.pdf") %>% paste0("plots/",.)
qpdf::pdf_combine(input = pdfs, output = "summary.pdf")
qsave(obj, "seurat.qs") # we added some more metadata while plotting

print("Done!")

# Make some plots
obj %<>% plot0("RNAcounts/metrics_summary.csv")             # cell ranger output
obj %<>% plot1("RNAcounts/molecule_info.h5", cb_whitelist)  # cell calling
obj %<>% plot2()                                            # UMAP + metrics
obj %<>% plot3(df, puckdf, f)                               # raw spatial data
obj %<>% plot4(df, puckdf)                                  # beadplot
