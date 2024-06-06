# Run positioning for each cell in input dataframe (matrix.csv -> coords.csv)
library(glue) ; g=glue ; len=length
library(magrittr)
library(dbscan)
library(dplyr)
library(purrr)

matrix_path = "output/matrix.csv"
out_path = "output"
setwd("~/spatial")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  matrix_path <- args[[1]]
  out_path <- args[[2]]
} else {
  stop("Usage: Rscript positioning.R matrix_path output_path", call. = FALSE)
}
stopifnot(file.exists(matrix_path))
if (!dir.exists(out_path)) { dir.create(out_path, recursive = T) }

# load data.list
df = read.table(matrix_path, header=T, sep=",")
stopifnot(names(df) == c("cb_index", "x_um", "y_um", "umi"))
xlims = range(df$x_um)
ylims = range(df$y_um)
data.list = split(df, df$cb_index)
rm(df) ; invisible(gc())
print(g("Running positioning on {len(data.list)} cells"))

# prepare for positioning
library(furrr)
library(future)
library(parallel)
ncores = parallel::detectCores()
plan(multisession, workers=ncores)
chunk_vector <- function(v, chunk_size) {return(split(v, ceiling(seq_along(v) / chunk_size)))}

### Run DBSCAN #################################################################

# Do a grid search to find the ideal DBSCAN parameters
opt_dbscan <- function(data.list) {
  eps.vec = c(50, 100) ; minPts.vec = c(3:42)
  res = data.frame()
  for (k in 0:26) {
    params = expand.grid(eps.vec, minPts.vec) %>% setNames(c("eps","minPts"))
    row_lists = chunk_vector(1:nrow(params), ceiling(nrow(params)/ncores))
    
    params$pct = furrr::future_map(row_lists, function(v) {
      map_dbl(v, function(i) {
        m = map_int(data.list, ~dbscan::dbscan(.[c("x_um","y_um")], eps=params$eps[[i]], minPts=params$minPts[[i]], weights=.$umi)$cluster %>% max)
        return(sum(m==1)/length(m))
      })
    }, .options=furrr_options(seed=T)) %>% flatten_dbl
    
    res = rbind(res, params)
    
    p = res %>% group_by(minPts) %>% summarize(pct=max(pct)) %>% ungroup %>% arrange(minPts) %>% pull(pct) %>% {which.max(.)/len(.)}
    if (p < 0.9) {
      break
    }
    minPts.vec %<>% add(40)
  }
  
  res$is.max = res$pct==max(res$pct)
  eps = res$eps[res$is.max][[1]] ; minPts = res$minPts[res$is.max][[1]]
  pct.placed = round(max(res$pct)*100, 2)
  print(g("Optimal eps: {eps} \t Optimal minPts: {minPts} \t %placed: {pct.placed}"))
  
  return(c(eps, minPts, pct.placed))
}

# assign centroid and record metadata
create_dbscan_coords <- function(data.list) {
  coords = lapply(data.list, function(df) {
    p = c(x_um = NA,
          y_um = NA,
          DBSCAN_clusters = max(df$cluster),
          SB_umi = sum(df$umi),
          SNR = NA,
          minPts = unique(df$minPts) %>% {ifelse(is.null(.), NA, .)} %T>% {stopifnot(len(.)==1)},
          eps = unique(df$eps) %>% {ifelse(is.null(.), NA, .)} %T>% {stopifnot(len(.)==1)},
          pct.placed = unique(df$pct.placed) %>% {ifelse(is.null(.), NA, .)} %T>% {stopifnot(len(.)==1)}
         )
    if (max(df$cluster) == 1) {
      sdf = dplyr::filter(df, cluster==1)
      p[["x_um"]] = weighted.mean(sdf$x_um, w=sdf$umi)
      p[["y_um"]] = weighted.mean(sdf$y_um, w=sdf$umi)
      p[["SNR"]] = sum(sdf$umi)/sum(df$umi)
    }
    return(p)
  }) %>% bind_rows %>% as.data.frame %>% mutate(cb_index=as.numeric(names(data.list))) %>% select(cb_index, everything())
  return(coords)
}

params = opt_dbscan(data.list)
data.list %<>% lapply(function(df){
  mutate(df,
         cluster = dbscan::dbscan(df[c("x_um","y_um")], eps=params[[1]], minPts=params[[2]], weights=df$umi, borderPoints=F)$cluster,
         eps = params[[1]],
         minPts = params[[2]],
         pct.placed = params[[3]])
})
dbscan_coords <- create_dbscan_coords(data.list)
invisible(gc())

### Plot DBSCAN ################################################################

# make a ggplot in there

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









### Positioning ################################################################


















head(data.list)

plan(multisession, workers=20L) # plan(sequential, workers=20L)
NULLto0 <- function(val) {return(ifelse(is.null(val)||is.na(val)||len(val)==0, 0, val))}
chunk_vector <- function(v, chunk_size) {return(split(v, ceiling(seq_along(v) / chunk_size)))}
identify <- function(obj,clusts = NA) {
  if (!is.na(clusts)) {obj = obj[,obj$seurat_clusters %in% clusts]}
  obj = obj[,!is.na(obj$x_um)]
  which(anno.polygon(obj$x_um,obj$y_um,obj$seurat_clusters))
}

puckdf = data[[15]][[2]]
df = data[[15]][[1]] %>% count_umis %>% merge(y=puckdf,all.x=T,by="sb_index")
cb.data = df %>% group_by(cb_index) %>% dplyr::summarize(umi=sum(umi),n=n()) %>% ungroup

obj=qread("seurat15.qs")
obj=obj[,!obj$seurat_clusters %in% c(0,6,13,26,18)] %>% process

DimPlot(obj,label=T)

################################################################################

# remove df dependence
plot.sb <- function(subdf) {
  stopifnot(c("x_um","y_um","umi") %in%names(df)) ; stopifnot(cb_i %in% df$cb_index)
  subdf %<>% arrange(umi)
  ggplot() + coord_fixed(ratio=1, xlim=range(df$x_um), ylim=range(df$y_um)) + theme_void() +
    geom_point(data=subdf, mapping=aes(x=x_um,y=y_um,col=umi), size=2, shape=16) +
    theme(legend.key.width=unit(0.5,"lines"), legend.position="right", legend.key.height=unit(1,"lines"), legend.title=element_blank(), legend.spacing.y=unit(0.2,"lines"), legend.margin=margin(0,0,0,0,"lines"), legend.box.margin=margin(0,0,0,0,"pt"), legend.box.background=element_blank(), legend.background=element_blank(), legend.direction="vertical", legend.justification="left",legend.box.just="left",legend.box.spacing=unit(0,"cm"))
}

# remove df dependence
# note: density of lone point is num umi, then add umis from surroundings weighted by umi
generate_kde <- function(subdf, r=200, bw=20, plot=F) {
  cb_i = unique(subdf$cb_index) ; stopifnot(len(cb_i)==1)
  subdf %<>% select(x_um,y_um,umi)
  res = c(cb_index=cb_i,
          x=NA, y=NA,
          d1=NA, d2=NA,
          totumi=sum(subdf$umi),totbeads=nrow(subdf),
          umi=NA,beads=NA)
  if(nrow(subdf)==1) {res[c("x","y","d1","d2","umi","beads")]=c(subdf$x_um,subdf$y_um,subdf$umi,0,subdf$umi,1) ; return(res)}
  subdf %<>% filter(umi>0.1*max(subdf$umi))
  
  xmu = pdist(subdf[,c("x_um","y_um")])
  subdf$val = exp(-xmu^2/(2*bw^2)) %>% sweep(MARGIN=1, STATS=subdf$umi, FUN="*") %>% colSums
  rowmax = subdf[which.max(subdf$val),]
  subdf$near = (subdf$x_um-rowmax$x_um)^2 + (subdf$y_um-rowmax$y_um)^2 < r^2
  rowmax2 = subdf %>% filter(!near) %>% {.[which.max(.$val),]}
  #if (nrow(rowmax2)==0) {rowmax2 = data.frame(x_um=NA,y_um=NA,umi=NA,val=NA)}
  
  sdf <- subdf %>% filter(near)
  x = matrixStats::weightedMedian(sdf$x_um,w=sdf$umi)
  y = matrixStats::weightedMedian(sdf$y_um,w=sdf$umi)
  beads = sum(sdf$umi)/max(sdf$umi)
  res[c("x","y","d1","d2","umi","beads")]=c(x,y,rowmax$val,NULLto0(rowmax2$val),sum(sdf$umi),beads)
  
  if (plot) {
    subdf %<>% arrange(val)
    p <- ggplot(subdf, aes(x=x_um,y=y_um,col=val))+geom_point()+
      annotate("path",x=rowmax$x_um+r*cos(seq(0,2*pi,length.out=100)),y=rowmax$y_um+r*sin(seq(0,2*pi,length.out=100))) +
      geom_point(aes(x=!!x,y=!!y), colour="red", shape=0, size=4) +
      coord_fixed(ratio=1) + theme_void() + ggtitle(g("{cb_i} {round(res[['d2']],2)}/{round(res[['d1']],2)} ({round(res[['d2']]/res[['d1']],2)}) umi={round(res[['umi']],2)} beads={round(res[['beads']],2)} bw={bw}"))
    if (nrow(rowmax2)>0){p <- p + geom_point(aes(x=rowmax2$x_um,y=rowmax2$y_um), colour="blue", shape=1, size=4)}
    return(p)
  }
  
  return(res)
}

debug_position <- function(df, cb_i) {
  print(cb_i)
  stopifnot(c("x_um","y_um","umi") %in% names(df))
  stopifnot(cb_i %in% df$cb_index)
  subdf <- df %>% filter(cb_index==cb_i)
  
  range = 1:30*2
  search = data.frame(a=range, b=map_dbl(range, function(bw) {generate_kde(subdf, bw=bw, plot=F) %>% {.[["d2"]]/.[["d1"]]} }))
  if (max(search$b)==0) {newbw=20} else {newbw = range[[which.min(search$b)[[1]]]]}
  
  p1 = plot.sb(subdf)
  p2 = generate_kde(subdf,plot=T)
  p4 = generate_kde(subdf,bw=newbw,plot=T)
  p3 = search %>% ggplot(aes(x=a,y=b))+geom_point()+ylim(0,1)
  plot_grid(p1,p2,p3,p4,ncol=2)
}

debug_position(df=df,cb_i=sample(1:ncol(obj),1))
debug_position(df=df,cb_i=2070)
debug_position(df=df,cb_i=7089)


data.list = split(df, df$cb_index) %>% map(~arrange(.,umi))
cb_is = map_int(data.list,~unique(.$cb_index))


v=chunk_vector(data.list,ceiling(len(data.list)/20))
res = furrr::future_map(v,function(b){map(b,~generate_kde(.))}, .options=furrr_options(seed=T)) %>% unlist(recursive=F)
coords <- res %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("x_um","y_um","val1","val2","snr","n"))
coords %<>% mutate(cb_index=cb_is, p=val2/val1)
stopifnot(res2$cb_index == sort(res2$cb_index))


coords$placed = !is.na(obj$x_um)
ggplot(coords,aes(x=val1 %>% log10,y=val2 %>% log10,col=placed))+geom_point(size=0.2)+theme_bw()

obj$x = coords$x_um
obj$y = coords$y_um
emb = obj@meta.data[,c("x","y")] ; colnames(emb) = c("s_1","s_2")
obj[["spatial2"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_")

DimPlot(obj[,coords$p<0.5],reduction="spatial",split.by="seurat_clusters",ncol=6)

sum(tidyr::replace_na(coords$p,0)<0.5)/nrow(coords)
sum(!is.na(obj$x_um))/ncol(obj)




which((obj@reductions$umap@cell.embeddings[,1]< -5)&(obj$seurat_clusters==14)&(obj$outlier==T))

# try knn?
coord_fixed(ratio=1)

# compute the KDE

kde$m = map2_lgl(kde$x,kde$y,function(x,y){(irow$x-x)**2+(irow$y-y)**2<r**2})
orow = filter(kde,!m) %>% {.[which.max(.$val),]}
pd.emp = orow$val/irow$val

df$m = map2_lgl(df$x_um,df$y_um,function(x,y){(x-irow$x)**2+(y-irow$y)**2<r**2})
pr.emp = sum(filter(df,m)$umi)/sum(df$umi)

res = c(irow[c("x","y")], orow[c("x","y")], pd.emp, pr.emp) %>% setNames(c("x1","y1","x2","y2","pd.emp","pr.emp")) %>% unlist
return(res)



head(df)
head(data.list)
data.list %<>% map(~mutate(.,cluster=0))

plot.sb(data.list[[1]])
plot.kde <- function(subdf) {
  subdf %<>% arrange(umi)
  if(nrow(subdf)==1) {return()}
  p = Nebulosa:::wkde2d(x=subdf$x_um, y=subdf$y_um, w=subdf$umi, adjust=c(1,xrange/yrange)/exp(1), h=1, n=200, lims=c(range(df$x_um), range(df$y_um))) %>% {transmute(reshape2::melt(as.matrix(.[[3]])), x_um=.[[1]][Var1], y_um=.[[2]][Var2], value=value)}
  p[which.max(p$value),]
  
  p1 = ggplot(p, aes(x=x_um,y=y_um,fill=value))+geom_tile()+coord_fixed(ratio=1)
  p2 = plot.sb(subdf)
  plot_grid(p1,p2,ncol=1)
}
plot.kde(data.list[[6]])
