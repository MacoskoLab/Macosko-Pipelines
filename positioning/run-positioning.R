library(glue) ; g=glue ; len=length
library(gridExtra)
library(magrittr)
library(jsonlite)
library(ggplot2)
library(cowplot)
library(viridis)
library(stringr)
library(Seurat)
library(rlist)
library(dplyr)
library(purrr)
library(qpdf)
library(qs)

# Helper files required for loading spatial data and positioning cells
stopifnot(file.exists("load_matrix.R", "positioning.R"))

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

# Check arguments
stopifnot(dir.exists(RNApath), len(list.files(RNApath)) > 0)
if (dir.exists(SBpath)) { SBpath = file.path(SBpath, "SBcounts.h5") }
stopifnot(str_sub(SBpath, -11, -1) == "SBcounts.h5", file.exists(SBpath))
if (!dir.exists(out_path)) { dir.create(out_path, recursive = T) }
stopifnot("Could not create output path" = dir.exists(out_path))
print(g("RNA dir: {normalizePath(RNApath)}"))
print(g("RNA files: {list.files(RNApath) %>% paste0(collapse=', ')}"))
print(g("SB file: {normalizePath(SBpath)}"))
print(g("Output dir: {normalizePath(out_path)}"))

# Determine the RNA method
if (file.exists(file.path(RNApath, "filtered_feature_bc_matrix.h5"))) {
  method = "10X"
} else if (dir.exists(file.path(RNApath, "filtered_feature_bc_matrix"))) {
  method = "10X"
} else {
  stop("ERROR: unknown RNA technique", call. = FALSE)
}

# Helper methods
gdraw <- function(text, s=14) {ggdraw()+draw_label(text, size=s)}
plot.tab <- function(df) {return(plot_grid(tableGrob(df,rows=NULL)))}
add.commas <- function(num){prettyNum(num, big.mark=",")}
make.pdf <- function(plots, name, w, h) {
  if ("gg" %in% class(plots) || class(plots)=="Heatmap") {plots = list(plots)}
  pdf(file=name, width=w, height=h)
  lapply(plots, function(x){print(x)})
  dev.off()
}

### Load the RNA ###############################################################

load_seurat <- function(matrix_path) {
  # Load the RNA count matrix
  if (str_sub(matrix_path,-3,-1) == ".h5") {
    obj <- matrix_path %>% Read10X_h5 %>% CreateSeuratObject
  } else {
    obj <- matrix_path %>% Read10X %>% CreateSeuratObject
  }
  
  # Add metadata
  obj[["cb"]] <- map_chr(colnames(obj), ~sub("-[0-9]*$", "", .)) %T>% {stopifnot(!duplicated(.))}
  obj[["cb_index"]] <- 1:ncol(obj)
  obj[["logumi"]] <- log10(obj$nCount_RNA+1)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^(MT-|mt-)")
  
  # PCA, Cluster, and UMAP
  process <- function(obj, res=0.8, n.epochs=NULL) {
    obj <- obj %>%
      Seurat::NormalizeData() %>%
      Seurat::FindVariableFeatures() %>%
      Seurat::ScaleData() %>%
      Seurat::RunPCA(npcs=50, verbose=F) %>%
      Seurat::FindNeighbors(dims=1:30) %>%
      Seurat::FindClusters(resolution=res) %>%
      Seurat::RunUMAP(dims=1:30, verbose=F, n.epochs=n.epochs)
  }
  obj %<>% process
  
  return(obj)
}

load_intronic <- function(obj, molecule_info_path) {
  library(rhdf5)
  fetch <- function(x){return(h5read(molecule_info_path, x))}
  # Add %intronic
  if ("barcode_idx" %in% h5ls(molecule_info_path)$name) {
    barcodes = fetch("barcodes")
    info = data.frame(barcode=fetch("barcode_idx")+1, umi_type=fetch("umi_type"))
    info %<>% group_by(barcode) %>% summarize(numi=n(), pct.intronic=sum(umi_type==0)/numi)
    obj$pct.intronic = info$pct.intronic[match(obj$cb, barcodes[info$barcode])] * 100
  } else {
    print("WARNING: no intronic information found in the molecule_info_path")
  }
  return(obj)
}

plot_metrics_summary <- function(summary_path) {
  if (nchar(summary_path)==0 || !file.exists(summary_path)) {
    make.pdf(gdraw("No metrics_summary.csv found"), file.path(out_path,"RNAmetrics.pdf"), 7, 8)
    return(c())
  }
  plotdf = read.table(summary_path, header=F, sep=",", comment.char="")
  plotdf %<>% t
  rownames(plotdf) = NULL
  
  plot = plot_grid(ggdraw()+draw_label(""),
                   ggdraw()+draw_label(g("{method} Metrics Summary")),
                   plot.tab(plotdf),
                   ggdraw()+draw_label(""),
                   ncol=1, rel_heights=c(0.1,0.1,0.7,0.2))
  make.pdf(plot, file.path(out_path, "RNAmetrics.pdf"), 7, 8)
  
  return(setNames(plotdf[,2], plotdf[,1]))
}

# get the data paths
if (method == "10X") {
  matrix_path = file.path(RNApath, "filtered_feature_bc_matrix.h5")
  molecule_info_path = file.path(RNApath, "molecule_info.h5")
  summary_path = file.path(RNApath, "metrics_summary.csv")
} else {
  stop(g("ERROR: unknown RNA technique ({method})"), call. = FALSE)
}

# Load the RNA data
obj <- load_seurat(matrix_path)
Misc(obj, "method") <- method
if (file.exists(molecule_info_path)) { obj %<>% load_intronic(molecule_info_path) }
Misc(obj, "RNA_metadata") = plot_metrics_summary(summary_path)

print(g("Loaded {ncol(obj)} cells"))


# Plot RNA curves
UvsI <- function(obj, molecule_info_path) {
  if (!file.exists(molecule_info_path) || nchar(molecule_info_path) == 0) {
    plot <- gdraw("No molecule_info.h5 found")
    return(plot)
  }
  if (!"barcode_idx" %in% h5ls(molecule_info_path)$name) {
    plot <- gdraw("Unrecognized molecule_info.h5")
    return(plot)
  }
  
  fetch <- function(x){return(h5read(molecule_info_path,x))}
  barcodes = fetch("barcodes")
  molecule_info = data.frame(barcode=fetch("barcode_idx"),
                             umi_type=fetch("umi_type"),
                             reads=fetch("count"))
  
  # Panel 3: downsampling curve
  tab = table(molecule_info$reads)
  downsampling = map_int(seq(0,1,0.05), function(p){
    map2_int(tab, as.numeric(names(tab)), function(v, k){
      length(unique(floor(sample(0:(k*v-1), round(k*v*p), replace=F)/k)))
    }) %>% sum
  })
  plotdf = data.frame(
    x = seq(0,1,0.05)*sum(molecule_info$reads)/1000/1000,
    y = downsampling/1000/1000
  )
  p0 = ggplot(plotdf, aes(x=x,y=y))+geom_line()+theme_bw()+xlab("Millions of reads")+ylab("Millions of filtered UMIs")+ggtitle("RNA Downsampling curve")
  
  df = molecule_info %>% group_by(barcode) %>% summarize(umi=n(), pct.intronic=sum(umi_type==0)/umi) %>% 
    ungroup %>% arrange(desc(umi)) %>% mutate(logumi=log10(umi))
  
  # Panel 2 and 4: intronic density
  if (!is.null(df$pct.intronic) && !all(df$pct.intronic==0)) {
    ct = 500
    if (any(df$umi>=ct)) {
      p1 = df %>% filter(umi>=ct) %>% ggplot(aes(x = logumi, y = pct.intronic)) + 
        geom_bin2d(bins=100) +
        scale_fill_viridis(trans="log", option="A", name="density") + 
        theme_minimal() +
        labs(title = g("Intronic vs. UMI droplets (>{ct} umi)"), x = "logumi", y = "%intronic") + NoLegend()
      
      max_density_x = density(filter(df,umi>=ct,pct.intronic>1/3)$pct.intronic) %>% {.$x[which.max(.$y)]}
      p2 = df %>% filter(umi>=ct) %>% ggplot(aes(x = pct.intronic)) + geom_density() + 
        theme_minimal() + labs(title = g("Intronic density (>{ct} umi)"), x = "%intronic", y = "Density") + 
        geom_vline(xintercept = max_density_x, color = "red", linetype = "dashed") +
        annotate(geom = 'text', label = round(max_density_x, 2), x = max_density_x+0.01, y = Inf, hjust = 0, vjust = 1, col="red")
    } else {
      p1 = gdraw(g("No cells with {ct}+ UMI"))
      p2 = gdraw(g("No cells with {ct}+ UMI"))
    }
  } else {
    p1 = gdraw(g("No intronic information"))
    p2 = gdraw(g("No intronic information"))
  }
  
  # Panel 1: cell barcode knee plot
  df %<>% mutate(index=1:nrow(df), called=barcodes[barcode+1] %in% obj$cb)
  p3 = ggplot(df,aes(x=index,y=umi,col=called))+geom_line()+theme_bw()+scale_x_log10()+scale_y_log10()+
    ggtitle("Barcode rank plot")+xlab("Cell barcodes")+ylab("UMI counts") +
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"))
  
  plot = plot_grid(p3, p1, p0, p2, ncol=2)
  return(plot)
}
plot <- UvsI(obj, molecule_info_path)
make.pdf(plot, file.path(out_path,"RNA.pdf"), 7, 8)

# Plot UMAP + metrics
plot_umaps <- function(obj) {
  mytheme <- function(){theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="top", legend.justification="center", legend.key.width=unit(2, "lines"))}
  umap <- DimPlot(obj,label=T) + ggtitle(g("UMAP")) + mytheme() + NoLegend()  + coord_fixed(ratio=1)
  logumi <- VlnPlot(obj,"logumi",alpha=0) + mytheme() + NoLegend() 
  mt <- FeaturePlot(obj,"percent.mt") + ggtitle("%MT") + mytheme() + coord_fixed(ratio=1) + 
    annotate("text", x = Inf, y = Inf, label = g("Median: {round(median(obj$percent.mt), 2)}%\nMean: {round(mean(obj$percent.mt), 2)}%"), hjust=1, vjust=1, size=2.5, color="black")
  if ("pct.intronic" %in% names(obj@meta.data)) {
    intronic <- FeaturePlot(obj,"pct.intronic") + ggtitle("%Intronic") + mytheme() + coord_fixed(ratio=1) +
      annotate("text", x = Inf, y = Inf, label = g("Median: {round(median(obj$pct.intronic), 2)}%\nMean: {round(mean(obj$pct.intronic), 2)}%"), hjust=1, vjust=1, size=2.5, color="black")
  } else {
    intronic <- gdraw("No intronic information")
  }
  
  plot <- plot_grid(umap, logumi, mt, intronic, ncol=2)
  return(plot)
}
plot <- plot_umaps(obj)
make.pdf(plot, file.path(out_path,"UMAP.pdf"), 7, 8)

### Load the Spatial ###########################################################

# Load the spatial barcode counts matrix and fuzzy match to our called-cells whitelist
print("Generating the matrix...")
cb_whitelist = unname(obj$cb)
writeLines(cb_whitelist, file.path(out_path, "cb_whitelist.txt"))
system(g("Rscript load_matrix.R {SBpath} {file.path(out_path, 'cb_whitelist.txt')} {out_path}"))
stopifnot(file.exists(file.path(out_path, "matrix.csv.gz"), file.path(out_path, "spatial_metadata.json")))
Misc(obj, "spatial_metadata") <- fromJSON(file.path(out_path, "spatial_metadata.json"))

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

# Create DimPlot
plot_clusters <- function(obj, reduction) {
  npucks = (max(obj$x_um,na.rm=T)-min(obj$x_um,na.rm=T))/(max(obj$y_um,na.rm=T)-min(obj$y_um,na.rm=T))
  nclusters = len(unique(obj$seurat_clusters))
  ncols = round(sqrt(npucks*nclusters/2)/npucks*2) 
  
  m = obj@reductions[[reduction]]@cell.embeddings %>% {!is.na(.[,1]) & !is.na(.[,2])}
  title = g("%placed: {round(sum(m)/len(m)*100,2)} ({sum(m)}/{len(m)}) [{reduction}]")
  p1 = DimPlot(obj, reduction=reduction) + coord_fixed(ratio=1) +
    ggtitle(title) + NoLegend() + xlab("x-position (\u00B5m)") + ylab("y-position (\u00B5m)") + 
    theme(axis.title.x=element_text(size=12), axis.text.x=element_text(size=10)) +
    theme(axis.title.y=element_text(size=12), axis.text.y=element_text(size=10))
  p2 = DimPlot(obj, reduction=reduction, split.by="seurat_clusters", ncol=ncols) + theme_void() + coord_fixed(ratio=1) + NoLegend()
  plot = plot_grid(p1, p2, ncol=1, rel_heights=c(0.4,0.6))
  return(plot)
}
plot <- plot_clusters(obj, reduction="dbscan")
make.pdf(plot, file.path(out_path,"DimPlotDBSCAN.pdf"), 7, 8)
plot <- plot_clusters(obj, reduction="kde")
make.pdf(plot, file.path(out_path,"DimPlotKDE.pdf"), 7, 8)
plot <- plot_clusters(obj, reduction="spatial")
make.pdf(plot, file.path(out_path,"DimPlot.pdf"), 7, 8)

# RNA vs SB metrics
plot_RNAvsSB <- function(obj) {
  if (!is.null(Misc(obj,"coords")$umi)) {
    obj$sb_umi <- Misc(obj,"coords")$umi %>% tidyr::replace_na(0)
  } else if (!is.null(Misc(obj,"coords")$umi_dbscan)) {
    obj$sb_umi <- Misc(obj,"coords")$umi_dbscan %>% tidyr::replace_na(0)
  } else {
    obj$sb_umi = 0
  }
  
  obj$clusters <- Misc(obj,"coords")$clusters %>% tidyr::replace_na(0)
  obj$placed <- !is.na(obj$x_um) & !is.na(obj$y_um)
  
  p1 <- ggplot(obj@meta.data, aes(x=log10(nCount_RNA), y=log10(sb_umi), col=placed)) + 
    geom_point(size=0.2) + theme_bw() + xlab("log10 RNA UMI") + ylab("log10 SB UMI") + ggtitle("SB UMI vs. RNA UMI") + 
    labs(color = "placed") +
    theme(legend.position = c(0.95, 0.05),
          legend.justification = c("right", "bottom"),
          legend.background = element_blank(),
          legend.title=element_text(size=10),
          legend.text=element_text(size=8),
          legend.margin=margin(0,0,0,0,"pt"),
          legend.box.margin=margin(0,0,0,0,"pt"),
          legend.spacing.y = unit(0.1,"lines"),
          legend.key.size = unit(0.5, "lines"))
  
  d = obj@meta.data %>% rowwise %>% mutate(x=min(clusters,5)) %>% ungroup
  p2 <- ggplot(d, aes(x=as.factor(x), y=log10(nCount_RNA))) + geom_violin(scale="count") + 
    scale_x_discrete(breaks=min(d$x):max(d$x), labels=(min(d$x):max(d$x)) %>% {ifelse(.==5, "5+", .)}) +
    xlab("DBSCAN clusters") + ylab("log10 RNA UMI") + ggtitle("RNA UMI vs. DBSCAN cluster") + theme_classic()
  
  d = obj@meta.data %>% group_by(seurat_clusters) %>% summarize(pct.placed=paste0(round(sum(placed)/n()*100,2),"%")) %>% setNames(c("cluster","placed"))
  m = ceiling(nrow(d)/2) ; d1 = d[1:m,] ; d2 = d[(m+1):nrow(d),]
  p3 <- plot_grid(plot.tab(d1), plot.tab(d2), ncol=2)
  
  plot = plot_grid(gdraw("RNA vs. SB metrics"),
                   plot_grid(plot_grid(p1,p2,ncol=1), p3, ncol=2, rel_widths=c(0.5,0.5)),
                   ncol=1, rel_heights=c(0.05,0.95))
  return(plot)
}
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
