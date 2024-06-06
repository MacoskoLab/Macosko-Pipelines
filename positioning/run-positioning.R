library(glue) ; g=glue ; len=length
library(matrixStats)
library(gridExtra)
library(magrittr)
library(jsonlite)
library(ggplot2)
library(cowplot)
library(stringr)
library(Seurat)
library(rlist)
library(dplyr)
library(purrr)
library(qpdf)
library(qs)

setwd("~/spatial")
RNApath = "RNAcounts"
SBpath = "SBcounts"
out_path <- "output"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  RNApath <- args[[1]]
  SBpath <- args[[2]]
  out_path <- "output"
} else if (length(args) == 3) {
  RNApath <- args[[1]]
  SBpath <- args[[2]]
  out_path <- args[[3]]
} else {
  stop("Usage: Rscript positioning.R RNApath SBpath [output_path]", call. = FALSE)
}
stopifnot(dir.exists(RNApath), len(list.files(RNApath)) > 0)
if (dir.exists(SBpath)) { SBpath = file.path(SBpath, "SBcounts.h5") }
stopifnot(str_sub(SBpath,-11,-1) == "SBcounts.h5", file.exists(SBpath))
if (!dir.exists(out_path)) { dir.create(out_path, recursive = T) }
print(g("RNA dir: {normalizePath(RNApath)}"))
print(g("RNA files: {list.files(RNApath) %>% paste0(collapse=', ')}"))
print(g("SB file: {normalizePath(SBpath)}"))

# Determine the RNA method
if (file.exists(file.path(RNApath,"filtered_feature_bc_matrix.h5"))) {
  method = "10X"
} else {
  stop("ERROR: unknown RNA technique", call. = FALSE)
}

### Helper methods #############################################################

gdraw <- function(text, s=14) {ggdraw()+draw_label(text, size=s)}
plot.tab <- function(df) {return(plot_grid(tableGrob(df)))}
add.commas <- function(num){prettyNum(num, big.mark=",")}
make.pdf <- function(plots, name, w, h) {
  if ("gg" %in% class(plots) || class(plots)=="Heatmap") {plots = list(plots)}
  pdf(file=name, width=w ,height=h)
  lapply(plots, function(x){print(x)})
  dev.off()
}

### Load the seurat ############################################################

load_seurat <- function(matrix_path) {
  # Load the RNA count matrix
  if (str_sub(matrix_path,-3,-1) == ".h5") {
    obj <- matrix_path %>% Read10X_h5 %>% CreateSeuratObject
  } else {
    obj <- matrix_path %>% Read10X %>% CreateSeuratObject
  }
  
  # Add metadata
  obj[["cb"]] <- map_chr(colnames(obj), ~sub("-[0-9]*$", "", .))
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

plot_metrics_summary <- function(summary_path, title, out_dir) {
  plotdf = read.table(summary_path, header=F, sep=",", comment.char="")
  plotdf %<>% t
  rownames(plotdf) = NULL
  
  plot = plot_grid(ggdraw()+draw_label(""),
                   ggdraw()+draw_label(title),
                   plot.tab(plotdf),
                   ggdraw()+draw_label(""),
                   ncol=1, rel_heights=c(0.1,0.1,0.7,0.2))
  make.pdf(plot, file.path(out_dir,"metrics_rna.pdf"), 7, 8)
  
  return(setNames(plotdf[,2], plotdf[,1]))
}

if (method == "10X") {
  matrix_path = file.path(RNApath, "filtered_feature_bc_matrix.h5")
  obj <- load_seurat(matrix_path)
  Misc(obj, "method") <- method
  
  molecule_info_path = file.path(RNApath, "molecule_info.h5")
  if (file.exists(molecule_info_path)) { obj %<>% load_intronic(molecule_info_path) }
  
  summary_path = file.path(RNApath, "metrics_summary.csv")
  if (file.exists(summary_path)) {
    Misc(obj, "RNA_metadata") = plot_metrics_summary(summary_path, "Cell Ranger Metrics Summary", out_path)
  } else {
    make.pdf(gdraw("No metrics_summary.csv found"), file.path(out_path,"metrics_rna.pdf"), 7, 8)
    rna_meta = list()
  }
} else {
  stop(g("ERROR: unknown RNA technique ({method})"), call. = FALSE)
}

print(g("Loaded {ncol(obj)} cells"))

### Generate the matrix ########################################################

print("Generating the matrix...")

cb_whitelist = unname(obj$cb)
writeLines(cb_whitelist, file.path(out_path, "cb_whitelist.txt"))
system(g("Rscript load_matrix.R {SBpath} {file.path(out_path, 'cb_whitelist.txt')} {out_path}"))
# df = read.table(file.path(out_path, "matrix.csv"), header=T, sep=",")
Misc(obj, "SB_metadata") <- fromJSON(file.path(out_path, "metadata.json"))

print("Matrix loading complete")

stopifnot(F)

### Positioning methods 2 ######################################################



# Plot 1: cell calling and %intronic
UvsI <- function(obj, molecule_info_path, cb_whitelist) {
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
      
      max_density_x = density(filter(df,umi>=ct,pct.intronic>0.35)$pct.intronic) %>% {.$x[which.max(.$y)]}
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
  df %<>% mutate(index=1:nrow(df), called=barcodes[barcode+1]%in%cb_whitelist)
  p3 = ggplot(df,aes(x=index,y=umi,col=called))+geom_line()+theme_bw()+scale_x_log10()+scale_y_log10()+
    ggtitle("Barcode rank plot")+xlab("Cell barcodes")+ylab("UMI counts") +
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"))
  
  plot = plot_grid(p3,p1,p0,p2,ncol=2)
  make.pdf(plot, "plots/1cellcalling.pdf", 7, 8)
  return(obj)
}
plot1 <- function(obj, molecule_info_path, cb_whitelist) {
  if (file.exists(molecule_info_path)) {
    obj %<>% UvsI(molecule_info_path, cb_whitelist)
  } else {
    make.pdf(gdraw("No molecule_info.h5 found"), "plots/1cellcalling.pdf", 7, 8)
  }
  return(obj)
}

# Plot 2: UMAP + metrics
plot2 <- function(obj) {
  mytheme <- function(){theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank())}
  plot = plot_grid(DimPlot(obj,label=T)+ggtitle(g("UMAP"))+NoLegend()  +mytheme()+coord_fixed(ratio=1),
                   VlnPlot(obj,"logumi",alpha=0)+NoLegend()            +mytheme(),
                   FeaturePlot(obj,"percent.mt")+ggtitle("%MT")        +mytheme()+coord_fixed(ratio=1)+theme(legend.position="top",legend.justification="center",legend.key.width=unit(2, "lines")),
                   if ("pct.intronic" %in% names(obj@meta.data)) {
                     FeaturePlot(obj,"pct.intronic")+ggtitle("%Intronic")+mytheme()+coord_fixed(ratio=1)+theme(legend.position="top",legend.justification="center",legend.key.width=unit(2,"lines"))
                   } else {
                     gdraw("No intronic information")
                   },
                   ncol=2)
  make.pdf(plot,"plots/2umap.pdf",7,8)
  return(obj)
}

# merge with seurat object
coords %<>% mutate(cb = cb_whitelist[cb_index])
rownames(coords) = paste0(coords$cb,"-1")
obj = AddMetaData(obj,coords)
Misc(obj,"pct.placed") = round(sum(!is.na(obj$x_um))/ncol(obj)*100,2)
emb = obj@meta.data[,c("x_um","y_um")] ; colnames(emb) = c("s_1","s_2")
obj[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_", assay = "RNA")

coords %>% select(cb,everything()) %>% write.csv("coords.csv",quote=F,row.names=F)
qsave(obj, "seurat.qs")

### Page 6: Spatial ############################################################

p1 = DimPlot(obj[,!is.na(obj$x_um)],reduction="spatial")+coord_fixed(ratio=1)+
  ggtitle(g("%placed: {round(sum(!is.na(coords$x_um))/nrow(coords)*100,2)} ({sum(!is.na(obj$x_um))}/{ncol(obj)})")) + 
  NoLegend() + xlab("x-position (\u00B5m)") + ylab("y-position (\u00B5m)")
p2 = DimPlot(obj[,!is.na(obj$x_um)], reduction="spatial",split.by="seurat_clusters",ncol=5) + theme_void() + coord_fixed(ratio=1) + NoLegend()
plot = plot_grid(p1, p2, ncol=1, rel_heights=c(1,1))
make.pdf(plot,"plots/6spatial.pdf",7,8)

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
