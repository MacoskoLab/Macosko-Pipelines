# options(warn=1)

gdraw <- function(text, s=14) {ggdraw()+draw_label(text, size=s)}
plot.tab <- function(df) {return(plot_grid(tableGrob(df)))}
add.commas <- function(num){prettyNum(num, big.mark=",")}
make.pdf <- function(plots, name, w, h) {
  if ("gg" %in% class(plots) || class(plots)=="Heatmap") {plots = list(plots)}
  pdf(file=name, width=w ,height=h)
  lapply(plots, function(x){print(x)})
  dev.off()
}

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

### Load the seurat ############################################################
load_seurat <- function(matrix_path, molecule_info_path="") {
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
  obj %<>% process
  
  # Add %intronic
  if (file.exists(molecule_info_path)) {
    fetch <- function(x){return(h5read(molecule_info_path, x))}
    if ("barcode_idx" %in% h5ls(molecule_info_path)$name) {
      barcodes = fetch("barcodes")
      info = data.frame(barcode=fetch("barcode_idx")+1, umi_type=fetch("umi_type"))
      info %<>% group_by(barcode) %>% summarize(numi=n(), pct.intronic=sum(umi_type==0)/numi)
      obj$pct.intronic = info$pct.intronic[match(obj$cb, barcodes[info$barcode])] * 100
    } else {
      print("WARNING: no intronic information found in the molecule_info_path")
    }
  }
  
  return(obj)
}

### Load the puck ##############################################################
get_scaling_factor <- function(bn) {
  if (bn < 150000) {
    k = 0.73
  } else if (bn < 600000) {
    k = 0.73 * 2
  } else {
    k = 0.645
  }
  return(k)
}

# scale the coordinates (to um)
# scaling_factors = map_dbl(num_beads, get_scaling_factor)
# puckdfs %<>% map2(scaling_factors, ~transmute(.x, sb=sb, x_um=x*.y, y_um=y*.y))

### Plots ######################################################################

# Plot 0: cellranger metrics_summary.csv
plot_metrics_summary <- function(obj, metrics_summary_path, title) {
  plotdf = read.table(metrics_summary_path, header=F, sep=",", comment.char="")
  if (nrow(plotdf)==2) { # count
    plotdf %<>% t
    Misc(obj, "RNA_metrics_summary") <- setNames(plotdf[,2], plotdf[,1])
  } else if (ncol(plotdf)==6) { # multi
    colnames(plotdf) = as.character(plotdf[1,])
    plotdf = plotdf[-1,c(5,6)]
  }
  rownames(plotdf) = NULL
  plot = plot_grid(ggdraw()+draw_label(""),
                   ggdraw()+draw_label(title),
                   plot.tab(plotdf),
                   ggdraw()+draw_label(""),
                   ncol=1,rel_heights=c(0.1,0.1,0.7,0.2))
  make.pdf(plot,"plots/0cellranger.pdf",7,8)
  return(obj)
}
plot0 <- function(obj, metrics_summary_path, title = "Cell Ranger Metrics Summary") {
  if (file.exists(metrics_summary_path)) {
    obj %<>% plot_metrics_summary(metrics_summary_path, title)
  } else {
    make.pdf(gdraw("No metrics_summary.csv found"),"plots/0cellranger.pdf",7,8)
  }
  return(obj)
}


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







# Run alignment
# system("~/pipseeker-v3.0.5-linux/pipseeker full --fastq /home/nsachdev/pipseeker/fluent_RNA_fastqs/T20/T20 --star-index-path /home/nsachdev/pipseeker/fluent_RNA_fastqs/pipseeker-gex-reference-GRCh38-2022.04 --output-path /home/nsachdev/pipseeker/fluent_RNA_fastqs/T20 --chemistry v4 --threads 1 --verbosity 2 --retain-barcoded-fastqs --skip-version-check")
# system("~/pipseeker-v3.0.5-linux/pipseeker full --fastq /home/nsachdev/pipseeker/fluent_RNA_fastqs/T2/T2 --star-index-path /home/nsachdev/pipseeker/fluent_RNA_fastqs/pipseeker-gex-reference-GRCh38-2022.04 --output-path /home/nsachdev/pipseeker/fluent_RNA_fastqs/T2 --chemistry v4 --threads 1 --verbosity 2 --retain-barcoded-fastqs --skip-version-check")

# cargo run /home/nsachdev/pipseeker/fluent_RNA_fastqs/T20/T20_S1_R1_001.fastq.gz /home/nsachdev/pipseeker/fluent_RNA_fastqs/T20/barcoded_fastqs/barcoded_R1.fastq.gz

library(ShortRead) # sread, quality, id

# df %<>% group_by(cb, umi, sb) %>% summarize(reads = n())
count_reads <- function(df) {
  stopifnot(names(df) == c("cb", "umi", "sb"))
  gdf = df %>% select(cb,umi,sb) %>% arrange(cb,sb,umi)
  bnds = (gdf$cb!=lead(gdf$cb) | gdf$sb!=lead(gdf$sb) | gdf$umi!=lead(gdf$umi)) %>% tidyr::replace_na(T) %>% which
  gdf %<>% distinct()
  gdf$reads = (bnds - lag(bnds)) %>% tidyr::replace_na(bnds[[1]])
  gdf %<>% arrange(desc(reads))
  return(gdf)
}

pipseq_process <- function(R1path, R2path, mappingpath, RNApath, puckpath, metricssummarypath) {
  obj = load_seurat(RNApath)
  cb_whitelist = unname(colnames(obj))
  stopifnot(!duplicated(cb_whitelist))
  stopifnot(len(unique(nchar(cb_whitelist))) == 1)
  stopifnot(map_lgl(strsplit(cb_whitelist, ""), ~all(. %in% c("A","C","G","T"))))
  
  mapping = read.csv(mappingpath,header=F) %>% loadpipseqcbmapping
  
  puckdf = read.table(puckpath,header=T,sep=",") %>% setNames(c("sb","x","y"))
  
  print("loading fastqs")
  R1 = R1path %>% readFastq %>% sread %>% as.data.frame
  R2 = R2path %>% readFastq %>% sread %>% as.data.frame
  stopifnot(ncol(R1)==1, ncol(R2)==1, nrow(R1)==nrow(R2))
  R = data.frame(R1=R1[[1]],R2=R2[[1]])
  num_SB_reads = nrow(R)
  rm(R1,R2) ; gc()
  
  if (nchar(R$R1[[1]]) < 28 && nchar(R$R1[[1]]) < 28) {
    print("ERROR: READS TOO SHORT")
    stopifnot(F)
  }
  if (sum(str_sub(R$R1,16,25)=="CTGTTTCCTG") > sum(str_sub(R$R2,16,25)=="CTGTTTCCTG")) {
    print("switching R1 and R2")
    R = data.frame(R1 = R$R2, R2 = R$R1)
  } else {print("not switching")}
  
  print("processing dataframe")
  pct_no_up <- round(sum(str_sub(R$R2,16,25)!="CTGTTTCCTG")/nrow(R)*100,2)
  print(g("{pct_no_up}% reads have no CTGTTTCCTG from 16-25"))
  R %<>% filter(str_sub(R2,16,25)=="CTGTTTCCTG")
  R %<>% transmute(R1=R1, sb=str_sub(R2,1,14))
  gc()
  
  ms <- map(0:3, function(i){
    str_sub(R$R1, 9+i, 11+i)=="ATG" & str_sub(R$R1, 18+i, 20+i)=="GAG" & str_sub(R$R1, 27+i, 31+i)=="TCGAG"
  })
  m = Reduce(`|`, ms)
  pct_no_spacer <- round(sum(!m)/nrow(R)*100,2)
  print(g("{pct_no_spacer}% reads did not exact match the cell barcode spacer regions"))
  R = R[m,]
  ms %<>% map(~.[m])
  R$R1[ms[[1]]] %<>% str_sub(1, 51)
  R$R1[ms[[2]]] %<>% str_sub(2, 52)
  R$R1[ms[[3]]] %<>% str_sub(3, 53)
  R$R1[ms[[4]]] %<>% str_sub(4, 54)
  stopifnot(nchar(R$R1) == 51)
  R %<>% tidyr::separate(R1,into=c("cb","umi"),sep=c(39))
  stopifnot(nchar(R$cb) == 39)
  stopifnot(nchar(R$umi) == 12)
  stopifnot(nchar(R$sb) == 14)
  rm(ms, m) ; gc()
  
  pct_no_cb <- round(sum(!R$cb %in% names(mapping))/nrow(R)*100,2)
  print(g("{pct_no_cb}% reads did not exact match a cell barcode"))
  R %<>% filter(cb %in% names(mapping))
  R$cb = mapping[R$cb]
  
  pct_no_sb <- round(sum(!R$sb %in% puckdf$sb)/nrow(R)*100,2)
  print(g("{pct_no_sb}% reads did not exact match a spatial barcode"))
  R %<>% filter(sb %in% puckdf$sb)
  
  R %<>% count_reads
  R %<>% setNames(c("cb_index", "umi_2bit", "sb_index", "reads"))
  R %<>% count_umis
  R %<>% setNames(c("cb","sb","umi"))
  gc()
  
  m = match(R$sb, puckdf$sb)
  R$x_um = puckdf$x[m]*2*0.73*100
  R$y_um = puckdf$y[m]*2*0.73*100
  R %<>% setNames(c("cb_index","sb_index","umi","x_um","y_um"))
  
  R_filtered = R %>% filter(cb_index %in% cb_whitelist)
  R_filtered$sb_index = match(R_filtered$sb_index, puckdf$sb)
  R_filtered$cb_index = match(R_filtered$cb_index, cb_whitelist)
  res = normal_positioning(R_filtered)
  coords = res[[1]] ; data.list = res[[2]] ; rm(res)
  
  coords %<>% mutate(cb = cb_whitelist[cb_index])
  rownames(coords) = coords$cb
  obj = AddMetaData(obj,coords)
  Misc(obj,"pct.placed") = round(sum(!is.na(obj$x_um))/ncol(obj)*100,2)
  emb = obj@meta.data[,c("x_um","y_um")] ; colnames(emb) = c("s_1","s_2")
  obj[["spatial"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "s_", assay = "RNA")
  
  qsave(obj, "plots/seurat.qs")
  coords %>% select(cb,everything()) %>% write.csv("plots/coords.csv",quote=F,row.names=F)
  qsave(data.list, "plots/data.list.qs")
  qsave(mutate(R,called=cb_index%in%cb_whitelist), "plots/df.qs")
  obj %<>% plot0(metricssummarypath, title="PIPseeker Metrics Summary")
  obj %<>% plot2()
  
  sb.data = R %>% mutate(sb=sb_index) %>% group_by(sb) %>% summarize(umi=sum(umi)) %>%
    ungroup %>% merge(y=puckdf, all.x=T, by="sb") %>% mutate(x_um=x*2*0.73*100,y_um=y*2*0.73*100) %>% filter(!is.na(x_um),!is.na(y_um)) %>% arrange(umi)
  p1 = beadplot(sb.data)
  sb.data = R %>% mutate(sb=sb_index) %>% filter(cb_index %in% cb_whitelist) %>% group_by(sb) %>% summarize(umi=sum(umi)) %>%
    ungroup %>% merge(y=puckdf, all.x=T, by="sb") %>% mutate(x_um=x*2*0.73*100,y_um=y*2*0.73*100) %>% filter(!is.na(x_um),!is.na(y_um)) %>% arrange(umi)
  p2 = beadplot(sb.data)
  make.pdf(plot_grid(p1,p2,ncol=1), "plots/4beadplot.pdf", 7, 8)
  
  plot <- plot_dbscan(obj, coords)
  make.pdf(plot,"plots/5DBSCAN.pdf",7,8)
  
  p1 = DimPlot(obj[,!is.na(obj$x_um)],reduction="spatial")+coord_fixed(ratio=1)+
    ggtitle(g("%placed: {round(sum(!is.na(coords$x_um))/nrow(coords)*100,2)} ({sum(!is.na(obj$x_um))}/{ncol(obj)})")) + 
    NoLegend() + xlab("x-position (\u00B5m)") + ylab("y-position (\u00B5m)")
  p2 = DimPlot(obj[,!is.na(obj$x_um)], reduction="spatial",split.by="seurat_clusters",ncol=7) + theme_void() + coord_fixed(ratio=1) + NoLegend()
  plot = plot_grid(p1, p2, ncol=1, rel_heights=c(1,1))
  make.pdf(plot,"plots/6spatial.pdf",7,8)
  
  p8 = c(num_SB_reads = num_SB_reads,
       pct_no_up = pct_no_up %>% paste0("%"),
       pct_no_spacer = pct_no_spacer %>% paste0("%"),
       pct_no_cb = pct_no_cb %>% paste0("%"),
       pct_no_sb = pct_no_sb %>% paste0("%")) %>% {data.frame(metric=names(.),value=unname(.))} %>% plot.tab
  make.pdf(p8, "plots/7metrics.pdf", 7, 8)
  
  plots <- sample_bead_plots(data.list, range(R_filtered$x_um), range(R_filtered$y_um))
  make.pdf(plots,"plots/SB.pdf",7,7)
  
  pdfs = c("0cellranger.pdf", "2umap.pdf", "4beadplot.pdf", "5DBSCAN.pdf", "6spatial.pdf","7metrics.pdf", "SB.pdf") %>% paste0("plots/",.)
  qpdf::pdf_combine(input = pdfs, output = "plots/summary.pdf")
  
  qsave(obj, "plots/seurat.qs")

  print("Done!")
}

loadpipseqcbmapping <- function(df) {
  stopifnot(ncol(df)==3, typeof(df[[1]])=="character", typeof(df[[2]])=="character", typeof(df[[3]])%in%c("integer","numeric"))
  df %<>% setNames(c("V1","V2","V3"))
  sumbefore = sum(df$V3)
  df %<>% arrange(V1, desc(V3))
  before_same = tidyr::replace_na(df$V1==lag(df$V1), FALSE) 
  after_same = tidyr::replace_na(df$V1==lead(df$V1) & df$V3==lead(df$V3), FALSE)
  df = df[!before_same & !after_same,]
  df %<>% arrange(V2, desc(V3))
  before_same = tidyr::replace_na(df$V2==lag(df$V2), FALSE) 
  after_same = tidyr::replace_na(df$V2==lead(df$V2) & df$V3==lead(df$V3), FALSE)
  df = df[!before_same & !after_same,]
  sumafter = sum(df$V3)
  print(g("Ignoring {sumbefore-sumafter} discordant reads ({round((sumbefore-sumafter)/(sumbefore),10) %>% format(nsmall=10)}%)"))
  df %<>% select(V1,V2)
  stopifnot(!duplicated(df[[1]]), !duplicated(df[[2]]), nrow(df)==nrow(distinct(df)))
  mapping = setNames(df$V2,df$V1)
  return(mapping) # converts 28bp -> 16bp cell barcode
}

# setwd("~/output")
# system("mkdir plots")
# R1path = "/home/nsachdev/pipseeker/fluent_spatial_fastqs/T20S/T20S_S1_R1_001.fastq.gz"
# R2path = "/home/nsachdev/pipseeker/fluent_spatial_fastqs/T20S/T20S_S1_R2_001.fastq.gz"
# mappingpath = "/home/nsachdev/pipseeker/pipseqcbmapping/T20mapping.csv"
# RNApath = "~/pipseeker/fluent_RNA_fastqs/T20/filtered_matrix/sensitivity_4"
# puckpath = "/home/nsachdev/pipseeker/V15A_recon_loc_30000_80_0.2_fbmax1000.csv"
# metricssummarypath = "/home/nsachdev/pipseeker/fluent_RNA_fastqs/T20/metrics/sensitivity_4/metrics_summary.csv"
# pipseq_process(R1path, R2path, mappingpath, RNApath, puckpath, metricssummarypath)
# system("mv plots T20")

# setwd("~/output")
# system("mkdir plots")
# R1path = "/home/nsachdev/pipseeker/fluent_spatial_fastqs/T2S/T2S_S2_R1_001.fastq.gz"
# R2path = "/home/nsachdev/pipseeker/fluent_spatial_fastqs/T2S/T2S_S2_R2_001.fastq.gz"
# mappingpath = "/home/nsachdev/pipseeker/pipseqcbmapping/T2mapping.csv"
# RNApath = "~/pipseeker/fluent_RNA_fastqs/T2/filtered_matrix/sensitivity_4"
# puckpath = "/home/nsachdev/pipseeker/V15A_recon_loc_30000_80_0.2_fbmax1000.csv"
# metricssummarypath = "/home/nsachdev/pipseeker/fluent_RNA_fastqs/T2/metrics/sensitivity_4/metrics_summary.csv"
# pipseq_process(R1path, R2path, mappingpath, RNApath, puckpath, metricssummarypath)
# system("mv plots T2")
