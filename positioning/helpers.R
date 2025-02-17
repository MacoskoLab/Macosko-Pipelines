### Helper methods #############################################################

gdraw <- function(text, s=14) {cowplot::ggdraw()+cowplot::draw_label(text, size=s)}
plot.tab <- function(df) {return(cowplot::plot_grid(gridExtra::tableGrob(df, rows=NULL)))}
dec2pct <- function(x, d=2) {return(paste0(round(x*100, digits=d), "%"))}
add.commas <- function(num) {prettyNum(num, big.mark=",")}
make.pdf <- function(plots, name, w, h) {
  if (any(c("gg", "ggplot", "Heatmap") %in% class(plots))) {plots = list(plots)}
  pdf(file=name, width=w, height=h)
  lapply(plots, print)
  dev.off()
}

remap_10X_CB <- function(vec) {
  stopifnot(class(vec) == "character")
  stopifnot(nchar(vec) == 16)
  basemap = setNames(c("AG","TC","CA","GT"), c("TC","AG","GT","CA"))
  stopifnot(unique(substr(vec,8,9)) %in% names(basemap))
  ret = paste0(substr(vec,1,7), basemap[substr(vec,8,9)], substr(vec,10,16))
  stopifnot(len(vec) == len(ret))
  stopifnot(nchar(vec) == nchar(ret))
  return(ret)
}

listHD1neighbors <- function(input_string) {
  nucleotides <- c('A','C','G','T','N')
  result <- c()
  for (i in 1:nchar(input_string)) {
    for (nuc in nucleotides[nucleotides != substr(input_string, i, i)]) {
      new_string <- paste0(substr(input_string, 1, i-1), nuc, substr(input_string, i+1, nchar(input_string)))
      result <- c(result, new_string)
    }
  }
  return(result)
}

count_umis <- function(df) {
  # return(df %>% group_by(cb_index, sb_index) %>% summarize(umi=n()) %>% arrange(desc(umi)))
  stopifnot(ncol(df) == 4, names(df) == c("cb_index", "umi_2bit", "sb_index", "reads"))
  gdf = df %>% filter(reads > 0) %>% select(cb_index, sb_index) %>% arrange(cb_index, sb_index) 
  bnds = (gdf$cb_index!=lead(gdf$cb_index) | gdf$sb_index!=lead(gdf$sb_index)) %>% tidyr::replace_na(T) %>% which
  gdf %<>% distinct()
  gdf$umi = (bnds - lag(bnds)) %>% tidyr::replace_na(bnds[[1]])
  gdf %<>% arrange(desc(umi))
  return(gdf)
}

### Loading methods ############################################################

ReadOptimus <- function(matrix_path) {
  fetch <- function(x){return(rhdf5::h5read(matrix_path, x))}
  stopifnot("Unrecognized .h5ad matrix" = c("data","indices","indptr") %in% rhdf5::h5ls(matrix_path)$name)
  
  data <- fetch("X/data") %>% as.numeric
  indices <- fetch("X/indices") %>% as.numeric
  indptr <- fetch("X/indptr") %>% as.numeric
  
  row_names <- fetch("var/_index") %>% as.character # var/Gene and var/gene_names were too short
  col_names <- fetch("obs/_index") %>% as.character # obs/CellID and obs/cell_names were the same
  star_IsCell <- fetch("obs/star_IsCell") %>% as.logical
  
  mat <- Matrix::sparseMatrix(
    j = indices,
    p = indptr,
    x = data,
    dimnames = list(row_names, col_names),
    index1 = F,
    check = T
  )
  
  stopifnot(ncol(mat) == len(star_IsCell))
  return(mat[,star_IsCell])
} # REMAP GENE NAMES

ReadIntronic <- function(intronic_path, cb_list) {
  if (!file.exists(intronic_path)) {
    return(rep(NA, length(cb_list)))
  }
  
  fields <- rhdf5::h5ls(intronic_path)$name
  fetch <- function(x){return(rhdf5::h5read(intronic_path, x))}
  
  if (all(c("barcodes", "barcode_idx", "umi_type") %in% fields)) { # 10X
    barcodes = fetch("barcodes")
    stopifnot(!any(duplicated(barcodes)), cb_list %in% barcodes)
    dt <- data.table(barcode_idx=fetch("barcode_idx")+1,
                     umi_type=fetch("umi_type"))
    pct.intronic <- dt[, .(total_umi=.N, intronic_umi=sum(umi_type==0)), barcode_idx][cb_list %>% match(barcodes) %>% match(barcode_idx), intronic_umi/total_umi*100]
  } else if (all(c("reads_mapped_intronic", "reads_mapped_exonic") %in% fields)) { # Optimus
    dt <- data.table(barcode = fetch("/obs/CellID"),
                     intronic = fetch("/obs/reads_mapped_intronic"),
                     exonic = fetch("/obs/reads_mapped_exonic"))
    stopifnot(!any(duplicated(dt[,barcode])), cb_list %in% dt[,barcode])
    pct.intronic <- dt[, .(pct.intronic=intronic/(intronic+exonic)*100), keyby=barcode][cb_list, pct.intronic]
  } else {
    print("WARNING: no intronic information found in the .h5")
    pct.intronic <- rep(NA, length(cb_list))
  }
  return(pct.intronic)
}

ReadLibraryMetrics <- function(metrics_path) {
  if (!file.exists(metrics_path)) {
    return(data.frame())
  }
  
  df = read.table(metrics_path, header=F, sep=",", comment.char="")
  if (ncol(df) > nrow(df)) { df %<>% t } # 10X
  df %<>% as.data.frame
  rownames(df) <- NULL
  colnames(df) <- NULL
  
  if (is.numeric(df[,2])) {
    df[,2] %<>% {case_when(
      . > 0 & . < 1 ~ dec2pct(.),
      . > 1 ~ add.commas(.),
      TRUE ~ as.character(.)
    )}
  }
  
  return(df)
}

ReadSpatialMetadata <- function() {
  
}

### Plotting methods ###########################################################

# Page 1: 10X/Optimus RNA metrics
plot_metrics_csv <- function(df) {
  if (nrow(df) == 0) {
    return(gdraw("No metrics_summary.csv found"))
  }
  
  hide_cols <- c("NHashID","keeper_cells","percent_keeper","percent_target","percent_usable","keeper_mean_reads_per_cell","keeper_median_genes","total_genes_unique_detected")
  plot <- plot_grid(ggdraw()+draw_label(""),
                    ggdraw()+draw_label(g("Metrics Summary")),
                    plot.tab(df[!df[,1] %in% hide_cols,]),
                    ggdraw()+draw_label(""),
                    ncol=1, rel_heights=c(0.1,0.1,0.7,0.2))
  return(plot)
}

# Page 2: Cell calling
plot_cellcalling <- function(intronic_path, cb_list, ct = 200) {
  # Check if the supplementary .h5 exists
  if (!file.exists(intronic_path) || nchar(intronic_path) == 0) {
    plot <- gdraw("No supplementary .h5 found")
    return(plot)
  }
  
  # Load the supplementary .h5
  fields <- rhdf5::h5ls(intronic_path)$name
  fetch <- function(x){return(rhdf5::h5read(intronic_path, x))}
  if (all(c("barcodes", "barcode_idx", "umi_type") %in% fields)) { # 10X
    dt <- data.table(barcode_idx=fetch("barcode_idx")+1,
                     umi_type=fetch("umi_type")
                    )[,.(umi=.N, pct.intronic=sum(umi_type==0)/.N*100), barcode_idx]
    dt[, called := barcode_idx %in% match(cb_list, fetch("barcodes"))]
    dt[, barcode_idx := NULL]
    dt2 <- data.table(x=fetch("count"))[,.(N=.N),x]
  } else if (all(c("reads_mapped_intronic", "reads_mapped_exonic") %in% fields)) { # Optimus
    dt <- data.table(barcode = fetch("/obs/CellID"),
                     umi = data.table(data=fetch("X/data"),indices=fetch("X/indices"))[,.(umi=sum(data)),indices][order(indices), umi],
                     pct.intronic = fetch("/obs/reads_mapped_intronic") %>% {./(.+fetch("/obs/reads_mapped_exonic"))*100}
                    )
    dt[, called := barcode %in% cb_list]
    dt[, barcode := NULL]
    dt2 <- data.table(x=c(), N=c())
  } else {
    plot <- gdraw("Unrecognized supplementary .h5 file")
    return(plot)
  }
  # dt: umi pct.intronic called
  # dt2: x N
  
  # Panel 1: cell barcode rank plot
  plotdf <- dt[order(-umi), .(umi,called)]
  plotdf[, index := .I]
  plotdf <- plotdf[plotdf$umi != lag(plotdf$umi,1,0) | plotdf$umi != lead(plotdf$umi,1,0)]
  p1 <- ggplot(plotdf, aes(x=index,y=umi,col=called)) + geom_line() + 
    theme_bw() + scale_x_log10() + scale_y_log10() +
    ggtitle("Cell barcode rank plot") + xlab("Barcode rank") + ylab("UMI count") +
    theme(legend.position = "inside",
          legend.position.inside = c(0.05, 0.05),
          legend.justification.inside = c("left", "bottom"),
          legend.background = element_blank(),
          legend.spacing.y = unit(0.1,"lines"))
  
  # Panel 2: cell calling cloud
  if (is.null(dt$pct.intronic) || all(dt$pct.intronic==0)) {
    p2 <- gdraw(g("No intronic information"))
  } else if (all(dt$umi < ct)) {
    p2 <- gdraw(g("No cells with {ct}+ UMI"))
  } else {
    p2 <- dt[umi >= ct] %>% ggplot(aes(x=log10(umi+1), y=pct.intronic, color=called)) + 
      ggrastr::rasterize(geom_point(alpha = 0.1, size = 0.5), dpi = 300) +
      theme_minimal() +
      labs(title = g("Cell calling"), x = "log10(UMI)", y = "%Intronic") +
      theme(legend.position = "inside",
            legend.position.inside = c(0.95, 0.05),
            legend.justification.inside = c("right", "bottom"),
            legend.background = element_blank(),
            legend.spacing.y = unit(0.1,"lines"))
  }
  
  # Panel 3: Read per UMI histogram
  if (nrow(dt2) > 0) {
    dt2[, x := pmin(x,10)]
    dt2 <- dt2[,.(N=sum(N)),x][order(x)]
    p3 <- ggplot(dt2, aes(x=as.factor(x),y=N)) + geom_col() +
      scale_x_discrete(breaks=min(dt2$x):max(dt2$x), labels=(min(dt2$x):max(dt2$x)) %>% {ifelse(.==10, "10+", .)}) +
      xlab("Reads per UMI") + ylab("UMI Count") + ggtitle("Reads per UMI distribution") +
      theme_minimal()
  } else {
    p3 <- gdraw("No Read per UMI data")
  }
  
  # Panel 4: Intronic density
  max_density_x = dt[called==T, pct.intronic] %>% density %>% {.$x[which.max(.$y)]}
  p4 <- dt[umi>=ct] %>% ggplot(aes(x=pct.intronic, color=called)) + geom_density() + 
    labs(title = g("Intronic density (UMI>={ct})"), x = "%Intronic", y = "Density") + 
    geom_vline(xintercept=max_density_x, color="red", linetype="dashed") +
    annotate(geom='text', label=round(max_density_x, 2), x=max_density_x-1, y=Inf, hjust=1, vjust=1, col="red") +
    theme_minimal() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.05, 0.95),
          legend.justification.inside = c("left", "top"),
          legend.background = element_blank(),
          legend.spacing.y = unit(0.1,"lines"))
  
  plot <- plot_grid(p1, p2, p3, p4, ncol=2)
  return(plot)
}

# Plot UMAP + metrics
plot_umaps <- function(obj) {
  mytheme <- function(){theme(plot.title=element_text(hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="top", legend.justification="center", legend.key.width=unit(2, "lines"))}
  p1 <- DimPlot(obj, label=T) + ggtitle(g("UMAP")) + mytheme() + NoLegend()  + coord_fixed(ratio=1)
  p2 <- VlnPlot(obj,"logumi", alpha=0) + mytheme() + NoLegend() 
  p3 <- FeaturePlot(obj,"percent.mt") + ggtitle("%MT") + mytheme() + coord_fixed(ratio=1) + 
    annotate("text", x = Inf, y = Inf, label = g("Median: {round(median(obj$percent.mt), 2)}%\nMean: {round(mean(obj$percent.mt), 2)}%"), hjust=1, vjust=1, size=2.5, color="black")
  if ("pct.intronic" %in% names(obj@meta.data)) {
    p4 <- FeaturePlot(obj, "pct.intronic") + ggtitle("%Intronic") + mytheme() + coord_fixed(ratio=1) +
      annotate("text", x = Inf, y = Inf, label = g("Median: {round(median(obj$pct.intronic), 2)}%\nMean: {round(mean(obj$pct.intronic), 2)}%"), hjust=1, vjust=1, size=2.5, color="black")
  } else {
    p4 <- gdraw("No intronic information")
  }
  plot <- plot_grid(p1, p2, p3, p4, ncol=2)
  return(plot)
}



















plot_rankplots <- function(df, f, out_path) {
  gdf = count_umis(df)
  
  # p1
  cb.data = gdf %>% group_by(cb_index) %>% summarize(umi=sum(umi)) %>% arrange(desc(umi)) %>% {mutate(.,index=1:nrow(.), filter="all cell barcodes")}
  cb.data2 = cb.data %>% filter(cb_index>0) %>% {mutate(.,index=1:nrow(.), filter="called cell barcodes only")}
  sb_pct_in_called_cells = round(sum(filter(cb.data,cb_index>0)$umi)/sum(cb.data$umi)*100,2) %>% paste0("%")
  p1 = ggplot(mapping=aes(x=index, y=umi,col=filter))+geom_line(data=cb.data)+geom_line(data=cb.data2) +
    scale_x_log10()+scale_y_log10()+theme_bw()+ggtitle("SB UMI per cell")+ylab("SB UMI counts")+xlab("Cell barcodes") +
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"), legend.title=element_blank()) +
    annotate("text", x = Inf, y = Inf, label = g("SB UMI in called cells: {sb_pct_in_called_cells}"), hjust = 1.02, vjust = 1.33)
  rm(cb.data, cb.data2) ; invisible(gc())
  
  # p2
  sb.data = gdf %>% group_by(sb_index) %>% summarize(umi=sum(umi)) %>% arrange(desc(umi)) %>% {mutate(.,index=1:nrow(.),filter="all cell barcodes")}
  sb.data2 = gdf %>% filter(cb_index > 0) %>% group_by(sb_index) %>% summarize(umi=sum(umi)) %>% arrange(desc(umi)) %>% {mutate(.,index=1:nrow(.),filter="called cell barcodes only")}
  p2 = ggplot(mapping=aes(x=index,y=umi,col=filter))+geom_line(data=sb.data)+geom_line(data=sb.data2)+
    scale_x_log10()+scale_y_log10()+theme_bw()+ggtitle("SB UMI per bead")+ylab("SB UMI counts")+xlab("Beads")+
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"), legend.title=element_blank())
  rm(sb.data, sb.data2) ; invisible(gc())
  
  # p3
  x = seq(0, 1, 0.05) * f("metadata/reads")/1000000
  plot.df = data.frame(x=x, y=f("metadata/downsampling")/1000000)
  p3 = ggplot(plot.df, aes(x=x,y=y)) + geom_point() + theme_bw() + 
    xlab("Millions of reads") + ylab("Millions of filtered SB UMIs") + ggtitle("SB downsampling curve")
  
  # p4
  d = table(df$reads) %>% as.data.frame %>% setNames(c("reads","count")) %>% mutate(reads = as.numeric(levels(reads))[reads])
  d %<>% filter(reads < 10) %>% bind_rows(data.frame(reads = 10, count = sum(d$count[d$reads >= 10])))
  total_reads = prettyNum(f('metadata/reads'), big.mark=',')
  sequencing_saturation = round((1 - nrow(df) / sum(df$reads)) * 100, 2) %>% paste0("%")
  p4 = ggplot(d, aes(x=reads, y=count/1000/1000)) + geom_col() +
    theme_bw() + xlab("Reads per UMI") + ylab("Millions of filtered SB UMIs") + ggtitle("SB read depth") + 
    annotate("text", x = Inf, y = Inf, label = g("sequencing saturation = {sequencing_saturation}\ntotal reads = {total_reads}"), hjust = 1.02, vjust = 1.33) +
    scale_x_continuous(breaks=min(d$reads):max(d$reads), labels=(min(d$reads):max(d$reads)) %>% {ifelse(.==10, "10+", .)})
  
  plot = plot_grid(p1, p2, p3, p4, ncol=2)
  
  make.pdf(plot, file.path(out_path, "SB.pdf"), 7, 8)
  
  meta <- metadata
  meta$SB_info$UMI_pct_in_called_cells = sb_pct_in_called_cells
  meta$SB_info$sequencing_saturation = sequencing_saturation
  metadata <<- meta
  return(T)
}

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

beadplot <- function(sb.data) {
  ggplot(sb.data, aes(x=x, y=y, col=umi)) +
    rasterize(geom_point(size=0.1, shape=16), dpi=200) +
    coord_fixed(ratio=1) +
    theme_classic() +
    labs(x="x (\u00B5m)", y="y (\u00B5m)") +
    scale_color_viridis(trans="log", option="B", name="UMI") + 
    ggtitle(g("SB UMI per bead"))
}
plot_beadplot <- function(df, puckdf, out_path) {
  sb.data = df %>% count_umis %>% group_by(sb_index) %>% summarize(umi=sum(umi), .groups="drop") %>%
    inner_join(y=puckdf, by="sb_index") %>% arrange(umi)
  p1 = beadplot(sb.data) + ggtitle(g("SB UMI per bead (total)"))
  sb.data = df %>% filter(cb_index>0) %>% count_umis %>% group_by(sb_index) %>% summarize(umi=sum(umi), .groups="drop") %>%
    inner_join(y=puckdf, by="sb_index") %>% arrange(umi)
  p2 = beadplot(sb.data) + ggtitle(g("SB UMI per bead (called cells only)"))
  plot = plot_grid(p1, p2, ncol=1)
  make.pdf(plot, file.path(out_path, "beadplot.pdf"), 7, 8)
  return(T)
}

plot_metrics <- function(metadata, out_path) {
  
  plot.df = list(
    c("Total Reads", metadata$SB_filtering[["reads_total"]] %>% add.commas),
    c("Final UMIs", metadata$SB_filtering[["UMIs_final"]] %>% add.commas),
    c("R1<->R2", metadata$SB_info$switch_R1R2),
    c("Remap 10X CB", metadata$SB_info$remap_10X_CB),
    c("Bead type", metadata$SB_info$bead_type)
  ) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("Metric", "Value"))
  p_sp = plot_grid(gdraw("Library information"), plot.tab(plot.df), ncol=1, rel_heights=c(0.1,0.6))
  
  header = c("Metric", metadata$puck_info$puck_name %>% str_remove("^Puck_") %>% str_remove("\\.csv$"))
  plot.df = list(c("Beads", metadata$puck_info$num_beads %>% add.commas),
                 c("Size", metadata$puck_info$puck_sizes),
                 c("Filtered beads", Reduce(`+`,metadata$puck_info[c("num_dup","num_N","num_degen","num_lowQ")]) %>% add.commas),
                 c("Scaling factor", metadata$puck_info$scaling_factors %>% round(2)),
                 c("Final UMIs", metadata$puck_info$umi_final %>% add.commas)
  ) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(header)
  p_puck = plot_grid(gdraw("Puck information"),
                     plot.tab(plot.df),
                     ncol=1, rel_heights=c(0.1,0.6))
  
  UP_matching = metadata$UP_matching
  SB_filtering = metadata$SB_filtering
  
  plot.df = data.frame(a=c("exact","fuzzy", "none", "GG"),
                       b=c(UP_matching[["exact"]],
                           UP_matching[["fuzzy"]],
                           UP_matching[["none"]],
                           UP_matching[["GG"]]
                       ) %>% {./sum(.)*100} %>% round(2) %>% paste0("%")
  ) %>% arrange(desc(b)) %>% unname
  p_up = plot_grid(gdraw("UP matching"), plot.tab(plot.df), ncol=1, rel_heights=c(0.1,0.4))
  
  plot.df = data.frame(a=c("exact","fuzzy","none","ambig"),b=metadata$SB_matching[c("exact","HD1","none","HD1ambig")] %>% {./sum(.)*100} %>% round(2) %>% unname %>% paste0("%")) %>% unname
  p_sb = plot_grid(gdraw("SB matching"), plot.tab(plot.df), ncol=1, rel_heights=c(0.1,0.4))
  
  plot.df = data.frame(a=c("exact","fuzzy","none","ambig"),b=metadata$CB_matching[c("exact","HD1","none","HD1ambig")] %>% {./sum(.)*100} %>% round(2) %>% unname %>% paste0("%")) %>% unname
  p_cb = plot_grid(gdraw("CB matching"), plot.tab(plot.df), ncol=1, rel_heights=c(0.1,0.4))
  
  plot.df = list(
    c("Invalid UMI", SB_filtering[["reads_noumi"]]),
    c("No UP", SB_filtering[["reads_noup"]]),
    c("No SB", SB_filtering[["reads_nosb"]]),
    c("Chimeric", SB_filtering[["reads_chimeric"]]),
    c("Invalid SB", SB_filtering[["reads_lowQsb"]]),
    c("Uncalled CB", SB_filtering[["reads_uncalled"]])
  ) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("Filter","Percent"))
  plot.df$Percent = as.numeric(plot.df$Percent) / (SB_filtering[["reads_total"]] - lag(cumsum(plot.df$Percent),1,0))
  plot.df$Percent = round(plot.df$Percent * 100, 2) %>% paste0("%")
  p_filter = plot_grid(gdraw("SB filtering"), plot.tab(plot.df), ncol=1, rel_heights=c(0.1,0.7))
  
  plot.df = metadata$SB_fuzzy_position %>% {data.frame(pos=as.numeric(names(.)),count=as.numeric(unname(.)))}
  switch(metadata$SB_info$bead_type, 
         "V10"={u1=9  ; u2=26},
         "V17"={u1=10 ; u2=27},
         "V15"={u1=16 ; u2=25},
         "V16"={u1=18 ; u2=27}
  )
  plot.df$pos[plot.df$pos>=u1] %<>% add(u2-u1+1)
  p_loc = ggplot(plot.df,aes(x=pos,y=count))+geom_col()+theme_bw() +
    geom_rect(aes(xmin=u1, xmax=u2, ymin=-Inf, ymax=Inf), fill="grey") +
    annotate("text", x=(u1+u2)/2, y=max(plot.df$count, na.rm=T)*0.1, label="UP Site", color="black") +
    xlab("Spatial barcode base position") + ylab("Fuzzy matches") + ggtitle("Location of spatial barcode fuzzy match")
  
  p_R = list(c("R1s", metadata$SB_info$R1s %>% basename %>% str_remove("\\.fastq\\.gz$") %>% str_remove("_R1_001")),
             c("R2s", metadata$SB_info$R2s %>% basename %>% str_remove("\\.fastq\\.gz$") %>% str_remove("_R2_001"))) %>%
    {do.call(rbind,.)} %>% as.data.frame %>% setNames(NULL) %>% plot.tab
  
  plot = plot_grid(
    gdraw("Spatial library metadata",16),
    plot_grid(p_sp, p_puck, ncol=2, rel_widths = c(0.38,0.62)),
    plot_grid(p_up, p_sb, p_cb, ncol=3),
    plot_grid(p_filter, p_loc, ncol=2, rel_widths = c(0.35,0.65)),
    p_R,
    ncol=1,
    rel_heights = c(0.04,0.15,0.11,0.17,0.1)
  )
  
  make.pdf(plot, file.path(out_path, "SBmetrics.pdf"), 7, 8)
  return(T)
}
