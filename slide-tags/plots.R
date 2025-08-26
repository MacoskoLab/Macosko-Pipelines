library(ggplot2)
library(ggrastr)

gdraw <- function(text, s=14) {cowplot::ggdraw()+cowplot::draw_label(text, size=s)}
plot.tab <- function(df) {return(cowplot::plot_grid(gridExtra::tableGrob(df, rows=NULL)))}
make.pdf <- function(plots, name, w, h) {
  if (any(c("gg", "ggplot", "Heatmap") %in% class(plots))) {plots = list(plots)}
  grDevices::pdf(file=name, width=w, height=h)
  lapply(plots, print)
  invisible(dev.off())
}

################################################################################
### RNA plots ##################################################################
################################################################################

# Page 1: 10X/Optimus RNA metrics
plot_metrics_csv <- function(df) {
  if (nrow(df) == 0) {
    return(gdraw("No GEX CSV metrics found"))
  }
  
  rownames(df) <- NULL
  colnames(df) <- NULL
  
  # Format decimals as percents, add commas to large numbers
  if (is.numeric(df[,2])) {
    df[,2] %<>% {case_when(
      . >= 0 & . <= 1 ~ dec2pct(.),
      . > 1 ~ add.commas(.),
      TRUE ~ as.character(.)
    )}
  }
  
  # Plot the table
  nslots <- (33/689*768 - 20/531*768)/(33 - 20) * (nrow(df) - 20) + 20/531*768
  hide_cols <- c("NHashID","keeper_cells","percent_keeper","percent_target","percent_usable","keeper_mean_reads_per_cell","keeper_median_genes")
  plot <- plot_grid(gdraw(""),
                    gdraw("GEX Metrics Summary"),
                    plot.tab(df[!df[,1] %in% hide_cols,]),
                    gdraw(""),
                    ncol=1, rel_heights=c((nslots-nrow(df)-1)/2-1, 1+1, nrow(df), (nslots-nrow(df)-1)/2)) #c(0.1,0.1,nrow(df)/11L,0.1))
  return(plot)
}

# Page 2: Cell calling
plot_cellcalling <- function(dt, ct=200) {
  stopifnot(names(dt) == c("umi","pct_intronic","called"))
  
  # Panel 1: cell barcode rank plot
  plotdf <- dt[order(-umi), .(umi,called)]
  plotdf[, index := .I]
  plotdf <- plotdf[plotdf$umi != lag(plotdf$umi,1,0) | plotdf$umi != lead(plotdf$umi,1,0)]
  p1 <- ggplot(plotdf, aes(x=index, y=umi, col=called)) + geom_line() + 
    theme_bw() + scale_x_log10() + scale_y_log10() +
    scale_color_manual(values=c("Cells"="#00BFC4", "Background"="#F8766D"),
                       breaks=c("Cells","Background")) +
    labs(title="Cell barcode rank plot", x="Barcode rank", y="UMI counts", color=NULL) + 
    annotate("text", x=1, y=10,
             label=g("Median UMI Counts per Cell: {dt[called=='Cells', median(umi)] %>% add.commas}"), 
             hjust=0, vjust=0, size=3)
  
  # Panel 2: cell calling cloud
  if (is.null(dt$pct_intronic) || all(dt$pct_intronic==0) || all(is.na(dt$pct_intronic))) {
    p2 <- gdraw(g("No intronic information"))
  } else if (all(dt$umi < ct)) {
    p2 <- gdraw(g("No cells with {ct}+ UMI"))
  } else {
    p2 <- dt[umi >= ct] %>% ggplot(aes(x=log10(umi), y=pct_intronic, color=called)) + 
      ggrastr::rasterize(geom_point(alpha=0.2, size=0.5), dpi=300) +
      theme_bw() +
      scale_color_manual(values=c("Cells"="#00BFC4", "Background"="#F8766D"),
                         breaks=c("Cells","Background")) +
      labs(title="Cell calling", x="log10(UMI)", y="%Intronic", color=NULL)
  }
  
  plot <- plot_grid(gdraw(g("Estimated Number of Cells: {sum(dt$called=='Cells')}")), p1, p2,
                    ncol=1, rel_heights=c(0.1,1,1))
  return(plot)
}

# Page 3: UMAP + QC
plot_umaps <- function(obj) {
  mytheme <- function(){theme(plot.title=element_text(hjust=0.5, size=12),
                              axis.title.x=element_blank(),
                              axis.title.y=element_blank(),
                              axis.text.x=element_text(size=10),
                              axis.text.y=element_text(size=10),
                              legend.position="top",
                              legend.justification="center",
                              legend.key.width=unit(2, "lines"),
                              legend.key.height=unit(0.5, "lines"))}
  p1 <- DimPlot(obj, label=TRUE) + ggtitle(g("UMAP")) + mytheme() + NoLegend() + coord_fixed(ratio=1)
  
  # Violin plots
  obj$log10UMI <- log10(obj$nCount_RNA)
  p2 <- VlnPlot(obj, "log10UMI", alpha=0) + mytheme() + NoLegend() 
  
  # MT% UMAP
  obj$pct.mt <- obj$pct_mt * 100
  p3 <- FeaturePlot(obj, "pct.mt") + ggtitle("%Mito") + mytheme() + coord_fixed(ratio=1) + 
    annotate("text", x = Inf, y = Inf, label = g("Median: {round(median(obj$pct.mt), 2)}%\nMean: {round(mean(obj$pct.mt), 2)}%"), hjust=1, vjust=1, size=2.5, color="black")
  
  # Intronic% UMAP
  if ("pct_intronic" %in% names(obj@meta.data)) {
    obj$pct.intronic <- obj$pct_intronic * 100
    p4 <- FeaturePlot(obj, "pct.intronic") + ggtitle("%Intronic") + mytheme() + coord_fixed(ratio=1) +
      annotate("text", x = Inf, y = Inf, label = g("Median: {round(median(obj$pct.intronic), 2)}%\nMean: {round(mean(obj$pct.intronic), 2)}%"), hjust=1, vjust=1, size=2.5, color="black")
  } else {
    p4 <- gdraw("No intronic information")
  }
  plot <- plot_grid(p1, p2, p3, p4, ncol=2)
  return(plot)
}

################################################################################
### SB plots ###################################################################
################################################################################

# Page 4: Spatial library
plot_SBlibrary <- function(dt, f) {
  sb_pct_in_called_cells <- round(dt[,sum(!is.na(cb))/.N]*100, digits=2) %>% paste0("%")
  
  # Panel 1: Spatial barcodes per cell
  cb.data <- dt[, .N, cr][order(-N), .(N)]
  cb.data[, index := .I]
  cb.data <- cb.data[N != lag(N,1,0) | N != lead(N,1,0)]
  cb.data[, filter := "all cell barcodes"]
  
  cb.data2 <- dt[!is.na(cb), .N, cb][order(-N), .(N)]
  cb.data2[, index := .I]
  cb.data2 <- cb.data2[N != lag(N,1,0) | N != lead(N,1,0)]
  cb.data2[, filter := "called cell barcodes only"]
  
  p1 <- ggplot(mapping=aes(x=index, y=N, col=filter)) + 
               geom_line(data=cb.data) + geom_line(data=cb.data2) +
               scale_x_log10() + scale_y_log10() + theme_bw() + 
               ggtitle("SB UMI per cell") + ylab("SB UMI counts") + xlab("Cell barcode rank") +
               theme(legend.position = "inside",
                     legend.position.inside = c(0.05, 0.05),
                     legend.justification.inside = c("left", "bottom"),
                     legend.background = element_blank(),
                     legend.spacing.y = unit(0.1,"lines"),
                     legend.title=element_blank()) +
               annotate("text", x=Inf, y=Inf, label=g("SB UMI in called cells: {sb_pct_in_called_cells}"), hjust=1.02, vjust=1.33)
  rm(cb.data, cb.data2) ; invisible(gc())
  
  # Panel 2: Spatial barcodes per bead
  sb.data <- dt[, .N, sb][order(-N), .(N)]
  sb.data[, index := .I]
  sb.data <- sb.data[N != lag(N,1,0) | N != lead(N,1,0)]
  sb.data[, filter := "all cell barcodes"]
  
  sb.data2 <- dt[!is.na(cb), .N, sb][order(-N), .(N)]
  sb.data2[, index := .I]
  sb.data2 <- sb.data2[N != lag(N,1,0) | N != lead(N,1,0)]
  sb.data2[, filter := "called cell barcodes only"]
  
  p2 <- ggplot(mapping=aes(x=index, y=N, col=filter)) + 
               geom_line(data=sb.data) + geom_line(data=sb.data2) +
               scale_x_log10() + scale_y_log10() + theme_bw() + 
               ggtitle("SB UMI per bead") + ylab("SB UMI counts") + xlab("Beads") +
               theme(legend.position = "inside",
                     legend.position.inside = c(0.05, 0.05),
                     legend.justification.inside = c("left", "bottom"),
                     legend.background = element_blank(),
                     legend.spacing.y = unit(0.1,"lines"),
                     legend.title=element_blank())
  rm(sb.data, sb.data2) ; invisible(gc())
  
  # Panel 3: Spatial barcode library downsampling curve
  p3 <- data.frame(x = seq(0, 1, 0.05) * f("metadata/reads")/1000000,
                   y = f("metadata/downsampling")/1000000) %>% 
        ggplot(aes(x=x,y=y)) + geom_point() + theme_bw() + 
        xlab("Millions of SB reads") + ylab("Millions of filtered SB UMIs") + ggtitle("SB downsampling curve")
  
  # Panel 4: Reads per UMI distribution
  d <- dt[, .N, .(reads=pmin(reads, 10))][order(reads)]
  total_reads = add.commas(f('metadata/reads'))
  sequencing_saturation = round((1 - nrow(dt) / sum(dt$reads)) * 100, 2) %>% paste0("%")
  
  p4 <- ggplot(d, aes(x=reads, y=N/1000/1000)) + geom_col() +
    theme_bw() + xlab("Reads per UMI") + ylab("Millions of filtered SB UMIs") + ggtitle("Reads per UMI") + 
    annotate("text", x = Inf, y = Inf, label = g("sequencing saturation = {sequencing_saturation}\ntotal reads = {total_reads}"), hjust = 1.02, vjust = 1.33) +
    scale_x_continuous(breaks=min(d$reads):max(d$reads), labels=(min(d$reads):max(d$reads)) %>% {ifelse(.==10, "10+", .)})
  
  return(plot_grid(p1, p2, p3, p4, ncol=2))
}

# Page 5: Bead plots
beadplot <- function(sb.data) {
  sb.data[order(umi)] %>% ggplot(aes(x=x, y=y, col=umi)) +
    ggrastr::rasterize(geom_point(size=0.1, shape=16), dpi=200) +
    coord_fixed(ratio=1) +
    theme_classic() +
    labs(x="", y="") +
    scale_color_viridis(trans="log", option="B", name="UMI") +
    ggtitle(g("SB UMI per bead"))
}
plot_SBplot <- function(dt, puckdf) {
  sbd <- dt[, .(umi=.N), sb] %>% merge(puckdf, by="sb", all=FALSE)
  p1 <- beadplot(sbd) + ggtitle(g("SB UMI per bead (total)"))
  
  sbd <- dt[!is.na(cb), .(umi=.N), sb] %>% merge(puckdf, by="sb", all=FALSE)
  p2 <- beadplot(sbd) + ggtitle(g("SB UMI per bead (called cells only)"))

  return(plot_grid(p1, p2, ncol=1))
}

# Page 6: Spatial library metrics
plot_SBmetrics <- function(metadata) {
  
  # Panel A1
  plot.df = list(
    c("Total Reads", metadata$SB_filtering[["reads_total"]] %>% add.commas),
    c("Final UMIs", metadata$SB_filtering[["UMIs_final"]] %>% add.commas),
    c("R1 <-> R2", metadata$SB_info$switch_R1R2),
    c("Remap CB", metadata$SB_info$remap_10X_CB),
    c("Bead type", metadata$SB_info$bead_type)
  ) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("Metric", "Value"))
  A1 <- plot_grid(gdraw("Library information", 13), plot.tab(plot.df), ncol=1, rel_heights=c(0.1, 0.6))
  
  # Panel A2
  header = c("Metric", metadata$puck_info$puck_name %>% trim_puck_name)
  plot.df = list(c("Beads total", metadata$puck_info$num_beads %>% add.commas),
                 c("Beads removed", Reduce(`+`,metadata$puck_info[c("num_dup","num_N","num_degen","num_lowQ")]) %>% add.commas),
                 # c("Size", metadata$puck_info$puck_sizes),
                 c("Diameter", metadata$puck_info$puck_boundaries %>% {map2_dbl(head(.,-1), tail(.,-1), ~round(.y-.x, 2))}),
                 #c("Scaling factor", metadata$puck_info$scaling_factors %>% round(2)),
                 c("Final UMIs", metadata$puck_info$umi_final %>% add.commas)
  ) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(header)
  A2 <- plot_grid(gdraw("Puck information", 13), plot.tab(plot.df), ncol=1, rel_heights=c(0.1 ,0.6))
  
  # Panel B1
  plot.df = data.frame(a=c("exact", "fuzzy", "GG", "none"),
                       b=metadata$UP_matching[c("exact", "fuzzy", "GG", "none")]) %>%
    mutate(b=(b/sum(b)*100) %>% round(2) %>% paste0("%")) %>% unname
  B1 <- plot_grid(gdraw("UP matching", 13), plot.tab(plot.df), ncol=1, rel_heights=c(0.1, 0.4))
  
  # Panel B2
  plot.df = data.frame(a=c("exact","fuzzy","ambig","none"),
                       b=metadata$SB_matching[c("exact","HD1","HD1ambig","none")] %>% {./sum(.)*100} %>% round(2) %>% unname %>% paste0("%")) %>% unname
  B2 <- plot_grid(gdraw("SB matching", 13), plot.tab(plot.df), ncol=1, rel_heights=c(0.1, 0.4))
  
  # Panel B3
  plot.df = data.frame(a=c("exact","fuzzy","ambig","none"),
                       b=metadata$CB_matching[c("exact","fuzzy","ambig","none")] %>% {./sum(.)*100} %>% round(2) %>% unname %>% paste0("%")) %>% unname
  B3 <- plot_grid(gdraw("CB matching", 13), plot.tab(plot.df), ncol=1, rel_heights=c(0.1, 0.4))
  
  # Panel C1
  plot.df = data.frame(a=c("LQ UMI", "No UP", "No SB", "LQ SB", "Uncalled CB", "Chimeric"),
                       b=metadata$SB_filtering[c("reads_noumi", "reads_noup", "reads_nosb", "reads_lqsb", "reads_nocb", "reads_chimeric")] %>% unname) %>% 
    setNames(c("Filter","Percent"))
  plot.df$Percent = as.numeric(plot.df$Percent) / (metadata$SB_filtering[["reads_total"]] - lag(cumsum(plot.df$Percent),1,0))
  plot.df$Percent = round(plot.df$Percent * 100, 2) %>% paste0("%")
  C1 <- plot_grid(gdraw("SB filtering", 13), plot.tab(plot.df), ncol=1, rel_heights=c(0.1,0.7))
  
  # Panel C2
  if (metadata$SB_info$bead_type %in% c("V10","V17","V15","V16")) {
    plot.df = metadata$SB_fuzzy_position %>% {data.frame(pos=as.numeric(names(.)),count=as.numeric(unname(.)))}
    switch(metadata$SB_info$bead_type, 
           "V10"={u1=9  ; u2=26},
           "V17"={u1=10 ; u2=27},
           "V19"={u1=10 ; u2=19},
           "V15"={u1=16 ; u2=25},
           "V16"={u1=18 ; u2=27}
    )
    plot.df$pos[plot.df$pos>=u1] %<>% add(u2-u1+1)
    C2 <- ggplot(plot.df, aes(x=pos,y=count)) + geom_col() + theme_bw() +
      geom_rect(aes(xmin=u1, xmax=u2, ymin=-Inf, ymax=Inf), fill="grey") +
      annotate("text", x=(u1+u2)/2, y=max(plot.df$count, na.rm=T)*0.1, label="UP Site", color="black") +
      xlab("Spatial barcode base position") + ylab("Fuzzy matches") + ggtitle("Location of spatial barcode fuzzy match")
  } else {
    C2 <- gdraw("Unrecognized bead type")
  }
  
  # Panel D
  plot.df <- data.frame(R1s=metadata$SB_info$R1s %>% basename %>% str_remove("_R1_001.fastq.gz$"),
                        R2s=metadata$SB_info$R2s %>% basename %>% str_remove("_R2_001.fastq.gz$"))
  if (nrow(plot.df) > 4) {
    plot.df <- rbind(plot.df[1:3,], data.frame(R1s=g("({nrow(plot.df)-3} more)"), R2s=g("({nrow(plot.df)-3} more)")))
  }
  D1 <- plot_grid(gdraw("Spatial library FASTQs", 13), plot.tab(plot.df), ncol=1, rel_heights=c(0.1, 0.5))
  
  plot <- plot_grid(
    gdraw("Spatial Metrics Summary", 15),
    plot_grid(A1, A2, ncol=2, rel_widths = c(0.38, 0.62)),
    plot_grid(B1, B2, B3, ncol=3),
    plot_grid(C1, C2, ncol=2, rel_widths = c(0.35, 0.65)),
    D1,
    ncol=1,
    rel_heights = c(0.1*1.1, 0.7, 0.5, 0.8, 0.6)
  )
  
  return(plot)
}

################################################################################
### DBSCAN plots ###############################################################
################################################################################

plot_dbscan_opt <- function(coords, mranges) {
  # Panel 1: DBSCAN cluster distribution
  d <- coords[, .(pct=.N/len(cb_whitelist)*100), pmin(clusters, 5)]
  d[, pmin := factor(pmin, levels=0:5, labels=c(0:4, "5+"))]
  p1 <- ggplot(d, aes(x=pmin, y=pct)) + geom_col() +
    geom_text(aes(label=sprintf("%1.0f%%",pct), y=pct), vjust=-0.5) +
    theme_classic() +
    scale_y_continuous(limits=c(0,100)) +
    labs(x="DBSCAN clusters", y="Percent", title="Cluster distribution")
  
  # Panel 2: minPts search
  d <- IRanges::IRanges(start=mranges$i2+1L, end=mranges$i1) %>% IRanges::coverage() %>% as.integer %>% {data.table(i=seq_along(.), p=.)}
  i_opt <- d[p==max(p), tail(i, 1)]
  minPts <- i_opt * ms
  pct.placed <- d[, round(max(p)/len(cb_whitelist)*100, 2)]
  d <- d[i >= 2 & i <= max(20, i_opt*2)]
  p2 <- ggplot(d, aes(x=i*ms, y=p/len(cb_whitelist)*100)) + geom_line() +
    geom_vline(xintercept=minPts, color="red", linetype="dashed") +
    annotate(geom='text',
             label=g("minPts: {minPts}\nplaced: {pct.placed}%\neps: {round(eps,2)} | k: {knn}"),
             x=minPts+ms, y=0,
             hjust=0, vjust=0, col="red") + 
    theme_bw() +
    labs(x="minPts", y="%Placed", title="Parameter optimization")
  
  # Panel 3: SB UMI distribution
  p3 <- ggplot(coords, aes(x=factor(pmin(clusters,5), levels=0:5, labels=c(0:4,"5+")), y=log10(umi))) + 
    geom_violin(scale="count") + 
    theme_classic() +
    labs(x="DBSCAN clusters", y="log10 SB UMI", title="SB UMI distribution") +
    annotate("text", x=Inf, y=-Inf,
             label=g("Mean SB UMI per cell: {round(coords[,sum(umi)/len(cb_whitelist)],1)}"), 
             hjust=1.05, vjust=-0.75, size=3)
  
  # Panel 4: Start vs. Width
  p4 <- ggplot() + 
    geom_point(data=mranges, size=0.1,
               mapping=aes(x=(i2+1L)*ms,
                           y=(i1-i2)*ms,
                           col=pmin(coords$clusters,2) %>% factor(levels=c(0,1,2), labels = c("0","1","2+")))) + 
    geom_vline(xintercept = minPts) +
    geom_line(data=data.table(x=seq_len(minPts))[,y:=minPts-x+1], mapping=aes(x=x,y=y)) +
    theme_bw() +
    scale_x_continuous(trans='log10') + scale_y_continuous(trans='log10') +
    labs(x="Minimum minPts", y="Placement Interval Width", col="DBSCAN") +
    scale_color_manual(values = c("#619CFF", "#F8766D", "#00BA38")) +
    theme(legend.position = "inside",
          legend.position.inside = c(0,1),
          legend.justification.inside = c("left", "top"),
          legend.box.background = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.key.height = unit(0.5, "lines"),
          legend.key.width = unit(0.5, "lines"),
          legend.spacing.y = unit(0.0, "lines"),
          legend.text = element_text(size=7),
          legend.title = element_text(size=8))
  
  plot <- plot_grid(gdraw("Global DBSCAN Optimization"),
                    plot_grid(p1, p2, p3, p4, ncol=2),
                    ncol=1, rel_heights=c(0.05,0.95)) %>% suppressWarnings
  return(plot)
}

plot_dbscan_1 <- function(coords) {
  coords <- coords[clusters==1]
  
  # Panel 1: placements
  p1 <- ggplot(coords, aes(x=x1, y=y1)) + geom_point(size=0.1) +
    coord_fixed(ratio=1) + theme_classic() + 
    xlim(range(puckdf$x)) + ylim(range(puckdf$y)) +
    labs(x=NULL, y=NULL, title="DBSCAN=1 Placements") + guides(col="none")
  
  # Helper
  mydensplot <- function(coords, var, name, percent) {
    stopifnot(var %in% names(coords))
    max_density_x = mean(coords[[var]], na.rm=TRUE)
    ggplot(coords, aes(x = .data[[var]])) + geom_density() + 
      theme_minimal() +
      labs(x=name, y="Density", title=g("{name} distribution")) + 
      geom_vline(xintercept=max_density_x, color="red", linetype="dashed") +
      annotate(geom='text',
               label=ifelse(percent, g(" Mean: {round(max_density_x*100, 1)}%"), g(" Mean: {round(max_density_x, 2)}")),
               x=max_density_x, y=Inf,
               hjust=0, vjust=1, col="red")
  }
  
  # Panel 2: SNR
  coords[, SNR := umi1/umi]
  p2 <- mydensplot(coords, "SNR", "SNR", TRUE)
  
  # Panel 3: RMS
  p3 <- mydensplot(coords, "rmsd1", "RMSD", FALSE)
  
  # Panel 4: h-index
  p4 <- ggplot(coords, aes(x=h1)) + geom_bar() + theme_minimal() +
    labs(x="h-index", y="Count", title="h-index distribution") + 
    geom_vline(xintercept=mean(coords$h1), color="red", linetype="dashed") +
    annotate(geom='text',
             label=g(" Mean: {round(mean(coords$h1), 2)}"),
             x=mean(coords$h1), y=Inf,
             hjust=0, vjust=1, col="red")
  
  plot <- plot_grid(gdraw("Global DBSCAN=1 Results"),
                    plot_grid(p1, p2, p3, p4, ncol=2),
                    ncol=1, rel_heights=c(0.05,0.95))
  return(plot)
}

plot_dbscan_2 <- function(coords, cmes) {
  coords <- coords[clusters==2]
  
  # Panel 1: DBSCAN=1 placements
  p1 <- ggplot(coords, aes(x=x1, y=y1)) + geom_point(size=0.1) +
    coord_fixed(ratio=1) + theme_classic() + 
    xlim(range(puckdf$x)) + ylim(range(puckdf$y)) +
    labs(x=NULL, y=NULL, title="DBSCAN=1/2 Positions")
  
  # Panel 2: DBSCAN=2 placements
  p2 <- ggplot(coords, aes(x=x2, y=y2)) + geom_point(size=0.1) +
    coord_fixed(ratio=1) + theme_classic() + 
    xlim(range(puckdf$x)) + ylim(range(puckdf$y)) +
    labs(x=NULL, y=NULL, title="DBSCAN=2/2 Positions")
  
  # Panel 3: UMI comparison
  p3 <- ggplot(coords, aes(x=log10(umi1),y=log10(umi2))) + geom_point(size=0.1) +
    coord_fixed(ratio=1) + theme_bw() + 
    labs(x="log10(UMI1)", y="log10(UMI2)", title="SB UMI counts")
  
  # Panel 4: DBSCAN 1-2 Distance
  p4 <- ggplot(coords, aes(x=sqrt((x1-x2)^2+(y1-y2)^2))) + geom_histogram(bins=30) + 
    theme_classic() +
    labs(x="Distance", y="Cells", title="DBSCAN 1-2 Distance") +
    geom_vline(xintercept=cmes*eps, color="red", linetype="dashed") +
    annotate(geom='text',
             label=g(" CMES: {round(cmes, 2)}"),
             x=cmes*eps, y=Inf,
             hjust=0, vjust=1, col="red")
  
  plot <- plot_grid(gdraw("Global DBSCAN=2 Results"),
                    plot_grid(p1, p2, p3, p4, ncol=2),
                    ncol=1, rel_heights=c(0.05,0.95))
  return(plot)
}

plot_dbscan_score <- function(coords2) {
  p1 <- ggplot(coords2[!is.na(x) & !is.na(y)], aes(x=x,y=y,col=score)) + geom_point(size=0.1) + 
    coord_fixed(ratio=1) + theme_classic() + 
    xlim(range(puckdf$x)) + ylim(range(puckdf$y)) + 
    labs(x=NULL, y=NULL, title="Scored DBSCAN Positions") +
    cellplottheme + theme(plot.title=element_text(size=12)) +
    scale_color_continuous(limits = c(0, 1))
  
  median_score <- coords2[,median(score, na.rm=TRUE)]
  p2 <- ggplot(coords2, aes(x=score)) + geom_histogram(bins=50) +
    theme_bw() +
    geom_vline(xintercept=median_score, color="red", linetype="dashed") +
    annotate(geom='text', label=g("Median: {round(median_score,2)}"),
             x=median_score, y=Inf, hjust=1.05, vjust=2, col="red") + 
    labs(x="Score", y="Count", title="Score distribution")
  
  p3 <- ggplot(coords2) + theme_bw() + 
    geom_density(aes(x=log10(umi1), color="DBSCAN=1"), bw=0.1, key_glyph="path") +
    geom_density(aes(x=log10(umi2), color="DBSCAN=2"), bw=0.1, key_glyph="path") + 
    scale_color_manual(values = c("DBSCAN=1"="red", "DBSCAN=2"="green")) +
    labs(x="log10(UMI)", y="Density", title="DBSCAN UMI Density") + 
    theme(legend.position = "inside",
          legend.position.inside = c(0.99, 0.99),
          legend.justification.inside = c("right", "top"),
          legend.background = element_blank(),
          legend.spacing.y = unit(0.1,"lines"),
          legend.title=element_blank())
  
  xline <- seq(max(coords2$umi1, na.rm = TRUE))
  yline <- quantile(coords2$umi1, probs=0.1*F1(xline), type=1, names=FALSE, na.rm = TRUE)
  p4 <- ggplot() + geom_point(aes(x=log10(umi1),y=log10(umi2),col=score), coords2, size=0.1) + 
    #geom_path(aes(x=log10(xline), y=log10(yline)),col="red", linewidth=0.2) +
    theme_bw() + scale_color_continuous(limits = c(0, 1)) + 
    labs(x="log10(UMI1)", y="log10(UMI2)", title="DBSCAN UMI Distribution") + 
    theme(legend.position = "inside",
          legend.position.inside = c(0.01, 0.99),
          legend.justification.inside = c("left", "top"),
          legend.background = element_blank(),
          legend.spacing.y = unit(0.1,"lines"))
  
  plot <- plot_grid(gdraw("Dynamic DBSCAN Results"),
                    plot_grid(p1, p2, p3, p4, ncol=2),
                    ncol=1, rel_heights=c(0.05,0.95))
  
  return(plot)
}

################################################################################
### Per-cell plots #############################################################
################################################################################

cellplottheme <- theme(plot.title=element_text(hjust=0.5, size=10),
                       legend.background=element_blank(),
                       legend.margin=margin(0,0,0,0),
                       legend.key=element_blank(),
                       legend.key.height=unit(4/3, "lines"),
                       legend.key.width=unit(0.5,"lines"),
                       legend.key.spacing=unit(0.25,"lines"), # .x, .y
                       legend.frame=element_blank(),
                       # legend.ticks*, legend.axis*, legend.text*
                       legend.title=element_blank(),
                       legend.position="right",
                       legend.direction="vertical"
                       #legend.byrow*, legend.justification*, legend.location*, legend.box*
)

cellplot <- function(CB) {
  subdf <- data.list[[CB]][order(umi)]
  subdf1 <- subdf[cluster==1]
  subdf2 <- subdf[cluster==2]
  
  plot <- ggplot() + coord_fixed(ratio=1, xlim=range(puckdf$x), ylim=range(puckdf$y)) + theme_void() +
    geom_point(data=subdf, mapping=aes(x=x, y=y, col=umi), size=1, shape=20) +
    labs(title=g("[{CB}]"), col=NULL) +
    cellplottheme
  
  # Add the DBSCAN=1 point
  if(nrow(subdf1) > 0) {
    plot <- plot + geom_point(aes(x=weighted.mean(subdf1$x, w=subdf1$umi),
                                  y=weighted.mean(subdf1$y, w=subdf1$umi)),
                              color="red", shape=0, size=1.5)
  }
  # Add the DBSCAN=2 point
  if(nrow(subdf2) > 0) {
    plot <- plot + geom_point(aes(x=weighted.mean(subdf2$x, w=subdf2$umi),
                                  y=weighted.mean(subdf2$y, w=subdf2$umi)),
                              color="green", shape=0, size=1.5)
  }
  return(plot)
}

plot_dbscan_cellplots <- function(data.list) {
  list0 <- map_lgl(data.list, ~nrow(.)>0 && max(.$cluster)==0) %>% {names(.)[.]}
  list1 <- map_lgl(data.list, ~nrow(.)>0 && max(.$cluster)==1) %>% {names(.)[.]}
  list2 <- map_lgl(data.list, ~nrow(.)>0 && max(.$cluster)==2) %>% {names(.)[.]}
  
  if(len(list0) > 0) {
    p0s <- sample(list0, min(12,len(list0)), replace=FALSE) %>% map(cellplot)
    p0 <- plot_grid(gdraw("DBSCAN=0"),
                    plot_grid(plotlist=p0s, ncol=3),
                    ncol=1, rel_heights=c(0.1,2))
  } else {p0 <- gdraw("No DBSCAN=0")}
  if(len(list1) > 0) {
    p1s <- sample(list1, min(12,len(list1)), replace=FALSE) %>% map(cellplot)
    p1 <- plot_grid(gdraw("DBSCAN=1"),
                    plot_grid(plotlist=p1s, ncol=3),
                    ncol=1, rel_heights=c(0.1,2))
  } else {p1 <- gdraw("No DBSCAN=1")}
  if(len(list2) > 0) {
    p2s <- sample(list2, min(12,len(list2)), replace=FALSE) %>% map(cellplot)
    p2 <- plot_grid(gdraw("DBSCAN=2"),
                    plot_grid(plotlist=p2s, ncol=3),
                    ncol=1, rel_heights=c(0.1,2))
  } else {p2 <- gdraw("No DBSCAN=2")}
  
  plots <- list(p0, p1, p2)
  
  return(plots)
}

debug.plot <- function(dl) {
  stopifnot(c("x","y","umi","cluster2") %in% names(dl))
  p1 <- ggplot() + theme_void() + coord_fixed(ratio=1) + cellplottheme +
    geom_point(aes(x=x,y=y,col=umi), dl[order(umi)], size=0.5) +
    geom_point(aes(x=x,y=y), dl[cluster2==1, .(x=weighted.mean(x,umi), y=weighted.mean(y,umi))], color="red", shape=0) + 
    geom_point(aes(x=x,y=y), dl[cluster2==2, .(x=weighted.mean(x,umi), y=weighted.mean(y,umi))], color="green", shape=0)
  p2 <- ggplot() + theme_void() + coord_fixed(ratio=1) + guides(col="none") +
    geom_point(aes(x=x,y=y), dl[cluster2==0], color="grey", size=0.5) + 
    geom_point(aes(x=x,y=y,col=as.factor(cluster2)), dl[cluster2>0], size=0.5)
  
  p12 <- plot_grid(p1, p2, ncol=2)
  p1$layers[[1]]$aes_params$size <- 2
  p1$layers[[2]]$aes_params$size <- 3
  p1$layers[[3]]$aes_params$size <- 3
  p2$layers[[1]]$aes_params$size <- 2
  p2$layers[[2]]$aes_params$size <- 2
  
  # DBSCAN=1
  xlims <- dl[cluster2==1,range(x)] %>% {2*(.-mean(.))+mean(.)}
  ylims <- dl[cluster2==1,range(y)] %>% {2*(.-mean(.))+mean(.)}
  p3 <- p1 + xlim(xlims) + ylim(ylims)
  p4 <- p2 + xlim(xlims) + ylim(ylims) +
    annotate("text", x=mean(xlims), y=ylims[[2]], label=g("UMI: {dl[cluster2==1,sum(umi)]}"))
  
  # DBSCAN=2
  xlims <- dl[cluster2==2,range(x)] %>% {2*(.-mean(.))+mean(.)}
  ylims <- dl[cluster2==2,range(y)] %>% {2*(.-mean(.))+mean(.)}
  p5 <- p1 + xlim(xlims) + ylim(ylims)
  p6 <- p2 + xlim(xlims) + ylim(ylims) +
    annotate("text", x=mean(xlims), y=ylims[[2]], label=g("UMI: {dl[cluster2==2,sum(umi)]}"))
  
  return(suppressWarnings(plot_grid(p12,
                                    plot_grid(p3,p4,p5,p6,ncol=2),
                                    ncol=1, rel_heights=c(1,2))))
}

plot_debug_cellplots <- function(data.list, coords2) {
  coords2[!is.na(score)][order(-score)][round(seq(0.1,0.9,0.1)*.N), .(cb, score)] %>% 
    pmap(function(cb, score) {
      plot_grid(plot_grid(gdraw(g("[{cb}]")), plot_grid(gdraw(g("({score})")), ncol=2)),
                debug.plot(data.list[[cb]]),
                ncol=1, rel_heights=c(0.05,0.95))
    })
}

################################################################################
### RNA+Spatial plots ##########################################################
################################################################################

# Create DimPlot
plot_clusters <- function(obj, reduction) {
  npucks = (max(obj$x, na.rm=TRUE) - min(obj$x, na.rm=TRUE))/(max(obj$y, na.rm=TRUE) - min(obj$y, na.rm=TRUE))
  nclusters = len(unique(obj$seurat_clusters))
  ncols = round(sqrt(npucks*nclusters/2)/npucks*2) 
  
  m = obj@reductions[[reduction]]@cell.embeddings %>% {!is.na(.[,1]) & !is.na(.[,2])}
  title = g("%placed: {round(sum(m)/len(m)*100,2)} ({sum(m)}/{len(m)}) [{reduction}]")
  p1 = DimPlot(obj, reduction=reduction) + coord_fixed(ratio=1) +
    ggtitle(title) + NoLegend() + xlab("x-position") + ylab("y-position") + 
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
  obj$placed <- !is.na(obj$x) & !is.na(obj$y)
  
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
