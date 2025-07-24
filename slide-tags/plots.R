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
    labs(title="Cell barcode rank plot", x="Barcode rank", y="UMI counts", color=NULL)
  
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
  
  plot <- plot_grid(gdraw(g("Called cells: {sum(dt$called=='Cells')}")), p1, p2,
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
  sb_pct_in_called_cells <- round(nrow(dt[!is.na(cb)])/nrow(dt) * 100, digits=2) %>% paste0("%")
  
  # Panel 1: Spatial barcodes per cell
  cb.data <- rbindlist(list(dt[is.na(cb), .N, cb_raw][,.(N)], dt[!is.na(cb), .N, cb][,.(N)]))[order(-N)]
  cb.data[, index := .I]
  cb.data <- cb.data[N != lag(N,1,0) | N != lead(N,1,0)]
  cb.data[, filter := "all cell barcodes"]
  
  cb.data2 <- dt[!is.na(cb), .N, cb][,.(N)][order(-N)]
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
    annotate("text", x = Inf, y = Inf, label = g("SB UMI in called cells: {sb_pct_in_called_cells}"), hjust = 1.02, vjust = 1.33)
  rm(cb.data, cb.data2) ; invisible(gc())
  
  # Panel 2: Spatial barcodes per bead
  sb.data <- dt[, .N, sb][, .(N)][order(-N)]
  sb.data[, index := .I]
  sb.data <- sb.data[N != lag(N,1,0) | N != lead(N,1,0)]
  sb.data[, filter := "all cell barcodes"]
  
  sb.data2 <- dt[!is.na(cb), .N, sb][, .(N)][order(-N)]
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
  d <- dt[, .N, reads]
  d[, reads := pmin(reads, 10)]
  d <- d[, .(N = sum(N)), reads][order(reads)]
  
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
  ggplot(sb.data, aes(x=x, y=y, col=umi)) +
    ggrastr::rasterize(geom_point(size=0.1, shape=16), dpi=200) +
    coord_fixed(ratio=1) +
    theme_classic() +
    labs(x="x", y="y") +
    scale_color_viridis(trans="log", option="B", name="UMI") + 
    ggtitle(g("SB UMI per bead"))
}
plot_SBplot <- function(dt, puckdf) {
  p1 <- merge(dt[, .(umi=.N), sb], puckdf, by = "sb")[order(umi)] %>% 
    beadplot() + ggtitle(g("SB UMI per bead (total)"))
  
  p2 <- merge(dt[!is.na(cb), .(umi=.N), sb], puckdf, by = "sb")[order(umi)] %>% 
    beadplot() + ggtitle(g("SB UMI per bead (called cells only)"))
  
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
  header = c("Metric", metadata$puck_info$puck_name %>% str_remove("^Puck_") %>% str_remove("\\.csv$"))
  plot.df = list(c("Beads total", metadata$puck_info$num_beads %>% add.commas),
                 c("Beads removed", Reduce(`+`,metadata$puck_info[c("num_dup","num_N","num_degen","num_lowQ")]) %>% add.commas),
                 # c("Size", metadata$puck_info$puck_sizes),
                 c("Diameter", metadata$puck_info$puck_boundaries %>% {map2_dbl(head(.,-1), tail(.,-1), ~round(.y-.x, 2))}),
                 #c("Scaling factor", metadata$puck_info$scaling_factors %>% round(2)),
                 c("Final UMIs", metadata$puck_info$umi_final %>% add.commas)
  ) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(header)
  A2 <- plot_grid(gdraw("Puck information", 13), plot.tab(plot.df), gdraw(""), ncol=1, rel_heights=c(0.1 ,0.5, 0.1))
  
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

plot_gdbscan_opt <- function(coords_global, mranges, k, eps) {
  # Panel 1: DBSCAN cluster distribution
  d <- coords_global[, .(pct=.N/len(cb_whitelist)*100), pmin(clusters, 5)]
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
             label=g("minPts: {minPts}\nplaced: {pct.placed}%\neps: {round(eps,2)} | k: {k}"),
             x=minPts+ms, y=0,
             hjust=0, vjust=0, col="red") + 
    theme_bw() +
    labs(x="minPts", y="%Placed", title="Parameter optimization")
  
  # Panel 3: SB UMI distribution
  p3 <- ggplot(coords_global, aes(x=factor(pmin(clusters,5), levels=0:5, labels=c(0:4,"5+")), y=log10(umi))) + 
    geom_violin(scale="count") + 
    theme_classic() +
    labs(x="DBSCAN clusters", y="log10 SB UMI", title="SB UMI distribution") +
    annotate("text", x=Inf, y=-Inf,
             label=g("Mean SB UMI per cell: {round(coords_global[,sum(umi)/len(cb_whitelist)],1)}"), 
             hjust=1.05, vjust=-0.75, size=3)
  
  # Panel 4: Start vs. End
  # p4 <- ggplot() + geom_density(data=mranges, aes(x=log10((i2+1L)*ms), color=factor("min",levels=c("min","max")))) +
  #                  geom_density(data=mranges, aes(x=log10(i1*ms), color=factor("max",levels=c("min","max")))) + 
  #                  geom_vline(xintercept=log10(minPts), color="blue", linetype="dashed") + 
  #                  theme_bw() +
  #                  labs(x="log10(minPts)", y="Density", title="minPts placement intervals") + 
  #                  scale_color_manual(name=NULL, values=c("min"="green3", "max"="red3")) +
  #                  theme(legend.position = "inside",
  #                        legend.position.inside = c(.99, .99),
  #                        legend.justification.inside = c("right", "top"),
  #                        legend.background = element_blank(),
  #                        legend.spacing.y = unit(0.1,"lines"))
  p4 <- ggplot() + 
    geom_point(data=mranges, size=0.5,
               mapping=aes(x=(i2+1L)*ms,
                           y=(i1-i2)*ms,
                           col=pmin(coords_global$clusters,2) %>% factor(levels=c(0,1,2), labels = c("0","1","2+")))) + 
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
          #legend.margin = margin(t=0, r=0, b=0, l=0),
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

plot_gdbscan_1 <- function(coords_global) {
  coords <- coords_global[clusters==1]
  
  # Panel 1: placements
  p1 <- ggplot(coords, aes(x=x1, y=y1)) + geom_point(size=0.1) +
    coord_fixed(ratio=1) + theme_classic() + 
    xlim(range(puckdf$x)) + ylim(range(puckdf$y)) +
    labs(x="x-position", y="y-position", title="DBSCAN=1 Placements")
  
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
    labs(x="h1", y="Count", title="h-index distribution") + 
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

plot_gdbscan_2 <- function(coords_global, cmes) {
  coords <- coords_global[clusters==2]
  
  # Panel 1: DBSCAN=1 placements
  p1 <- ggplot(coords, aes(x=x1, y=y1)) + geom_point(size=0.1) +
    coord_fixed(ratio=1) + theme_classic() + 
    xlim(range(puckdf$x)) + ylim(range(puckdf$y)) +
    labs(x="x-position", y="y-position", title="DBSCAN=1 Positions")
  
  # Panel 2: DBSCAN=2 placements
  p2 <- ggplot(coords, aes(x=x2, y=y2)) + geom_point(size=0.1) +
    coord_fixed(ratio=1) + theme_classic() + 
    xlim(range(puckdf$x)) + ylim(range(puckdf$y)) +
    labs(x="x-position", y="y-position", title="DBSCAN=2 Positions")
  
  # Panel 3: UMI comparison
  p3 <- ggplot(coords, aes(x=log10(umi1),y=log10(umi2))) + geom_point() +
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

plot_ddbscan_xy <- function(coords_dynamic) {
  d <- coords_dynamic[!is.na(x1) & !is.na(y1)]
  p1 <- ggplot(d, aes(x=x1, y=y1, col=log10(umi1))) + geom_point(size=0.2) + 
    theme_classic() + coord_fixed(ratio=1) + 
    labs(title="1st DBSCAN", x="x", y="y", col="UMI")
  
  d <- coords_dynamic[!is.na(x2) & !is.na(y2)]
  p2 <- ggplot(d, aes(x=x2, y=y2,col=log10(umi2))) + geom_point(size=0.2) + 
    theme_classic() + coord_fixed(ratio=1) + 
    labs(title="2nd DBSCAN", x="x", y="y", col="UMI")
  
  d <- coords_dynamic[, .(r=sqrt((x2-x1)^2+(y2-y1)^2))][is.finite(r)]
  p3 <- ggplot(d, aes(x=r)) + geom_histogram(bins=30) +
    theme_classic() +
    labs(title="Distance between top-2 densities", x="Distance", y="Cells")
  
  plot <- plot_grid(gdraw("Dynamic DBSCAN Positions"),
                    plot_grid(p1, p2, p3, ncol=1),
                    ncol=1, rel_heights=c(0.05,0.95))
  return(plot)
}

# plot_ddbscan_stats <- function(coords_dynamic) {
#   d <- melt(coords_dynamic[,.(sumi1,sumi2)], measure.vars = c("sumi1","sumi2"), variable.name="var", value.name="val")
#   ggplot(d, aes(x=log10(val), col=var)) + geom_density()
# }

# plot_ddbscan_1 <- function() {}
# plot_ddbscan_opt <- function() {}

# plot_dbscan_comp <- function(gcoords, dcoords) {
#   stopifnot(nrow(gcoords) == nrow(dcoords),
#             "cb" %in% names(gcoords), "cb" %in% names(dcoords),
#             gcoords$cb == dcoords$cb)
# }
# 

################################################################################
### Per-cell plots #############################################################
################################################################################

cellplottheme <- theme(plot.title=element_text(hjust=0.5, size=12),
                       legend.background=element_blank(),
                       legend.margin=margin(0,0,0,0),
                       legend.key=element_blank(),
                       legend.key.height=unit(1,"lines"),    # colorbar height
                       legend.key.width=unit(0.5,"lines"),   # colorbar width
                       legend.key.spacing=unit(0.25,"lines"),# colorbar-lab spacing (x/y?)
                       legend.frame=element_blank(),
                       # legend.ticks*, legend.axis*, legend.text*
                       legend.title=element_blank(),
                       legend.position="right",
                       legend.direction="vertical"
                       #legend.byrow*, legend.justification*, legend.location*, legend.box*
)

my_geom_circle <- function(cx, cy, r) {
  geom_path(
    data = data.frame(x = cx + r*cos(seq(0, 2*pi, length.out=100)),
                      y = cy + r*sin(seq(0, 2*pi, length.out=100))),
    aes(x = x, y = y),
    color = "blue",
    linewidth = 0.3,
    linetype = "dashed"
  )
}

global_cellplot <- function(CB) {
  subdf <- data.list[[CB]][!is.na(x) & !is.na(y)][order(umi)]
  subdf1 <- subdf[cluster==1]
  subdf2 <- subdf[cluster==2]
  
  plot <- ggplot() + coord_fixed(ratio=1 ,xlim=range(puckdf$x), ylim=range(puckdf$y)) + theme_void() +
    geom_point(data=subdf, mapping=aes(x=x, y=y, col=umi), size=1, shape=20) +
    labs(title=g("[{CB}]"), col=NULL) +
    cellplottheme
  
  
  # Add the DBSCAN=1 point
  if(nrow(subdf1) > 0) {
    plot <- plot + geom_point(aes(x=weighted.mean(subdf1$x, w=subdf1$umi),
                                  y=weighted.mean(subdf1$y, w=subdf1$umi)),
                              color="red", shape=0, size=2)
  }
  # Add the DBSCAN=2 point
  if(nrow(subdf2) > 0) {
    plot <- plot + geom_point(aes(x=weighted.mean(subdf2$x, w=subdf2$umi),
                                  y=weighted.mean(subdf2$y, w=subdf2$umi)),
                              color="green", shape=0, size=2)
  }
  return(plot)
}

plot_gdbscan_cellplots <- function(data.list) {
  list0 <- map_lgl(data.list, ~max(.$cluster)==0 & nrow(.)>0) %>% {names(.)[.]}
  list1 <- map_lgl(data.list, ~max(.$cluster)==1 & nrow(.)>0) %>% {names(.)[.]}
  list2 <- map_lgl(data.list, ~max(.$cluster)==2 & nrow(.)>0) %>% {names(.)[.]}
  
  if(len(list0) > 0) {
    p0s <- sample(list0, min(12,len(list0)), replace=FALSE) %>% map(global_cellplot)
    p0 <- plot_grid(gdraw("DBSCAN=0"),
                    plot_grid(plotlist=p0s, ncol=3),
                    ncol=1, rel_heights=c(0.1,2))
  } else {p0 <- gdraw("No DBSCAN=0")}
  if(len(list1) > 0) {
    p1s <- sample(list1, min(12,len(list1)), replace=FALSE) %>% map(global_cellplot)
    p1 <- plot_grid(gdraw("DBSCAN=1"),
                    plot_grid(plotlist=p1s, ncol=3),
                    ncol=1, rel_heights=c(0.1,2))
  } else {p1 <- gdraw("No DBSCAN=1")}
  if(len(list2) > 0) {
    p2s <- sample(list2, min(12,len(list2)), replace=FALSE) %>% map(global_cellplot)
    p2 <- plot_grid(gdraw("DBSCAN=2"),
                    plot_grid(plotlist=p2s, ncol=3),
                    ncol=1, rel_heights=c(0.1,2))
  } else {p2 <- gdraw("No DBSCAN=2")}
  
  plots <- list(p0, p1, p2)
  
  return(plots)
}

dynamic_cellplot <- function(CB) {
  dl <- data.list[[CB]][!is.na(x) & !is.na(y)][order(umi)]
  row <- coords_dynamic[cb == CB] %T>% {stopifnot(nrow(.) == 1)}
  
  # Zoomed out
  p1 <- ggplot() + geom_point(data=dl, mapping=aes(x=x, y=y, col=umi)) + 
    geom_point(data=row, mapping=aes(x=x1,y=y1),col="red",shape=0) + 
    geom_point(data=row, mapping=aes(x=x2,y=y2),col="green",shape=0) +
    geom_point(data=row, mapping=aes(x=x,y=y),col="red",shape=5) +
    theme_void() +
    coord_fixed(ratio=1) +
    xlim(range(puckdf$x)) + ylim(range(puckdf$y)) +
    cellplottheme
  
  # Zoomed in
  s <- row[, max(abs(x1-x2), abs(y1-y2))]
  ggplot() + geom_point(data=puckdf, mapping=aes(x=x,y=y), shape=1, col="grey") + 
    geom_point(data=dl, mapping=aes(x=x, y=y, col=umi)) + 
    geom_point(data=dl[clusters2==1], mapping=aes(x=x, y=y), col="red",   shape=1, stroke=0.2, size=1.5) +
    geom_point(data=dl[clusters2==2], mapping=aes(x=x, y=y), col="green", shape=1, stroke=0.2, size=1.5) + 
    geom_point(data=row, mapping=aes(x=x1,y=y1),col="red",shape=0) + 
    geom_point(data=row, mapping=aes(x=x2,y=y2),col="green",shape=0) +
    geom_point(data=row, mapping=aes(x=x,y=y),col="red",shape=5) + 
    my_geom_circle(row$x, row$y, eps) +
    theme_void() +
    coord_fixed(ratio=1) +
    xlim(row[,c(min(x1,x2)-s, max(x1,x2)+s)]) + 
    ylim(row[,c(min(y1,y2)-s, max(y1,y2)+s)]) + 
    cellplottheme
  
  p2 <- ggplot()
  
}

# data.list requires x/y/umi and clusters2
# coords_dynamic requires
plot_dynamic_cellplots <- function(data.list, coords_dynamic) {
  stopifnot(len(data.list) == nrow(coords_dynamic),
            names(data.list) == coords_dynamic$cb)
  
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

# dbscan_vs_kde <- function(coords) {
#   coords %<>% mutate(dist=sqrt((x_um_dbscan-x_um_kde)^2+(y_um_dbscan-y_um_kde)^2))
#   
#   # Panel 1: Distance between DBSCAN and KDE assignments
#   d = filter(coords, !is.na(dist))
#   max_density_x = median(d$dist+1, na.rm=T) %>% log10
#   p1 <- ggplot(d, aes(x=log10(dist+1)))+geom_histogram(bins=30)+theme_bw()+
#     labs(title = "DBSCAN vs KDE distance", x = "log1p Distance (\u00B5m)", y = "Frequency") + 
#     geom_vline(xintercept = max_density_x, color = "red", linetype = "dashed") +
#     annotate(geom = 'text', label = round(10^max_density_x-1, 2) %>% paste0("\u00B5m"), x = max_density_x+0.1, y = Inf, hjust = 0, vjust = 1.3, col="red")
#   
#   # Panel 2: Distribution of KDE ratio for each DBSCAN cluster
#   d = coords %>% rowwise %>% mutate(clusters=min(clusters,5)) %>% ungroup
#   p2 <- ggplot(d, aes(x=clusters %>% as.factor, y=ratio)) + geom_violin() +
#     theme_classic() + xlab("DBSCAN clusters") + ylab("KDE ratio") +
#     scale_x_discrete(breaks=min(d$clusters):max(d$clusters), labels=(min(d$clusters):max(d$clusters)) %>% {ifelse(.==5, "5+", .)}) +
#     ggtitle("KDE ratio per DBSCAN cluster") + 
#     geom_hline(yintercept = 1/3, color = "red", linetype = "dashed")
#   
#   # Panel 3: KDE ratio for disagreeing placements
#   p3 <- coords %>% filter(clusters==1) %>% ggplot(aes(x=dist, y=ratio))+geom_point(size=0.5)+theme_bw()+
#     xlab("Distance between assignments")+ylab("KDE ratio")+ggtitle("Density ratio of disagreements")+
#     geom_hline(yintercept = 1/3, color = "red", linetype = "dashed")
#   
#   # Panel 4: Contingency table of placements
#   d = coords %>% mutate(dbscan_pass=clusters==1, kde_pass=ratio<1/3) %>% group_by(dbscan_pass, kde_pass) %>% summarize(n=n(), .groups="drop") %>% mutate(pct=g("{round(n/sum(n)*100,2)}%\n{n}"))
#   p4 <- ggplot(d, aes(x=dbscan_pass,y=kde_pass,fill=n))+geom_tile()+geom_text(label=d$pct)+theme_bw()+
#     xlab("DBSCAN=1")+ylab("KDE < 1/3")+theme(legend.position="none")+coord_fixed()+ggtitle("Placement table")
#   
#   plot = plot_grid(gdraw("DBSCAN vs. KDE Comparison"),
#                    plot_grid(p1, p2, p3, p4, ncol=2),
#                    ncol=1, rel_heights=c(0.05,0.95))
#   return(plot)
# }
