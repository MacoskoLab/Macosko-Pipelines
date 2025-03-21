library(glue) ; g=glue ; len=length
library(dplyr)
library(purrr)
library(rlist)
library(rhdf5)
library(stringr)
library(ggplot2)
library(ggrastr)
library(cowplot)
library(viridis)
library(magrittr)
library(jsonlite)
library(gridExtra)

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
if (!dir.exists(out_path)) { dir.create(out_path, recursive = T) }
stopifnot(system(g("h5ls {sb_path}/matrix"), intern=T) %>% strsplit(split = "\\s+") %>% map_chr(pluck(1)) == c("cb_index", "reads", "sb_index", "umi"))
f <- function(p){return(h5read(sb_path, p))}

# load the CB whitelist
cb_whitelist = readLines(cb_path)
stopifnot(class(cb_whitelist) == "character")
stopifnot(!duplicated(cb_whitelist))
stopifnot(len(unique(nchar(cb_whitelist))) == 1)
stopifnot(map_lgl(strsplit(cb_whitelist, ""), ~all(. %in% c("A","C","G","T"))))
print(g("{len(cb_whitelist)} cell barcodes loaded"))

# Load the SB count matrix
df = data.frame(cb_index=f("matrix/cb_index"),
                umi_2bit=f("matrix/umi"),
                sb_index=f("matrix/sb_index"),
                reads=f("matrix/reads"))
print(g("{sum(df$reads)} spatial barcode reads loaded"))

# load metadata
metadata = list(SB_info = list(R1s=f("lists/R1_list"),
                               R2s=f("lists/R2_list"),
                               pucks=f("lists/puck_list"),
                               UMI_downsampling=f("metadata/downsampling"),
                               switch_R1R2=f("metadata/switch") %>% as.logical,
                               bead_type=f("metadata/bead")),
                UMI_filtering = setNames(f("metadata/UMI/count"), f("metadata/UMI/type")),
                UP_matching = setNames(f("metadata/UP/count"), f("metadata/UP/type")),
                SB_matching = setNames(f("metadata/SB/count"),f("metadata/SB/type")),
                SB_fuzzy_position = setNames(f("metadata/SB_HD/count"), f("metadata/SB_HD/type")) %>% {.[order(as.integer(names(.)))]},
                puck_info = list(num_lowQ=f("metadata/num_lowQbeads"))
)
metadata$SB_filtering = c(reads_total=f("metadata/reads"),
                          reads_tooshort=f("metadata/R1_tooshort")+f("metadata/R2_tooshort"),
                          reads_noumi=sum(metadata$UMI[c("N","homopolymer")]),
                          reads_noup=sum(metadata$UP[c("none","GG")]),
                          reads_nosb=sum(metadata$SB_matching[c("none","HD1ambig")]))

### Helper methods #############################################################

remap_10X_CB <- function(vec) {
  stopifnot(class(vec) == "character")
  basemap = setNames(c("AG","TC","CA","GT"), c("TC","AG","GT","CA"))
  stopifnot(substr(vec,8,9) %in% names(basemap))
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

gdraw <- function(text, s=13) {ggdraw()+draw_label(text, size=s)}
add.commas <- function(num){prettyNum(num, big.mark=",")}
plot.tab <- function(df) {return(plot_grid(tableGrob(df,rows=NULL)))}
make.pdf <- function(plots, name, w, h) {
  if ("gg" %in% class(plots) || class(plots)=="Heatmap") {plots = list(plots)}
  pdf(file=name, width=w, height=h)
  lapply(plots, function(x){print(x)})
  dev.off()
}

### Workflow ###################################################################

# Convert cb_index from a cb_list index to a cb_whitelist index
print("Performing fuzzy cell-barcode matching")
fuzzy_matching <- function(df, cb_list, cb_whitelist) {
  # Check the lists
  stopifnot(!any(duplicated(cb_list)))
  stopifnot(!any(duplicated(cb_whitelist)))
  stopifnot(names(df) == c("cb_index","umi_2bit","sb_index","reads"))
  
  # Remap cb_whitelist
  cbs = cb_list[df$cb_index]
  reads_noremap = df$reads[cbs %in% cb_whitelist] %>% sum
  reads_remap = df$reads[cbs %in% remap_10X_CB(cb_whitelist)] %>% sum
  remap = reads_remap > reads_noremap ; rm(cbs)
  if (remap) {
    print("Remapping CB whitelist")
    cb_whitelist %<>% remap_10X_CB
  }
  stopifnot(!duplicated(cb_whitelist))
  
  # Exact matching dictionary
  exact_dict = match(cb_list, cb_whitelist)
  
  # HD1 fuzzy matching dictionary
  dict = imap(cb_whitelist, ~data.frame(original=.y, neighbor=listHD1neighbors(.x))) %>% {do.call(rbind,.)}
  dict %<>% filter(neighbor %in% cb_list) # remove barcodes that don't exist
  dict %<>% filter(!neighbor %in% cb_whitelist) # remove exact matches
  HD1ambig = unique(dict$neighbor[duplicated(dict$neighbor)])
  dict %<>% filter(!neighbor %in% HD1ambig) # remove ambiguous matches
  stopifnot(!any(duplicated(dict$neighbor)))
  fuzzy_dict = dict$original[match(cb_list, dict$neighbor)]
  HD1ambig_dict = !is.na(match(cb_list, HD1ambig))
  
  # Perform matching
  df %<>% mutate(exact = exact_dict[cb_index], HD1 = fuzzy_dict[cb_index], HD1ambig = HD1ambig_dict[cb_index])
  stopifnot((!is.na(df$exact)) + (!is.na(df$HD1)) + df$HD1ambig <= 1)
  rm(dict, exact_dict, fuzzy_dict, HD1ambig_dict) ; invisible(gc())
  
  # Record metadata
  cb_matching_type = c("exact", "HD1", "HD1ambig","none")
  cb_matching_count = c(df$reads[!is.na(df$exact)] %>% sum,
                        df$reads[!is.na(df$HD1)] %>% sum,
                        df$reads[df$HD1ambig] %>% sum,
                        df$reads[is.na(df$exact) & is.na(df$HD1) & !df$HD1ambig] %>% sum)
  stopifnot(sum(df$reads) == sum(cb_matching_count))
  
  # Perform the cb_index conversion
  df1 = df %>% filter(!is.na(exact)) %>% mutate(cb_index = exact) %>% select(1:4)
  df2 = df %>% filter(!is.na(HD1)) %>% mutate(cb_index = HD1) %>% select(1:4)
  df3 = df %>% filter(is.na(exact) & is.na(HD1)) %>% mutate(cb_index = -cb_index) %>% select(1:4)
  df2 %<>% group_by(cb_index, umi_2bit, sb_index) %>% summarize(reads=sum(reads), .groups="drop")
  df12 <- full_join(df1, df2, by = c("cb_index","umi_2bit","sb_index"))
  df12$reads.x %<>% tidyr::replace_na(0) ; df12$reads.y %<>% tidyr::replace_na(0)
  df12 %<>% mutate(reads = reads.x + reads.y) %>% select(-reads.x, -reads.y)
  stopifnot(colnames(df12) == colnames(df3))
  df = rbind(df12, df3)
  stopifnot(df$cb_index != 0)
  stopifnot(sum(df$cb_index > 0) > 0)
  stopifnot(sum(df$reads) == sum(cb_matching_count))
  
  meta <- metadata
  meta$SB_info$remap_10X_CB = remap
  meta$CB_matching = setNames(cb_matching_count, cb_matching_type)
  metadata <<- meta
  
  return(df)
}
df %<>% fuzzy_matching(f("lists/cb_list"), cb_whitelist)
invisible(gc())

# Remove chimeric reads
print("Removing chimeras")
remove_chimeras <- function(df) {
  df %<>% arrange(cb_index, umi_2bit, desc(reads))
  before_same = tidyr::replace_na(df$cb_index==lag(df$cb_index) & df$umi_2bit==lag(df$umi_2bit), FALSE) 
  after_same = tidyr::replace_na(df$cb_index==lead(df$cb_index) & df$umi_2bit==lead(df$umi_2bit) & df$reads==lead(df$reads), FALSE)
  chimeric = before_same | after_same
  meta <- metadata ; meta$SB_filtering %<>% c(reads_chimeric=sum(df[chimeric,]$reads)) ; metadata <<- meta
  return(df[!chimeric,])
}
df %<>% remove_chimeras
invisible(gc())

# Plot raw spatial data
print("Creating barcode rank plots")
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
plot_rankplots(df, f, out_path)
invisible(gc())

### Load puck ##################################################################

# load the puck information
print("Loading the puck")
load_puckdf <- function(f) {
  puckdf = data.frame(sb=f("puck/sb"),
                      x=f("puck/x"),
                      y=f("puck/y"),
                      puck_index=as.integer(f("puck/puck_index")))
  dups = unique(puckdf$sb[duplicated(puckdf$sb)])
  Ns = puckdf$sb[grepl("N", puckdf$sb)]
  
  # Add degeneracy statistics
  most_character_count <- function(vec) {
    stopifnot(typeof(vec) == "character")
    degenA = str_count(vec, "A")
    degenC = str_count(vec, "C")
    degenG = str_count(vec, "G")
    degenT = str_count(vec, "T")
    degenN = str_count(vec, "N")
    return(pmax.int(degenA, degenC, degenG, degenT, degenN))
  }
  puckdf$mc = most_character_count(puckdf$sb)
  longest_run_length <- function(vec) {
    stopifnot(typeof(vec) == "character", nchar(vec) > 0)
    ret = stringr::str_extract_all(vec, "(.)\\1*")
    return(map_int(ret, ~max(nchar(.))))
  }
  puckdf$lr = longest_run_length(puckdf$sb)
  
  puckdfs = map(sort(unique(puckdf$puck_index)), ~filter(puckdf, puck_index==.) %>% select(-puck_index))
  
  # load puck metadata
  meta <- metadata
  meta$puck_info$puck_name = as.character(f("lists/puck_list"))
  meta$puck_info$num_beads = map_int(puckdfs, nrow)
  
  # remove duplicated or low-quality beads
  sb_len = nchar(puckdf$sb[1]) ; mc_tol = round(sb_len*0.75) ; lr_tol = round(sb_len*0.5)
  meta$puck_info$num_dup = map_int(puckdfs, ~sum(.$sb %in% dups))
  meta$puck_info$num_N = map_int(puckdfs, ~sum(.$sb %in% Ns))
  meta$puck_info$num_degen = map_int(puckdfs, ~sum(.$mc > mc_tol | .$lr > lr_tol))
  puckdfs %<>% map(~filter(., !sb %in% dups))
  puckdfs %<>% map(~filter(., !sb %in% Ns))
  puckdfs %<>% map(~filter(., mc <= mc_tol, lr <= lr_tol))
  
  # scale the coordinates (to um)
  get_scaling_factor <- function(bn) {
    if (bn < 150000) {
      k = 0.73
      s = "3mm"
    } else if (bn < 700000) {
      k = 0.73 * 2
      s = "5.5mm"
    } else {
      k = 0.645
      s = "1cm"
    }
    return(c(k, s))
  }
  get_recon_scaling_factor <- function(bn, d) {
    if (bn < 300000) {
      k = 6000 / d
      s = "6mm"
    } else if (bn < 1000000) {
      k = 12000 / d
      s = "1.2cm"
    } else if (bn < 3200000) {
      k = 20000 / d
      s = "2cm"
    } else {
      k = 30000 / d
      s = "3cm"
    }
    return(c(k, s))
  }
  
  res = map2(puckdfs, meta$puck_info$num_beads, function(puckdf, bn) {
    if (mean(c(puckdf$x,puckdf$y)) > 100) {
      # Regular puck
      get_scaling_factor(bn)
    } else {
      # Recon puck
      d = mean(c(max(puckdf$x)-min(puckdf$x), max(puckdf$y)-min(puckdf$y)))
      get_recon_scaling_factor(bn, d)
    }
  })
  meta$puck_info$scaling_factors = map_dbl(res, ~as.numeric(.[[1]]))
  meta$puck_info$puck_sizes = map_chr(res, ~as.character(.[[2]]))
  puckdfs %<>% map2(meta$puck_info$scaling_factors, ~transmute(.x, sb=sb, x=x*.y, y=y*.y))
  
  # center the coordinates
  maxs = map_dbl(puckdfs, ~max(.$x))
  mins = map_dbl(puckdfs, ~min(.$x))
  starts = lag(cumsum(maxs-mins)) %>% tidyr::replace_na(0)
  puckdfs %<>% map2(starts, ~mutate(.x, x = x-min(x)+.y)) # line up from y-axis across
  puckdfs %<>% map(~mutate(., y = y-min(y))) # line up on x-axis
  
  puckdf = do.call(rbind, puckdfs)
  stopifnot(!any(duplicated(puckdf$sb)))
  meta$puck_info$puck_boundaries = c(starts, max(puckdf$x))
  
  # add sb_index
  sb_list = f("lists/sb_list")
  stopifnot(!any(duplicated(sb_list)))
  puckdf$sb_index = match(puckdf$sb, sb_list)
  puckdf %<>% arrange(sb_index) %>% select(sb_index, x, y)
  
  metadata <<- meta
  return(puckdf)
}
puckdf <- load_puckdf(f)

print("Making beadplot")
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
plot_beadplot(df, puckdf, out_path)

### Process results ############################################################

# remove reads with a filtered spatial barcode
print("Removing low-quality spatial barcodes")
m = df$sb_index %in% puckdf$sb_index
metadata$SB_filtering %<>% c(reads_lowQsb=sum(df$reads[!m]))
df = df[m,]

# remove reads that didn't match a called cell
print("Removing non-whitelist cells")
metadata$SB_filtering %<>% c(reads_uncalled=sum(df$reads[df$cb_index < 0]))
df %<>% filter(cb_index > 0)

# count umis
print("Counting UMIs")
metadata$SB_filtering %<>% c(reads_final=sum(df$reads))
df %<>% count_umis
metadata$SB_filtering %<>% c(UMIs_final=sum(df$umi))
invisible(gc())

# join tables (cb_index, sb_index, umi) -> (cb_index, x, y, umi)
stopifnot(df$sb_index %in% puckdf$sb_index)
df %<>% left_join(puckdf, by="sb_index") %>% select(cb_index, x, y, umi) %>% arrange(cb_index, desc(umi))

metadata$puck_info$umi_final = map_int(1:len(metadata$puck_info$puck_name), ~filter(df, x >= metadata$puck_info$puck_boundaries[[.]],
                                                                                        x <= metadata$puck_info$puck_boundaries[[.+1]])$umi %>% sum)

# plot metadata
print("Plotting metadata")
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
plot_metrics(metadata, out_path)

# write output
print("Writing results")
df %>% setNames(c("cb_index","x_um","y_um","umi")) %>% write.table(file.path(out_path, "matrix.csv.gz") %>% gzfile, sep=",", col.names=T, row.names=F, quote=F)
metadata %>% map(as.list) %>% toJSON(pretty=T) %>% writeLines(file.path(out_path, "spatial_metadata.json"))
