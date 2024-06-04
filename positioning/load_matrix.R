# Load the spatial barcode count matrix and filter to a cell-barcode whitelist
library(glue) ; g=glue ; len=length
library(dplyr)
library(purrr)
library(rlist)
library(rhdf5)
library(ggplot2)
library(cowplot)
library(magrittr)

# sb_path = "SBcounts/SBcounts.h5"
# cb_path = "CBwhitelist.txt"
# out_path <- "."
# setwd("~/spatial")

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
stopifnot(system(g("h5ls {sb_path}/matrix"), intern=T) %>% strsplit(split = "\\s+") %>% map_chr(pluck(1)) == c("cb_index", "reads", "sb_index", "umi"))
f <- function(p){return(h5read(sb_path, p))}
metadata = list()

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

### Workflow ###################################################################

# Fuzzy match and convert cb_index from a cb_list index to a cb_whitelist index
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
  df2 %<>% group_by(cb_index, umi_2bit, sb_index) %>% summarize(reads=sum(reads)) %>% ungroup
  df12 <- full_join(df1, df2, by = c("cb_index","umi_2bit","sb_index"))
  df12$reads.x %<>% tidyr::replace_na(0) ; df12$reads.y %<>% tidyr::replace_na(0)
  df12 %<>% mutate(reads = reads.x + reads.y) %>% select(-reads.x, -reads.y)
  stopifnot(colnames(df12) == colnames(df3))
  df = rbind(df12, df3)
  stopifnot(df$cb_index != 0)
  stopifnot(sum(df$reads) == sum(cb_matching_count))
  
  metadata = list(remap_10X_CB = remap,
                  cb_matching=setNames(cb_matching_count, cb_matching_type)
  )
  res = list(df, metadata)
  return(res)
}
res <- fuzzy_matching(df, f("lists/cb_list"), cb_whitelist)
df <- res[[1]]
metadata %<>% c(res[[2]])
rm(res) ; invisible(gc())

# Remove chimeric reads
print("Removing chimeras")
remove_chimeras <- function(df) {
  df %<>% arrange(cb_index, umi_2bit, desc(reads))
  before_same = tidyr::replace_na(df$cb_index==lag(df$cb_index) & df$umi_2bit==lag(df$umi_2bit), FALSE) 
  after_same = tidyr::replace_na(df$cb_index==lead(df$cb_index) & df$umi_2bit==lead(df$umi_2bit) & df$reads==lead(df$reads), FALSE)
  chimeric = before_same | after_same
  metadata = list("SB_reads_filtered_chimeric" = sum(df[chimeric,]$reads))
  return(list(df[!chimeric,], metadata))
}
res <- remove_chimeras(df)
df <- res[[1]]
metadata %<>% c(res[[2]])
rm(res) ; invisible(gc())

# Plot raw spatial data
print("Creating barcode rank plots")
plot_rankplots <- function(df, f, out_path) {
  gdf = count_umis(df)
  
  # p1
  cb.data = gdf %>% group_by(cb_index) %>% summarize(umi=sum(umi)) %>% arrange(desc(umi)) %>% {mutate(.,index=1:nrow(.), filter="all cell barcodes")}
  cb.data2 = cb.data %>% filter(cb_index>0) %>% {mutate(.,index=1:nrow(.), filter="called cell barcodes only")}
  sb_pct_in_called_cells = round(sum(filter(cb.data,cb_index>0)$umi)/sum(cb.data$umi)*100,2)
  p1 = ggplot(mapping=aes(x=index, y=umi,col=filter))+geom_line(data=cb.data)+geom_line(data=cb.data2) +
    scale_x_log10()+scale_y_log10()+theme_bw()+ggtitle("SB UMI per cell")+ylab("SB UMI counts")+xlab("Cell barcodes") +
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"), legend.title=element_blank()) +
    annotate("text", x = Inf, y = Inf, label = g("SB UMI in called cells: {sb_pct_in_called_cells}%"), hjust = 1.02, vjust = 1.33)
  rm(cb.data, cb.data2) ; invisible(gc())
  
  # p2
  sb.data = gdf %>% group_by(sb_index) %>% summarize(umi=sum(umi)) %>% arrange(desc(umi)) %>% {mutate(.,index=1:nrow(.),filter="all cell barcodes")}
  sb.data2 = gdf %>% filter(cb_index > 0) %>% group_by(sb_index) %>% summarize(umi=sum(umi)) %>% arrange(desc(umi)) %>% {mutate(.,index=1:nrow(.),filter="called cell barcodes only")}
  p2 = ggplot(mapping=aes(x=index,y=umi,col=filter))+geom_line(data=sb.data)+geom_line(data=sb.data2)+
    scale_x_log10()+scale_y_log10()+theme_bw()+ggtitle("SB UMI per bead")+ylab("SB UMI counts")+xlab("Beads")+
    theme(legend.position = c(0.05, 0.05), legend.justification = c("left", "bottom"), legend.background = element_blank(), legend.spacing.y = unit(0.1,"lines"), legend.title=element_blank())
  rm(sb.data, sb.data2) ; invisible(gc())
  
  # p3
  x = seq(0, 1, 0.05) * f("metadata/num_reads")/1000000
  plot.df = data.frame(x=x, y=f("metadata/downsampling")/1000000)
  p3 = ggplot(plot.df, aes(x=x,y=y)) + geom_point() + theme_bw() + 
    xlab("Millions of reads") + ylab("Millions of filtered SB UMIs") + ggtitle("SB downsampling curve")
  
  # p4
  d = table(df$reads) %>% as.data.frame %>% setNames(c("reads","count")) %>% mutate(reads = as.numeric(levels(reads))[reads])
  d %<>% filter(reads < 10) %>% bind_rows(data.frame(reads = 10, count = sum(d$count[d$reads >= 10])))
  sequencing_saturation = round((1 - nrow(df) / sum(df$reads)) * 100, 2) %>% paste0("%")
  p4 = ggplot(d, aes(x=reads, y=count/1000/1000)) + geom_col() +
    theme_bw() + xlab("Reads per UMI") + ylab("Millions of filtered SB UMIs") + ggtitle("SB read depth") + 
    annotate("text", x = Inf, y = Inf, label = g("sequencing saturation = {sequencing_saturation}"), hjust = 1.02, vjust = 1.33) +
    scale_x_continuous(breaks=min(d$reads):max(d$reads), labels=(min(d$reads):max(d$reads)) %>% {ifelse(.==10, "10+", .)})

  plot = plot_grid(p1, p2, p3, p4, ncol=2)
  
  make.pdf <- function(plots, name, w, h) {
    if ("gg" %in% class(plots) || class(plots)=="Heatmap") {plots = list(plots)}
    pdf(file=name, width=w ,height=h)
    lapply(plots, function(x){print(x)})
    dev.off()
  }
  make.pdf(plot, file.path(out_path,"spatial_rankplots.pdf"), 7, 8)
  
  metadata = list(sb_pct_in_called_cells = sb_pct_in_called_cells, sequencing_saturation = sequencing_saturation)
  return(metadata)
}
res = plot_rankplots(df, f, out_path)
metadata %<>% c(res)
rm(res) ; invisible(gc())

# remove reads that didn't match a called cell
print("Removing non-whitelist cells")
df %<>% filter(cb_index > 0)

# count umis
print("Counting UMIs")
df %<>% count_umis

# Compute metrics
# Misc(obj, "SB_reads_final") <- sum(df$reads)
# Misc(obj, "SB_umi_final") <- df %>% count_umis %>% pull(umi) %>% sum
# qsave(df, "df.qs")
# invisible(gc())

print("Writing results")
write.table(df, "matrix.csv", sep=",", col.names=T, row.names=F, quote=F)
