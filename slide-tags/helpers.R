library(glue) ; g=glue ; len=length
library(magrittr)
library(purrr)
library(tidyr)
library(dplyr)
library(data.table)
library(stringr)
library(rdist)
library(cowplot)
library(viridisLite)
library(viridis)
library(rhdf5)
library(qpdf)

### Helper methods #############################################################

dec2pct <- function(x, d=2) {return(paste0(round(x*100, digits=d), "%"))}

add.commas <- function(num) {prettyNum(num, big.mark=",")}

h_index <- function(vec) {
  return(data.table(x=vec)[,.(n=.N),x][order(-x), max(pmin(x, cumsum(n)), default=0)])
}

remap_10X_CB <- function(vec) {
  basemap <- setNames(c("AG","TC","CA","GT"), c("TC","AG","GT","CA"))
  return(paste0(substr(vec,1,7), basemap[substr(vec,8,9)], substr(vec,10,16)))
}

determine_remap_10X_CB <- function(vec, dt) {
  if (class(vec) != "character") {
    return(FALSE)
  } else if (!all(nchar(vec) == 16)) {
    return(FALSE)
  } else if (!all(unique(substr(vec,8,9)) %in% c("TC","AG","GT","CA"))) {
    return(FALSE)
  }
  
  # Determine which mapping has more reads
  reads_noremap <- dt[cb %in% cb_whitelist, sum(reads)]
  reads_remap <- dt[cb %in% remap_10X_CB(cb_whitelist), sum(reads)]
  return(reads_remap > reads_noremap)
}

listHD1neighbors <- function(input_string) {
  bases <- c('A','C','G','T','N')
  lapply(1:nchar(input_string), function(i){
    paste0(substr(input_string, 1, i-1),
           setdiff(bases, substr(input_string, i, i)),
           substr(input_string, i+1, nchar(input_string)))
  }) %>% list_c
}

trim_10X_CB <- function(vec) {
  map_chr(vec, ~sub("-[0-9]*$", "", .))
}

trim_puck_name <- function(name) {
  name %<>% str_remove(regex("^recon_", ignore_case=TRUE)) %>% 
            str_remove(regex("^reconstruction_", ignore_case=TRUE)) %>% 
            str_remove(regex("^pucks_", ignore_case=TRUE)) %>% 
            str_remove(regex("^puck_", ignore_case=TRUE)) %>% 
            str_remove(regex(".csv$", ignore_case=TRUE)) %>% 
            str_remove(regex("_Puck$", ignore_case=TRUE)) %>% 
            str_remove(regex("_UMAP[0-9]_.*", ignore_case=TRUE))
  name %<>% map_chr(~paste0(str_split_i(.,"_",-2) %>% str_sub(-floor(30/len(name))),
                            "\n",
                            str_split_i(.,"_",-1) %>% str_sub(-floor(30/len(name)))))
  
  return(name)
}

################################################################################
### Loading methods ############################################################
################################################################################
# DropSift::parseH5ad(matrix_path,
#                     expression_matrix_path = "X",
#                     gene_id_path = "/var/ensembl_ids",
#                     cell_id_path = "/obs/CellID")
# getOptimusGeneSymbols(h5ad_file)

# Supports Optimus and BICAN .h5ad formats
ReadAnnDataX <- function(matrix_path, calledonly=TRUE) {
  fetch <- function(x){return(rhdf5::h5read(matrix_path, x))}
  fields <- rhdf5::h5ls(matrix_path)$name
  
  data <- fetch("X/data") %>% as.numeric
  indices <- fetch("X/indices") %>% as.numeric
  indptr <- fetch("X/indptr") %>% as.numeric
  
  if ("star_IsCell" %in% fields) { # Optimus
    row_names <- fetch("var/_index") %>% as.character %>% make.unique
    col_names <- fetch("obs/_index") %>% as.character %>% make.unique
    mat <- Matrix::sparseMatrix(
      j = indices,
      p = indptr,
      x = data,
      dimnames = list(row_names, col_names),
      index1 = FALSE,
      check = TRUE
    )
    if (calledonly == TRUE) {
      star_IsCell <- fetch("obs/star_IsCell") %>% as.logical
      stopifnot(ncol(mat) == len(star_IsCell))
      mat <- mat[,star_IsCell==TRUE]
    }
  } else if ("CELL_BARCODE" %in% fields) { # BICAN
    mat <- DropSift::parseH5ad(matrix_path,
                               expression_matrix_path = "X",
                               gene_id_path = "/var/gene_symbol", # /var/gene_ids
                               cell_id_path = "/obs/CELL_BARCODE")$dge
    stopifnot(calledonly == TRUE)
  } else {
    stop("Unknown .h5ad format")
  }
  
  return(mat)
}

# Supports 10X, Optimus, and BICAN formats
ReadIntronic <- function(intronic_path, cb_list) {
  fetch <- function(x){return(rhdf5::h5read(intronic_path, x))}
  fields <- rhdf5::h5ls(intronic_path)$name
  cb_list %<>% trim_10X_CB
  
  if (all(c("barcodes", "barcode_idx", "umi_type") %in% fields)) { # 10X
    barcodes <- fetch("barcodes")
    stopifnot(!any(duplicated(barcodes)), cb_list %in% barcodes)
    dt <- data.table(barcode_idx=fetch("barcode_idx")+1,
                     umi_type=fetch("umi_type"))
    pct_intronic <- dt[, .(intronic_umi=sum(umi_type==0), total_umi=.N), barcode_idx][cb_list %>% match(barcodes) %>% match(barcode_idx), intronic_umi/total_umi]
  } else if (all(c("reads_mapped_intronic", "reads_mapped_exonic") %in% fields)) { # Optimus
    dt <- data.table(barcode = fetch("/obs/CellID"),
                     intronic = fetch("/obs/reads_mapped_intronic"),
                     exonic = fetch("/obs/reads_mapped_exonic"))
    stopifnot(!any(duplicated(dt$barcode)), cb_list %in% dt$barcode)
    pct_intronic <- dt[, .(pct_intronic=intronic/(intronic+exonic)), keyby=barcode][cb_list, pct_intronic]
  } else if (all(c("pct_intronic","CELL_BARCODE") %in% fields)) { # BICAN
    pct_intronic <- fetch("/obs/pct_intronic")
    CELL_BARCODE <- fetch("/obs/CELL_BARCODE")
    stopifnot(cb_list %in% CELL_BARCODE)
    pct_intronic <- pct_intronic[match(cb_list, CELL_BARCODE)]
  } else {
    print("Intronic information not found!")
    pct_intronic <- rep(NA, len(cb_list))
  }
  
  return(pct_intronic)
}


ReadSpatialMatrix <- function(f) {
  cb_list <- f("lists/cb_list")
  sb_list <- f("lists/sb_list")
  
  dt <- data.table(cb = f("matrix/cb_index") %>% as.integer %>% factor(levels=seq_along(cb_list), labels=cb_list),
                   umi = f("matrix/umi") %>% as.integer,
                   sb = f("matrix/sb_index") %>% as.integer %>% factor(levels=seq_along(sb_list), labels=sb_list),
                   reads = f("matrix/reads") %>% as.integer)
  
  invisible(gc())
  return(dt)
}


ReadSpatialMetadata <- function(f) {
  metadata = list(SB_info = list(R1s=f("lists/R1_list"),
                                 R2s=f("lists/R2_list"),
                                 pucks=f("lists/puck_list"),
                                 UMI_downsampling=f("metadata/downsampling"),
                                 switch_R1R2=f("metadata/switch") %>% as.logical,
                                 bead_type=f("metadata/bead")),
                  UMI_filtering = setNames(f("metadata/UMI/count"), f("metadata/UMI/type")),
                  UP_matching = setNames(f("metadata/UP/count"), f("metadata/UP/type")),
                  SB_matching = setNames(f("metadata/SB/count"), f("metadata/SB/type")),
                  SB_fuzzy_position = setNames(f("metadata/SB_HD/count"), f("metadata/SB_HD/type")) %>% {.[order(as.integer(names(.)))]},
                  puck_info = list(num_lowQ=f("metadata/num_lowQbeads"))
  )
  metadata$SB_filtering = c(reads_total=f("metadata/reads"),
                            reads_tooshort=f("metadata/R1_tooshort")+f("metadata/R2_tooshort"),
                            reads_noumi=sum(metadata$UMI[c("N","homopolymer")]),
                            reads_noup=sum(metadata$UP[c("none","GG")]),
                            reads_nosb=sum(metadata$SB_matching[c("none","HD1ambig")]))
  return(metadata)
}


ReadPuck <- function(f) {
  # Read matrix from .h5 file
  puckdf <- data.table(sb=f("puck/sb") %>% as.character,
                       x=f("puck/x") %>% as.numeric,
                       y=f("puck/y") %>% as.numeric,
                       puck_index=f("puck/puck_index") %>% as.integer)
  
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
  puckdf[, mc := most_character_count(sb)]
  
  longest_run_length <- function(vec) {
    stopifnot(typeof(vec) == "character", nchar(vec) > 0)
    ret = stringr::str_extract_all(vec, "(.)\\1*")
    return(map_int(ret, ~max(nchar(.))))
  }
  puckdf[, lr := longest_run_length(sb)]
  
  sb_len = max(nchar(puckdf$sb))
  mc_tol = round(sb_len*4/5)
  lr_tol = round(sb_len*2/3)
  
  dups = unique(puckdf$sb[duplicated(puckdf$sb)])
  Ns = puckdf[grepl("N", sb), sb]
  
  # split into pucks
  puckdfs <- split(puckdf, by = "puck_index", sorted = TRUE)
  for (puck in puckdfs) {puck[, puck_index := NULL]}
  
  # save metadata
  puckmeta <- list()
  puckmeta$puck_name = as.character(f("lists/puck_list")) %>% unname
  puckmeta$num_beads = map_int(puckdfs, nrow) %>% unname
  puckmeta$num_dup = map_int(puckdfs, ~.[sb %in% dups, .N]) %>% unname
  puckmeta$num_N = map_int(puckdfs, ~.[sb %in% Ns, .N]) %>% unname
  puckmeta$num_degen = map_int(puckdfs, ~.[mc > mc_tol | lr > lr_tol, .N]) %>% unname
  
  # remove duplicated or low-quality beads
  puckdfs %<>% map(~.[!sb %in% dups])
  puckdfs %<>% map(~.[!sb %in% Ns])
  puckdfs %<>% map(~.[mc <= mc_tol & lr <= lr_tol])
  for (puck in puckdfs) {puck[, mc := NULL]}
  for (puck in puckdfs) {puck[, lr := NULL]}
  
  # for multiplexed pucks, stack along "x" to prevent overlapping
  if (len(puckdfs) > 1) {
    maxs <- map_dbl(puckdfs, ~max(.$x, na.rm=TRUE)) %>% unname
    mins <- map_dbl(puckdfs, ~min(.$x, na.rm=TRUE)) %>% unname
    starts <- cumsum(maxs-mins) %>% lag(1, 0)
    puckdfs %>% map2(starts, ~.x[, x := x-min(x)+.y])
  } else {
    starts <- min(puckdfs[[1]]$x)
  }
  
  puckdf <- rbindlist(puckdfs)
  stopifnot(!any(duplicated(puckdf[, sb])))
  puckmeta$puck_boundaries = c(starts, max(puckdf$x))
  puckdf[, sb := factor(sb)]
  
  return(list(puckdf, puckmeta))
}
