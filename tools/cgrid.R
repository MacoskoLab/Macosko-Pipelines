library(glue) ; g=glue ; len=length
library(magrittr)
library(ggplot2)
library(stringr)
library(Seurat)
library(shiny)
library(rlist)
library(rdist)
library(purrr)
library(sf)
library(qs)

### Helper methods #############################################################

plot.new() # needed to initialize make.spline()

unfactor <- function(vec) {
  if(!is.factor(vec)) {return(vec)}
  vec = levels(vec)[vec]
  if (!any(is.na(suppressWarnings(as.numeric(vec))))) {
    vec %<>% as.numeric
  }
  return(vec)
}

plot.tab <- function(df) {return(cowplot::plot_grid(gridExtra::tableGrob(df)))}

scalevec <- function(x){stopifnot(!is.na(x)) ; vec=(x-min(x))/(max(x)-min(x)) ; stopifnot(max(vec)==1&min(vec)==0) ; return(vec)}

# Given a data vector, return a color vector
vec2col <- function(vec, cont=F) {
  if (len(unique(vec)) == 1) {
    col = rep("#000000",len(vec))
  } else if (!cont) {
    clusters = purrr::discard(gtools::mixedsort(unique(vec)), is.na)
    pal = scales::hue_pal()(len(clusters))
    col = pal[match(vec, clusters)]
  } else {
    stopifnot(is.na(vec) | !is.na(suppressWarnings(as.numeric(vec))))
    vec = as.numeric(vec)
    vec = vec - min(vec, na.rm=T)
    vec = vec / max(vec, na.rm=T)
    col = map_chr(vec, function(i){ if (is.na(i)) {return(NA)} else {return(viridis::magma(n=1,begin=i,end=i))} })
  }
  col[is.na(vec)] <- "grey"
  stopifnot(!is.na(col))
  return(col)
}

# assert that the first two columns of a df (x,y) are numeric and finite
validate.xydf <- function(df) {
  stopifnot(class(df) == "data.frame", ncol(df) >= 2)
  stopifnot(names(df)[1:2] == c("x","y"))
  stopifnot(is.numeric(df[[1]]), is.numeric(df[[2]]))
  stopifnot(is.finite(df[[1]]), is.finite(df[[2]])) # catch NA, NaN, Inf
}

# ensure that the cgrid (x,y,i,r) is coherent
validate.cgrid <- function(cgrid) {
  stopifnot(class(cgrid) == "data.frame", ncol(cgrid) >= 4)
  stopifnot(names(cgrid)[1:4] == c("x","y","i","r"))
  stopifnot(range(cgrid$i) == c(0,1), range(cgrid$r) == c(0,1))
  stopifnot(len(unique(table(cgrid$i))) == 1, len(unique(table(cgrid$r))) == 1)
  stopifnot(unique(table(cgrid$i)) * unique(table(cgrid$r)) == nrow(cgrid))
  map(cgrid, ~stopifnot(is.numeric(.)))
  map(cgrid, ~stopifnot(is.finite(.)))
}

# get the area enclosed by an (x,y) df
polygon.area <- function(polygon_df) {
  validate.xydf(polygon_df)
  polygon_df %<>% dplyr::select(1,2) %>% setNames(c("x","y"))
  poly <- rbind(polygon_df, polygon_df[1,]) %>% as.matrix %>% list %>% st_polygon
  return(st_area(poly))
}

# points_df and polygon_df are 2-column dataframes of (x,y), return inclusion mask
in.bounds <- function(points_df, polygon_df) { 
  validate.xydf(points_df) ; validate.xydf(polygon_df)
  points_df %<>% dplyr::select(1,2) ; points_df %<>% setNames(c("x","y"))
  polygon_df %<>% dplyr::select(1,2) ; polygon_df %<>% setNames(c("x","y"))
  
  points_df %<>% dplyr::mutate(i = dplyr::row_number())
  m <- is.na(points_df[[1]]) | is.na(points_df[[2]])
  xy <- points_df %>% dplyr::filter(!is.na(x) & !is.na(y)) %>% st_as_sf(coords=c("x","y"))
  poly <- rbind(polygon_df, polygon_df[1,]) %>% as.matrix %>% list %>% st_polygon
  
  inds <- points_df$i %in% xy$i[unlist(st_contains(poly, xy))]
  inds[m] <- NA
  return(inds)
}

# input (x,y) guide points df, return (x,y) spline df
make.spline <- function(df, shape=-0.75) {
  stopifnot(class(df) == "data.frame", ncol(df) == 2)
  stopifnot(is.numeric(df[[1]]), is.numeric(df[[2]]))
  stopifnot(is.finite(df[[1]]), is.finite(df[[2]])) # catch NA, NaN, Inf
  if (nrow(df) >= 2) {
    spl = graphics::xspline(x=df[[1]], y=df[[2]], shape=shape, open=T, repEnds=T, draw=F)
    spl <- as.data.frame(spl)
  } else {
    spl = df %>% dplyr::select(1,2) %>% setNames(c("x","y"))
  }
  stopifnot(ncol(spl)==2)
  validate.xydf(spl)
  return(spl)
}

# add splines to an existing ggplot, each param is a df or list of dfs
plot.splines <- function(plot, ...) {
  dfs <- list(...) ; dfs %<>% map(~if(class(.)=="data.frame") {list(.)} else {.}) # turn into list of lists
  s = purrr::map_depth(dfs, .depth=2, ~geom_path(data=., aes(x,y), linewidth=.3))
  s = purrr::flatten(s)
  return(Reduce(`+`, s, init=plot))
}

# add points to an existing ggplot, each param is a df or list of dfs
plot.points <- function(plot, ...) {
  dfs <- list(...) ; dfs %<>% map(~if(class(.)=="data.frame") {list(.)} else {.}) # turn into list of lists
  s = purrr::map_depth(dfs, .depth=2, ~geom_point(data=., aes(x,y)))
  s = purrr::flatten(s)
  return(Reduce(`+`, s, init=plot))
}

# input x,y df, get list of segment lengths
path.dists <- function(df) {
  validate.xydf(df)
  dists = sqrt((df[,1]-dplyr::lag(df[,1]))**2+(df[,2]-dplyr::lag(df[,2]))**2) %>% tidyr::replace_na(0)
  return(dists)
}

# input x,y path df and distances along the path, get x,y df with coord at each distance
dist2xy <- function(df, ds) {
  validate.xydf(df) ; stopifnot(len(ds) > 0, is.numeric(ds), is.finite(ds))
  dists = cumsum(path.dists(df))
  res <- purrr::map(ds, function(d) {
    if (d > max(dists)) {print(g("setting {d} to {max(dists)}")) ; d = max(dists)}
    if (d < 0) {print(g("setting {d} to {0}")) ; d = 0}
    if (d %in% dists) { # return the point if it matches exactly
      return(unlist(df[which(d==dists)[[1]],c(1,2)]))
    } else { # return a weighted average if not
      i = sum(dists < d)
      p = (d-dists[[i]])/(dists[[i+1]]-dists[[i]])
      avg = c(x=df[i,1]*(1-p)+df[i+1,1]*(p), y=df[i,2]*(1-p)+df[i+1,2]*(p))
      return(avg) 
    }
  })
  res <- do.call(rbind, res) %>% as.data.frame %>% setNames(c("x","y"))
  validate.xydf(res) ; stopifnot(ncol(res) == 2)
  return(res)
}

# input x,y df, return evenly-spaced x,y df of n points
equalize <- function(df, n) {
  validate.xydf(df)
  dists <- seq(0, sum(path.dists(df)), length.out=n)
  return(dist2xy(df, dists))
}

# take in cgrid (x,y,i,r), return list of splines (x,y)
cgrid2splines <- function(cgrid) {
  dfs1 = purrr::map(unique(cgrid$i), ~dplyr::filter(cgrid, i==.) %>% dplyr::arrange(r) %>% dplyr::select("x","y"))
  dfs2 = purrr::map(unique(cgrid$r), ~dplyr::filter(cgrid, r==.) %>% dplyr::arrange(i) %>% dplyr::select("x","y"))
  return(append(dfs1,dfs2))
}

# get the median width of a cgrid (x,y,i,r)
cgrid.width <- function(cgrid) {
  validate.cgrid(cgrid)
  bot = cgrid %>% dplyr::filter(i==min(i)) %>% dplyr::arrange(r)
  top = cgrid %>% dplyr::filter(i==max(i)) %>% dplyr::arrange(r)
  stopifnot(bot$r == top$r)
  avg <- mean(sqrt((bot$x-top$x)**2 + (bot$y-top$y)**2))
  return(avg)
}

# take in list of splines, num layers, num columns, return cgrid (x,y,i,r)
makegrid1 <- function(splines, num.i, num.r) {
  stopifnot(class(splines) == "list", len(splines) >= 2) ; names(splines) <- NULL
  map(splines, function(spline){validate.xydf(spline) ; stopifnot(ncol(spline) == 2)})
  stopifnot(is.numeric(num.i), is.numeric(num.r), is.finite(num.i), is.finite(num.r), num.i > 1, num.r > 1)
  splines %<>% purrr::map(~equalize(.,num.r))
  dfs <- purrr::map(splines, ~split(., seq(nrow(.)))) %>% purrr::transpose() %>% purrr::map(~do.call("rbind", .)) # convert n m-point dfs into m n-point dfs
  dfs %<>% map(~equalize(.,num.i))
  dfs %<>% purrr::imap(~dplyr::mutate(.x, i=dplyr::row_number(), r=as.numeric(.y)))
  cgrid = do.call(rbind, dfs) ; cgrid$i %<>% scalevec ; cgrid$r %<>% scalevec
  validate.cgrid(cgrid)
  return(cgrid)
}

# take in two splines, num layers, num columns, return cgrid (x,y,i,r)
makegrid2 <- function(spl1, spl2, num.i, num.r) {
  validate.xydf(spl1) ; validate.xydf(spl2)
  stopifnot(ncol(spl1) == 2, ncol(spl2) == 2)
  spl1 %<>% dplyr::select(1,2) %>% setNames(c("x","y"))
  spl2 %<>% dplyr::select(1,2) %>% setNames(c("x","y"))
  stopifnot(is.numeric(num.i), is.numeric(num.r), is.finite(num.i), is.finite(num.r), num.i > 1, num.r > 1)
  
  d.end1 = sum(path.dists(spl1)) ; d.end2 = sum(path.dists(spl2))
  d = (d.end1+d.end2)/(num.r-1) # spline step distance
  df = data.frame(a=0, b=0) # distances along each spline to mark points
  
  a = 0 ; b = 0 # these must start at 0 - represents the distance moved along the spline
  sr = 100 ; z = d/5 # sr is the search resolution, z is is min distance moved
  for (it in 1:(num.r-1)) {
    if (abs(d.end1+d.end2-a-b) < 1.5*d) { a = d.end1 ; b = d.end2 } # only one step left - set to end
    else if (a == d.end1) { b %<>% add(d) } 
    else if (b == d.end2) { a %<>% add(d) }
    else if (a+z >= d.end1) { b %<>% add(d-(d.end1-a)) ; a = d.end1 }
    else if (b+z >= d.end2) { a %<>% add(d-(d.end2-b)) ; b = d.end2 }
    else {
      d1 = a+seq(z, d-z, length.out=sr) ; d2 = b+seq(d-z, z, length.out=sr)
      m = (d1 <= d.end1) & (d2 <= d.end2) ; stopifnot(sum(m)>0)
      d1 = d1[m] ; d2 = d2[m] # filter moves that go beyond the end
      p1 = dist2xy(spl1, d1) ; p2 = dist2xy(spl2, d2)
      i = which.min((p1[[1]]-p2[[1]])**2 + (p1[[2]]-p2[[2]])**2)
      a = d1[[i]] ; b = d2[[i]]
    }
    df %<>% rbind(c(a,b))
  }
  
  pts1 <- dist2xy(spl1, df$a) ; pts2 <- dist2xy(spl2, df$b) ; stopifnot(nrow(pts1)==nrow(pts2))
  dfs <- map(1:nrow(pts1), ~rbind(pts1[.,], pts2[.,])) # transpose
  dfs %<>% map(~equalize(.,num.i))
  dfs %<>% purrr::imap(~dplyr::mutate(.x, i=dplyr::row_number(), r=as.numeric(.y)))
  cgrid = do.call(rbind, dfs) ; cgrid$i %<>% scalevec ; cgrid$r %<>% scalevec
  validate.cgrid(cgrid)
  return(cgrid) 
}

# x,y,r of shorter is set to match x,y,r of longer (default)
# i is concatenated and rescaled, weighted by cgrid thickness
cgridmerge <- function(cgrid1, cgrid2, flip=F) {
  fetch <- function(cgrid, i, r) {return(unlist(dplyr::filter(cgrid, i==!!i, r==!!r)[c("x","y")]))}
  cgrid2above1 <- function(cgrid1, cgrid2) {
      all(fetch(cgrid1,1,0)==fetch(cgrid2,0,0) & fetch(cgrid1,1,1)==fetch(cgrid2,0,1))
  }
  
  # Validate both cgrids, and assert that they are adjacent
  validate.cgrid(cgrid1) ; validate.cgrid(cgrid2)
  if(!cgrid2above1(cgrid1, cgrid2)) {tmp=cgrid1; cgrid1=cgrid2; cgrid2=tmp}
  stopifnot(cgrid2above1(cgrid1, cgrid2))
  
  # Get average cgrid widths
  avg1 <- cgrid.width(cgrid1)
  avg2 <- cgrid.width(cgrid2)
  
  # Adjust r
  update_cgrid <- function(cgrid1, cgrid2) { # first argument gets r changed to match the second
      # Get the seams
      if (cgrid2above1(cgrid1, cgrid2)) {
          path1 = cgrid1 %>% dplyr::filter(i==max(i)) %>% dplyr::arrange(r)
          path2 = cgrid2 %>% dplyr::filter(i==min(i)) %>% dplyr::arrange(r)
      } else if (cgrid2above1(cgrid2, cgrid1)) {
          path1 = cgrid1 %>% dplyr::filter(i==min(i)) %>% dplyr::arrange(r)
          path2 = cgrid2 %>% dplyr::filter(i==max(i)) %>% dplyr::arrange(r)
      }
      pdists = cdist(path1[,c("x","y")], path2[,c("x","y")])
      
      # compute p and (i,j) for each cgrid2 seam point
      ijpr = purrr::map_dfr(1:ncol(pdists), function(i) {
          twosmallest <- data.frame(i=1:nrow(pdists), d=pdists[,i]) %>% dplyr::arrange(d)
          p = twosmallest$d[[1]] / (twosmallest$d[[1]]+twosmallest$d[[2]]) ; if (!is.finite(p)) { p = 0 }
          return(c(i=twosmallest$i[[1]], j=twosmallest$i[[2]], p=p))
      }) %>% dplyr::mutate(r = path2$r)
      stopifnot(!duplicated(ijpr$r))
      
      cgrid1 = pmap_dfr(ijpr, function(i, j, p, r) {
          col1 <- dplyr::filter(cgrid1, r==path1$r[[!!i]])
          col2 <- dplyr::filter(cgrid1, r==path1$r[[!!j]])
          stopifnot(col1$i == col2$i)
          res <- col1[,c("x","y")]*(1-p) + col2[,c("x","y")]*(p)
          res$i = col1$i ; res$r = r
          return(res)
      })
      validate.cgrid(cgrid1)
      return(cgrid1)
  }
  d1 = sum(path.dists(cgrid1 %>% dplyr::filter(i==min(i)) %>% dplyr::arrange(r)))
  d2 = sum(path.dists(cgrid2 %>% dplyr::filter(i==max(i)) %>% dplyr::arrange(r)))
  if (xor(d2 > d1, flip)) { # change cgrid1
      cgrid1 = update_cgrid(cgrid1, cgrid2)
  } else { # change cgrid2
      cgrid2 = update_cgrid(cgrid2, cgrid1)
  }
  stopifnot(sort(unique(cgrid1$r)) == sort(unique(cgrid2$r)))
  stopifnot(range(cgrid1$i)==c(0,1), range(cgrid2$i)==c(0,1))
  
  # Adjust i
  cgrid2 %<>% dplyr::filter(i > 0) # arbitrary, could pick cgrid1 i < 1
  cgrid1$i = cgrid1$i * avg1/(avg1+avg2)
  cgrid2$i = cgrid2$i * avg2/(avg1+avg2) + avg1/(avg1+avg2)
  
  # rbind the cgrids together
  cgrid <- rbind(cgrid1, cgrid2) %>% dplyr::mutate(i=scalevec(i))
  cgrid %<>% dplyr::arrange(r, i)
  validate.cgrid(cgrid)
  return(cgrid)
}



# assign i,r to x,y using a cgrid
xy2ir <- function(xy, cgrid) {
  xy %<>% dplyr::select(1:2) ; stopifnot(names(xy)==c("x","y"))
  cgrid %<>% dplyr::select(1:4) ; validate.cgrid(cgrid)
  xy %<>% dplyr::mutate(row=dplyr::row_number())
  
  # split in and out points
  dfs = list(cgrid %>% dplyr::filter(i==1) %>% dplyr::arrange(r),
             cgrid %>% dplyr::filter(r==1) %>% dplyr::arrange(dplyr::desc(i)),
             cgrid %>% dplyr::filter(i==0) %>% dplyr::arrange(dplyr::desc(r)),
             cgrid %>% dplyr::filter(r==0) %>% dplyr::arrange(i)
  )
  stopifnot(dfs[[1]] %>% dplyr::slice_tail(n=1) == dfs[[2]] %>% dplyr::slice_head(n=1))
  stopifnot(dfs[[2]] %>% dplyr::slice_tail(n=1) == dfs[[3]] %>% dplyr::slice_head(n=1))
  stopifnot(dfs[[3]] %>% dplyr::slice_tail(n=1) == dfs[[4]] %>% dplyr::slice_head(n=1))
  stopifnot(dfs[[4]] %>% dplyr::slice_tail(n=1) == dfs[[1]] %>% dplyr::slice_head(n=1))
  poly <- do.call(rbind,dfs)
  m <- in.bounds(xy,poly[,c("x","y")])
  indf <- xy[m,] ; outdf <- xy[!m,]
  
  # assign <out> points
  avg <- cgrid.width(cgrid[,1:4])
  index = cdist(outdf[,c("x","y")],poly[,c("x","y")]) %>% apply(1,which.min)
  outdf$d <- sqrt((outdf$x - poly[index,]$x)^2 + (outdf$y - poly[index,]$y)^2)
  outdf$r <- poly[index,]$r
  outdf$sign <- poly[index,]$i * 2 - 1
  outdf$i = poly[index,]$i + outdf$sign * outdf$d / avg
  outdf %<>% dplyr::select(row,i,r)
  
  point_line_distance <- function(pt, l1, l2) {
    return(abs((l2[[1]]-l1[[1]])*(l1[[2]]-pt[[2]])-(l1[[1]]-pt[[1]])*(l2[[2]]-l1[[2]]))/sqrt((l2[[1]]-l1[[1]])^2 + (l2[[2]]-l1[[2]])^2))
  }
  
  # assign <in> points
  cgrid %<>% dplyr::mutate(ii=match(i,sort(unique(i))),ri=match(r,sort(unique(r))))
  row_polys <- map(2:max(cgrid$ii), ~rbind(cgrid[cgrid$ii==.-1,] %>% dplyr::arrange(r) %>% dplyr::select(c("x","y")),
                                          cgrid[cgrid$ii==.,]%>% dplyr::arrange(dplyr::desc(r)) %>% dplyr::select(c("x","y"))))
  col_polys <- map(2:max(cgrid$ri), ~rbind(cgrid[cgrid$ri==.-1,] %>% dplyr::arrange(i) %>% dplyr::select(c("x","y")),
                                           cgrid[cgrid$ri==.,]%>% dplyr::arrange(dplyr::desc(i)) %>% dplyr::select(c("x","y"))))
  row_i <- do.call(cbind,map(row_polys,~in.bounds(indf,.))) %>% apply(1, which) ; stopifnot(class(row_i)=="integer") # point is in exactly one polygon
  col_i <- do.call(cbind,map(col_polys,~in.bounds(indf,.))) %>% apply(1, which) ; stopifnot(class(col_i)=="integer") # point is in exactly one polygon
  res <- pmap(list(row_i,col_i),~cgrid[(cgrid$ii==.x|cgrid$ii==.x+1)&(cgrid$ri==.y|cgrid$ri==.y+1),c("x","y","i","r")] %>% dplyr::arrange(r,i))
  stopifnot(len(res)==nrow(indf), names(indf)[1:2]==c("x","y")) ; invisible(map(res,~stopifnot(nrow(.)==4,names(.)[1:2]==c("x","y"))))
  rs <- map(1:nrow(indf), ~c(res[[.]][1,"r"], point_line_distance(indf[.,], res[[.]][1,], res[[.]][2,]),
                             res[[.]][3,"r"], point_line_distance(indf[.,], res[[.]][3,], res[[.]][4,])
                             )) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("r1","d1","r2","d2"))
  is <- map(1:nrow(indf), ~c(res[[.]][1,"i"], point_line_distance(indf[.,], res[[.]][1,], res[[.]][3,]),
                             res[[.]][2,"i"], point_line_distance(indf[.,], res[[.]][2,], res[[.]][4,])
                             )) %>% {do.call(rbind,.)} %>% as.data.frame %>% setNames(c("i1","d1","i2","d2"))
  indf$r <- (rs$r1*rs$d2+rs$r2*rs$d1)/(rs$d1+rs$d2)
  indf$i <- (is$i1*is$d2+is$i2*is$d1)/(is$d1+is$d2)
  indf %<>% dplyr::select(row,i,r)
  
  # Return
  res <- rbind(indf,outdf) %>% dplyr::arrange(row) %>% dplyr::select(i,r) ; row.names(res) <- NULL
  return(res)
}

# given an 
obj2ir <- function(obj, cgrid) {
  df <- obj@meta.data %>% dplyr::mutate(row=dplyr::row_number())
  if (all(c("x","y") %in% names(df))) {xy = df[,c("x","y","row")]
  } else if (all(c("x_um","y_um") %in% names(df))) {xy = df[,c("x_um","y_um","row")] %>% setNames(c("x","y","row"))
  } else {stop("Error: no coordinates found in object")}
  xy %<>% dplyr::filter(!is.na(x),!is.na(y))
  xy %<>% cbind(xy2ir(xy, cgrid))
  m <- match(df$row, xy$row)
  df$i <- xy$i[m] ; df$r <- xy$r[m]
  obj %<>% AddMetaData(df$i, col.name="i")
  obj %<>% AddMetaData(df$r, col.name="r")
  emb = obj@meta.data[,c("i","r")]
  colnames(emb) = c("f_1","f_2")
  obj[["flat"]] <- CreateDimReducObject(embeddings = as.matrix(emb), key = "f_")
  return(obj)
}

### Shiny Servers ##############################################################

# Return a seurat object with only the selected clusters
cluster_selector <- function(seurat_obj, col="seurat_clusters") {
  # Pre-process
  stopifnot(col%in%names(seurat_obj@meta.data))
  obj = AddMetaData(seurat_obj, seurat_obj@meta.data[[col]], "col") ; obj$col %<>% unfactor
  clusters = gtools::mixedsort(unique(obj$col)) ; stopifnot(len(clusters) <= 100)
  stopifnot(typeof(clusters) %in% c("integer","double","character"))
  if (is.character(obj$col) || !all(obj$col %% 1 == 0)) {
    wholenums=F
    obj$col2 = g("[{match(obj$col,clusters)}] {obj$col}")
  } else {
    wholenums=T
    obj$col2 = obj$col
  }
  obj$col2 %<>% factor(gtools::mixedsort(unique(obj$col2)))
  
  # Plot and select
  plot = DimPlot(obj, reduction= "spatial", split.by="col2", ncol=ceiling(sqrt(len(clusters)))) & theme_void() & coord_fixed(ratio=1) & NoLegend()
  print(plot) ; print(g("Input clusters: "))
  input = readline() %>% str_split_1(" ") %>% keep(nchar(.)>=1)
  stopifnot(len(input) > 0) ; stopifnot(!any(is.na(suppressWarnings(as.numeric(input)))))
  input %<>% as.numeric ; stopifnot(input %% 1 == 0)
  if (wholenums) {
    m = obj$col %in% input
  } else {
    m = obj$col %in% clusters[input]
  }
  return(seurat_obj[,m])
}

# Return a (x,y,col) dataframe with only the selected cols
layer_selector <- function(df) {
  # Pre-process
  stopifnot(ncol(df) == 3) ; names(df) = c("x","y","col")
  validate.xydf(df) ; df$col %<>% unfactor
  clusters <- gtools::mixedsort(unique(df$col)) ; stopifnot(len(clusters) <= 100)
  ret.df <- NULL
  
  # Shiny server
  ui <- fluidPage(
    titlePanel("Layer Selector"),
    sidebarLayout(
      sidebarPanel(uiOutput("checkboxes"),actionButton("select","Select")),
      mainPanel(plotOutput("plot"))
    )
  )
  server <- function(input, output) {
    output$checkboxes <- renderUI({checkboxGroupInput("selected_values", "Select Clusters:", choices=clusters, selected=clusters)})
    output$plot <- renderPlot({
      df$active[df$col %in% input$selected_values] <- T
      df$active[!df$col %in% input$selected_values] <- F
      plot <- ggplot(dplyr::filter(df,active), aes(x=x,y=y,col=as.factor(col)))+geom_point()+
        coord_fixed(ratio=1)+theme_void()+lims(x=range(df$x),y=range(df$y))+labs(color="")+
        theme(legend.text=element_text(size=11),legend.title=element_text(size=12))+
        guides(color = guide_legend(override.aes=list(size=5)))
      print(plot)
      ret.df <<- dplyr::select(dplyr::filter(df,active), x,y,col)
    })
    observeEvent(input$select, {stopApp(0)})
  }
  print(shinyApp(ui=ui, server=server))
  validate.xydf(ret.df)
  return(ret.df)
}

# Returns an x,y dataframe of clicked points
anno_points <- function(df, lines=F, splines=F, previous_lines=list(), previous_splines=list()) {
  # Pre-process
  if(ncol(df)==2){df%<>%cbind(1)} # add cluster if missing
  stopifnot(ncol(df)==3) ; df %<>% setNames(c("x","y","col"))
  validate.xydf(df)
  df = df[sample(1:nrow(df), size=nrow(df), replace=F),] # shuffle the rows
  df$colors <- vec2col(df$col) # assign a color to each cluster
  ret.df <- NULL ; x <- df$x ; y <- df$y ; colors = df$colors
  
  # Shiny server
  ui <- fluidPage(
    plotOutput("plot", click = "plot_click", height = "600px"),
    actionButton("quit", "Quit"),
    actionButton("restart", "Restart"),
    actionButton("undo", "Undo"),
    actionButton("select", "Select")
  )
  server <- function(input, output) {
    df <- reactiveVal(data.frame(matrix(ncol=2, nrow=0)))
    output$plot <- renderPlot({
      plot(x=x, y=y, col=colors, pch=19, asp=1, cex=1)
      d <- df() ; points(d)
      if (lines) {lines(d,lwd=2)} ; if (splines && nrow(d)>=2) {lines(make.spline(d),lwd=2)}
      for (prev in previous_lines) {lines(prev, lwd=2, col="red")}
      for (prev in previous_splines) {lines(make.spline(prev), lwd=2, col="red")}
    })
    observeEvent(input$quit, {stopApp(0)})
    observeEvent(input$restart, {df(data.frame(matrix(ncol=2, nrow=0)))})
    observeEvent(input$undo, {d <- df() ; d <- d[-nrow(d),] ; df(d)})
    observeEvent(input$select, {ret.df <<- setNames(df(), c("x","y")) ; stopApp(0)})
    observeEvent(input$plot_click, {d <- rbind(df(), c(input$plot_click$x,input$plot_click$y)) ; df(d)})
  }
  print(shinyApp(ui, server))
  stopifnot(!is.na(ret.df$x), !is.na(ret.df$y), is.finite(ret.df$x), is.finite(ret.df$y))
  return(ret.df)
}

# Returns a mask of points inside the polygon
select_points <- function(df) {
  poly <- anno_points(df, lines=T, splines=F)
  print(paste0("Area: ",st_area(rbind(poly, poly[1,]) %>% as.matrix %>% list %>% st_polygon)))
  m <- in.bounds(df, poly)
  print(plot(df[[1]], df[[2]], col=c("blue","red")[m+1], pch=20, asp=1))
  print(g("{sum(m)}/{len(m)} selected"))
  return(m)
}

# Returns spline x,y dataframe
anno_spline <- function(df) {
  spl <- anno_points(df, lines=F, splines=T)
  if (is.null(spl)) {return(data.frame(matrix(ncol=2, nrow=0)) %>% setNames(c("x","y")))}
  if (nrow(spl)<2) {return(spl)}
  return(make.spline(spl))
}

# Returns list of spline x,y dataframes
anno_splines <- function(df) {
  dfs <- list()
  repeat {
    spl <- anno_points(df, lines=F, splines=T, previous_splines=dfs)
    if (is.null(spl) || nrow(spl)<2) {break}
    dfs %<>% list.append(make.spline(spl))
  }
  return(dfs)
}

### Workflow ###################################################################

obj = qread("~/analyses/9positioning/seurat15.qs")
obj %<>% cluster_selector
df <- data.frame(x=obj$x_um,y=obj$y_um,col=obj$seurat_clusters) %>% dplyr::filter(is.finite(x),is.finite(y))
df %<>% layer_selector
splines <- anno_splines(df)
DimPlot(obj,reduction="spatial") %>% plot.splines(splines)

# Spline drawing tips:
# - all splines must be drawn in order and with the same orientation
# - all interior points must be included between the spline edges

# Method 1 (simple)
cgrid <- makegrid1(splines, 10, 100)
ggplot(cgrid,aes(x=x,y=y,col=i)) + geom_point()
ggplot(cgrid,aes(x=x,y=y,col=r)) + geom_point()
DimPlot(obj,reduction="spatial") %>% plot.splines(cgrid2splines(cgrid))

# Method 2 (advanced)
cgrids <- map(2:len(splines), ~makegrid2(splines[[.-1]], splines[[.]], 20, 200))
cgrid <- Reduce(cgridmerge,cgrids)
#ggplot(cgrid, aes(x=x,y=y,col=i)) + geom_point()
#ggplot(cgrid, aes(x=x,y=y,col=r)) + geom_point()
DimPlot(obj,reduction="spatial") %>% plot.splines(cgrid2splines(cgrid))

# Add coords to object
obj %<>% obj2ir(cgrid)
DimPlot(obj,reduction="flat")



# todo: outpoint r

# First spline is always at 0
# Last spline is always at 1
a <- map(sort(unique(cgrid$r)),function(rr){cgrid %>% dplyr::filter(abs(r-rr)<1E-10) %>% dplyr::arrange(i) %>% path.dists %>% {abs(dplyr::lead(.)-.)>1e-10}})
a <- map(a,function(vec){vec %>% head(-1) %>% tail(-1) %>% which %>% add(1)})
sort(unique(cgrid$i))[Reduce(intersect,a)]

# store the boundaries in the object

# unused
# select_points
# 0.75 and flip
