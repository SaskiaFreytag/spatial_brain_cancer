library(vegan)
library(Seurat)
library(tidyverse)
library(parallel)
library(pals)
library(parallelDist)

cell_to_tile <- function(df_final, sce,  tile_shift, cores=10, prob=FALSE){
  
  index <- split(1:nrow(df_final), df_final$fov)
  df_final_split <- lapply(index, function(x) df_final[x,]) 
  
  index <- split(1:ncol(sce), sce$fov)
  sce_split <- lapply(index, function(x) as.data.frame(colData(sce)[x,
                                                                    c("CenterX_local_px", "CenterY_local_px", "Barcode")]))
  df_final_split <- df_final_split[names(sce_split)]
  
  all_cluster <- mapply(Y = df_final_split, X = sce_split, 
          function(X,Y) identify_cell_in_tile(X, Y, tile_shift, cores=10, prob=prob))
  
  
  return(all_cluster)
  
}

find_tile <- function(x, df_final_1, tile_shift) {
   
  which(df_final_1[,c("x")]<x[1] & x[1]<df_final_1[,c("x")]+tile_shift &
  x[2]> df_final_1[,c("y")] & x[2]<df_final_1[,c("y")]+tile_shift)
  
}

identify_cell_in_tile <- function(df_spe_1, df_final_1, tile_shift, cores=10, prob=FALSE){
  
  df_tmp <- df_spe_1[, c("CenterX_local_px", "CenterY_local_px")]
  df_tmp <- as.list(as.data.frame(t(df_tmp)))
  names(df_tmp) <- df_spe_1$Barcode
  
  if(prob){
    
    cluster_tmp <- mclapply(df_tmp, function(x) {
      tmp <- find_tile(x, df_final_1, tile_shift)
      df_final_1[tmp, -c(1:4)]
    }, mc.cores=cores)
    
  } else {
    
    cluster_tmp <- mclapply(df_tmp, function(x) {
      tmp <- find_tile(x, df_final_1, tile_shift)
      df_final_1$Cluster[tmp]
    }, mc.cores=cores)
    
    names(cluster_tmp)<- names(df_tmp)
    cluster_tmp <- unlist(cluster_tmp)
    
  }
  
  return(cluster_tmp)
  
}




identify_cells_within_tile <- function(start_tile, df, tile_height){
  
  tile <- rbind(start_tile, start_tile+tile_height)
  
  df %>% dplyr::filter(x_local_px>=tile[1,1] & x_local_px<tile[2,1] & 
                         y_local_px>=tile[1,2] & y_local_px<tile[2,2]) %>%
    dplyr::count(Anno, .drop=FALSE) 
  
}

average_clusters_within_shift <- function(start_tile, clusters, tile_shift, tile_height){
  
  tile <- rbind(as.numeric(start_tile), as.numeric(start_tile)+tile_shift)
  
  clusters %>%
    dplyr::filter(tile[1,1]>=x & tile[2,1]<=x+tile_height & 
                         tile[1,2]>=y & tile[2,2]<=y+tile_height) %>%
    dplyr::count(cluster, .drop=FALSE) 
  
}


find_niches <- function(df_polygons, tile_height=1000, tile_shift=100, 
  image_size=c(4400, 4400), num_clusters=9, cores=10, prob=FALSE) {
  
  start=0-(tile_height-tile_shift)
  image_size=image_size+tile_height-tile_shift
    
  all_tiles <- expand.grid(seq(start,image_size[1],tile_shift), 
                           seq(start,image_size[2],tile_shift))
  all_tiles_list <- as.list(as.data.frame(t(all_tiles)))
  
  index <- split(1:nrow(df_polygons), df_polygons$fov)
  df_polygons_split <- lapply(index, function(x) df_polygons[x, ])  
  
  cells_in_tiles <- lapply(df_polygons_split, function(df) mclapply(all_tiles_list, function(x) 
    identify_cells_within_tile(x, df, tile_height)[,2], mc.cores=cores))
  # per fov find all cells within the defined tile (starting point defined in all_tiles)
  cells_in_tiles <- unlist(cells_in_tiles, recursive = FALSE)
  cells_in_tiles <- do.call(rbind, cells_in_tiles)
  colnames(cells_in_tiles) <- levels(df_polygons$Anno)
  
  all_tiles <- data.frame(x=rep(all_tiles[,1], length(index)), y=rep(all_tiles[,2], length(index)),
      fov=rep(names(index), each=nrow(all_tiles)))
  index_remove <- rowSums(cells_in_tiles) < 3  
  all_tiles <- all_tiles[!index_remove, ]  
  cells_in_tiles <- cells_in_tiles[!index_remove, ]  
  rownames(cells_in_tiles) <- paste0(all_tiles[,1], "_", all_tiles[,2], "_", all_tiles[,3])
  # find name for each tile 
  
  dist_tiles <- parDist(cells_in_tiles, method = 'bray', threads = cores)
  # calculate distance between tiles
  clusters <- hclust(dist_tiles, method="ward.D2")
  clusters <- cutree(clusters, k = num_clusters)
  # cluster tiles
  clusters_df <- as.data.frame(all_tiles)
  clusters_df$cluster <- as.factor(clusters)
  colnames(clusters_df) <- c("x", "y", "fov", "cluster")
  # make data frame for each tile starting point and cluster
  
  all_tiles$x <- as.numeric(all_tiles$x)
  all_tiles$y <- as.numeric(all_tiles$y)
  index <- split(1:nrow(all_tiles), all_tiles$fov)
  all_tiles_list <- lapply(index, function(x) all_tiles[x, ])
  all_tiles_list <- lapply(all_tiles_list, function(x) as.list(as.data.frame(t(x))))
  # split all tiles by fov
  clusters_df <- lapply(index, function(x) clusters_df [x, ])
  # split dataframe for each tile starting point and cluster by fov
  
  cluster_in_shift <- mapply(X=all_tiles_list, Y=clusters_df, function(X,Y) mclapply(X, function(x) 
    average_clusters_within_shift(x, Y, tile_shift, tile_height)[,2], mc.cores=cores))
  # find all tiles that overlap starting point and collect clusters do so by fov separately
  cluster_in_shift_x <- lapply(cluster_in_shift, function(x) do.call(rbind, x))
  
  if(prob){
    
    final_cluster <- lapply(cluster_in_shift_x, function(y) 
      t(apply(y, 1, function(x) x/sum(x))))
    # identify most common cluster for the tile do so by fov separately 
    all_tiles <- lapply(all_tiles_list, function(x) do.call(rbind, x))
    df_final <- mapply(X=all_tiles, Y=final_cluster, function(X,Y) 
      cbind(data.frame(x=as.numeric(X[,1]),
                 y=as.numeric(X[,2]), 
                 fov=X[,3], id=rownames(Y)), Y), SIMPLIFY = FALSE)
    df_final <- do.call(rbind, df_final)
    
    
  } else {
  
    final_cluster <- lapply(cluster_in_shift_x, function(y) 
      apply(y, 1, function(x) c(1:num_clusters)[which.max(x)]))
    # identify most common cluster for the tile do so by fov separately 
    all_tiles <- lapply(all_tiles_list, function(x) do.call(rbind, x))
    df_final <- mapply(X=all_tiles, Y=final_cluster, function(X,Y) 
      data.frame(x=as.numeric(X[,1]),
             y=as.numeric(X[,2]), 
             fov=X[,3],
             Cluster=as.factor(Y), id=names(Y)), SIMPLIFY = FALSE)
    df_final <- do.call(rbind, df_final)
    
  }
  
  return(df_final)
  
}




