library(vegan)
library(Seurat)
library(tidyverse)
library(parallel)
library(pals)



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


find_niches <- function(df_polygons, tile_height=1000, tile_shift=100, image_size=4400, num_clusters=9) {
  
  all_tiles <- expand.grid(seq(0,image_size,tile_shift), seq(0,image_size,tile_shift))
  all_tiles_list <- as.list(as.data.frame(t(all_tiles)))
  
  index <- split(1:nrow(df_polygons), df_polygons$fov)
  df_polygons_split <- lapply(index, function(x) df_polygons[x, ])  
  
  cells_in_tiles <- lapply(df_polygons_split, function(df) mclapply(all_tiles_list, function(x) 
    identify_cells_within_tile(x, df, tile_height)[,2], mc.cores=10))
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
  
  dist_tiles <- vegdist(cells_in_tiles)
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
    average_clusters_within_shift(x, Y, tile_shift, tile_height)[,2], mc.cores=10))
  # find all tiles that overlap starting point and collect clusters do so by fov separately
  cluster_in_shift_x <- lapply(cluster_in_shift, function(x) do.call(rbind, x))
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
  
  return(df_final)
  
}




