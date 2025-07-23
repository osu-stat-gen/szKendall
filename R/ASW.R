#' Calculate the average silhouette width (ASW)
#'
#' This function computes the average silhouette width (ASW), given the clustering labels and the distance matrix on which the clustering was performed.
#'
#' @param clust_result A vector of length `n_single`, with each element indicating the index of cluster where the cell belongs.
#' @param distance_matrix  An `n_single`x`n_single` symmetric dissimilarity matrix that `cluster_result` is based on.
#' @return A numeric value of ASW, between -1 and 1.
#' @examples
#' set.seed(1234)
#' km1 <- kmeans(t(scale(sim1.data)), centers=3, nstart=10, iter.max=100, algorithm="Hartigan-Wong")
#' euclid.diss <- as.matrix(dist( t(scale(sim1.data)) ))
#' ASW(km1$cluster, euclid.diss)
#' @export
ASW <- function(clust_result, distance_matrix){
  ss <- cluster::silhouette(clust_result, dmatrix=distance_matrix)
  return(mean(ss[,3], na.rm=T))
}

