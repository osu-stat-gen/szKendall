#' Calculate the within-cluster sum of squares percentage (WSS%)
#'
#' This function computes the within-cluster sum of squares percentage (WSS%), given the clustering labels (from any clustering algorithm) and the distance matrix on which the clustering was performed.
#'
#' @param clust_result A vector of length `n_single`, with each element indicating the index of cluster where the cell belongs.
#' @param distance_matrix  An `n_single`x`n_single` symmetric dissimilarity matrix that `cluster_result` is based on.
#' @return A numeric value of WSS%, between 0 and 1.
#' @examples
#' set.seed(1234)
#' km1 <- kmeans(t(scale(sim1.data)), centers=3, nstart=10, iter.max=100, algorithm="Hartigan-Wong")
#' euclid.diss <- as.matrix(dist( t(scale(sim1.data)) ))
#' WSS_percentage(km1$cluster, euclid.diss)
#' @export
WSS_percentage <- function(clust_result, distance_matrix){
  K <- length(unique(clust_result))
  if(K == 1){
    return(1)
  } else {
    within <- c()
    distance_matrix <- as.matrix(distance_matrix)
    for(k in 1:K){
      clustk <- distance_matrix[clust_result == k, clust_result == k]
      within <- c(within, sum(clustk^2)/(2*sum(clust_result == k)))
    }

    return(sum(within)/(sum(distance_matrix^2)/(2*nrow(distance_matrix))))

  }
}

