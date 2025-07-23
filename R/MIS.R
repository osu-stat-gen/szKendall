#' Calculate the minimum isolation score (MIS)
#'
#' This function computes the minimum isolation score (MIS), given the clustering labels and the distance matrix on which the clustering was performed.
#'
#' @param clust_result A vector of length `n_single`, with each element indicating the index of cluster where the cell belongs.
#' @param distance_matrix  An `n_single`x`n_single` symmetric dissimilarity matrix that `cluster_result` is based on.
#' @return A numeric value of MIS across all clusters, excluding the clusters that contain only a single cell (as their isolation score is infinity).
#' @examples
#' set.seed(1234)
#' km1 <- kmeans(t(scale(sim1.data)), centers=3, nstart=10, iter.max=100, algorithm="Hartigan-Wong")
#' euclid.diss <- as.matrix(dist( t(scale(sim1.data)) ))
#' MIS(km1$cluster, euclid.diss)
#' @export
MIS <- function(clust_result, distance_matrix){

  K <- length(unique(clust_result))
  isol <- rep(0, K)

  # the isolation score is calculated for each cluster. The larger the better.
  # when isolation score > 1, the cluster is a L(or L*) cluster defined in PAM.
  # when isolation score < 1, the cluster is a non-isolated cluster defined in PAM.
  for(k in 1:K){
    within <- distance_matrix[clust_result == k, clust_result == k]
    diameter <- max(within)

    between <- distance_matrix[clust_result == k, clust_result != k]
    separation <- min(between)

    isol[k] <- separation/diameter
  }

  # isol[is.infinite(isol)] <- NA
  # isol <- sort(isol, decreasing=TRUE, na.last=TRUE)

  MIS <- min(isol, na.rm=T)

  return(MIS)
}

