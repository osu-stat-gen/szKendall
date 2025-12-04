#' Calculate the adjusted Rand index (ARI)
#'
#' This function computes the Adjusted Rand Index (ARI) between two clustering label vectors.
#'
#' @importFrom mclust adjustedRandIndex
#' @param true_label A vector of length `n_single`, with each element indicating the true class (cell type or subtype) of each observation.
#' @param clust_result A vector of length `n_single`, with each element indicating the index of the cluster where the corresponding observation is assigned.
#' @return A numeric value of the Adjusted Rand Index (ARI), bounded above by 1
#' @examples
#' set.seed(1234)
#' n_single <- 90
#' true_label <- factor(rep(1:3, each = n_single/3))
#' km1 <- kmeans(matrix(rnorm(n_single*10), nrow = n_single), centers = 3,
#'               nstart = 10, iter.max = 100, algorithm = "Hartigan-Wong")
#' ARI(true_label, km1$cluster)
#' @export
ARI <- function(true_label, clust_result) {
  mclust::adjustedRandIndex(true_label, clust_result)
}