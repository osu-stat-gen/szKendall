#' Calculate the within-cluster sum of squares percentage (WSS%)
#'
#' This function computes the within-cluster sum of squares percentage (WSS%), given the K-means clustering result.
#'
#' @param kmeans_object An object of class "`kmeans`", which is returned by the `kmeans()` function.
#' @return A numeric value of WSS%, between 0 and 1.
#' @examples
#' set.seed(1234)
#' km1 <- kmeans(t(scale(sim1.data)), centers=3, nstart=10, iter.max=100, algorithm="Hartigan-Wong")
#' WSS_percentage_kmeans(km1)
#' @export
WSS_percentage_kmeans <- function(kmeans_object){
  return(kmeans_object$tot.withinss/(kmeans_object$totss))
}

