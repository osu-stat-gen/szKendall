#' Calculate szKendall-based dissimilarity among a set of single cells using scHi-C data.
#'
#' @description
#' This function takes scHi-C data input and computes dissimilarity between each pair of single cells.
#' In addition to the selected szKendall method, dissimilarity based on Euclidean and Kendall will also be provided for comparison. 
#' 
#' 
#' @param contacts Single-cell Hi-C input data.This package supports three input formats:
#' \itemize{
#'  \item a list of scHi-C 2D matrices of size \eqn{n} (e.g. the number of bins in a chromosome);
#'  \item a 3D array of dimension \eqn{n \times n \times \#cells};
#'  \item a \eqn{\#locus\text{-}pairs \times \#cells} matrix.
#' }
#' 
#' @param dist.method Character string specifying which szKendall dissimilarity measure to use.
#'   Supported options: \code{"szKendall"} (default), \code{"szKendall1"}. See @detail below.
#' 
#'
#' @return A list of three squared symmetric dissimilarity matrices will be returned, one for szKendall(1),
#'  one for Euclidean, and one for Kendall (the last two was for the comparison purpose when drawing the heatmaps).
#' 
#' @section Warning:
#' \itemize{
#'   \item Register a parallel backend if it is not done yet. Other wise this function uses \code{NumCores() - 1} cores by default.
#'   \item 
#' }
#'
#' @seealso \code{\link{bind_to_matrix}},
#'          \code{\link{szKendall.diss}},
#'          \code{\link{szKendall.diss1}}
#' 
#' @details
#' This function first calls \code{bind_to_matrix()} to prepare the input data into the \eqn{\#locus\text{-}pairs \times \#cells} matrix format if necessary.
#' Then it calls the function corresponding to the \code{dist.method} specified.
#' 
#' When \code{dist.method = "szKendall"}, pairwise contributions are weighted
#' according to the genomic distance difference using
#' \deqn{
#'   W^{K}_{\,|\!|j-i|\!|-|\!|v-u|\!|} =
#'   \left( n - 1 - \bigl|\,|j-i| - |v-u|\,\bigr| \right)^{-0.4}
#' }
#' where \eqn{i,j} index loci of one contact pair and \eqn{u,v} index loci of the
#' other contact pair.
#' 
#' For \code{dist.method = "szKendall1"}, 
#' \deqn{
#'   W^{K1}_{\,|\!|j-i|\!|-|\!|v-u|\!|} = 1
#' }
#' 
#' 
#'
#' @examples
#' foreach::registerDoSEQ()
#' 
#' 
#' @export
szKendall <- function(counts,
                      dist.method = "szKendall") {

  # Validate selected method
  if (!dist.method %in% names(method_dispatch)) {
    stop(
      sprintf(
        "Unsupported dist.method: '%s'.\nAvailable options: %s",
        dist.method,
        paste(names(method_dispatch), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Normalize input types into matrix form
  mat_count <- bind_to_matrix(counts)

  # Compute kendall.dist
  kendall.dist <- kendall.diss(count_mat)

  # Compute euclidean.dist; TODO parallelize this
  euclidean.dist <- as.matrix(dist(t(count_mat),method="euclidean"))
  # #Standardize?
  # euclidean.dist <- as.matrix(dist(t(scale(count_mat)),method="euclidean"))
  
  # Compute szkendall.dist (NOTE kendall.dist is reused if dist.method=="szkendall1")
  method_dispatch <- list(
    szKendall  = function(mat_count,...) szKendall.diss(mat_count,...),
    szKendall1 = function(mat_count,...) szKendall1.diss(mat_count,...)
  )
  szkendall.dist <- method_dispatch[[dist.method]](mat_count, kendall.dist)
  
  res=list(
    szKendall = szkendall.dist,
    Euclid = euclidean.dist,
    Kendall = kendall.dist
  )
  return(res)
}