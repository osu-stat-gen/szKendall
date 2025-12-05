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
#' @param heatmap A logical. If TRUE, draw heatmaps of
#'   (i) the szKendall dissimilarity,
#'   (ii) the Euclidean distance, and
#'   (iii) the Kendall dissimilarity,
#'   in that order. Default is FALSE. Note that this may take
#'   some time if the number of cells is large.
#' 
#'
#' @return A list of three squared symmetric dissimilarity matrices will be returned, one for szKendall(1),
#'  one for Euclidean, and one for Kendall (the last two was for the comparison purpose when drawing the heatmaps).
#' 
#' @section Warning: Register a parallel backend if it is not done yet. Other wise this function uses \code{NumCores() - 1} cores by default.
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
#' res <- szKendall(szKendall::PFC_subset_matrix,dist.method="szKendall")
#' # or res <- szKendall(szKendall::PFC_subset_array,dist.method="szKendall")
#' # or res <- szKendall(szKendall::PFC_subset_list,dist.method="szKendall")
#' 
#' 
#' 
#' @export
szKendall <- function(contacts, dist.method = "szKendall", heatmap=FALSE) {

  # Normalize input types into matrix form
  mat_count <- bind_to_matrix(contacts)

  # Compute kendall.dist
  kendall.dist <- kendall.diss(mat_count)

  # Compute euclidean.dist; TODO parallelize this
  euclidean.dist <- as.matrix(dist(t(scale(mat_count)),method="euclidean"))
  
  # Compute szkendall.dist (NOTE kendall.dist is reused if dist.method=="szkendall1")
  method_dispatch <- list(
    szKendall  = function(mat_count,...) szKendall.diss(mat_count,...),
    szKendall1 = function(mat_count,...) szKendall1.diss(mat_count,...)
  )
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
  szkendall.dist <- method_dispatch[[dist.method]](mat_count, kendall.dist)

  if(heatmap){
    par(mfrow=c(1,3), mar = c(4,4,2,5), oma = c(0,0,0,2))
    scHiC_hm_dist(res[[1]],title=dist.method)
    scHiC_hm_dist(res[[2]],title="Euclidean")
    scHiC_hm_dist(res[[3]],title="Kendall")
  }
  
  res=list(
    szKendall = szkendall.dist,
    Euclid = euclidean.dist,
    Kendall = kendall.dist
  )
  return(res)
}