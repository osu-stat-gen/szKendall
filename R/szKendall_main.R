#' Calculate szKendall-based dissimilarity among a set of single cells using scHi-C data.
#'
#' @description
#' This function takes scHi-C data input and computes dissimilarity between each pair of single cells.
#' In addition to the selected szKendall method, dissimilarity based on Euclidean and Kendall will also be provided for comparison. 
#' 
#' 
#' @param contacts Single-cell Hi-C input data.This package supports three input formats:
#' \itemize{
#'  \item \strong{(1)} a list of scHi-C 2D matrices of size \dqn{n} (e.g. the number of bins in a chromosome);
#'  \item \strong{(2)} a 3D array of dimension \deqn{n \times n \times \#cells};
#'  \item \strong{(3)} a \deqn{\#locus\text{-}pairs \times \#cells} matrix.
#' }
#' 
#' @param dist.method Character string specifying which szKendall dissimilarity measure to use.
#'   Supported options: \code{"szKendall"} (default), \code{"szKendall1"}. See @detail below.
#'
#' @return A list of three squared symmetric dissimilarity matrices will be returned, one for szKendall(1),
#'  one for Euclidean, and one for Kendall (the last two was for the comparison purpose when drawing the heatmaps).
#'
#' @details
#' This function first calls \code{bind_to_matrix()} to prepare the input data into the \deqn{\#locus\text{-}pairs \times \#cells} matrix format if necessary.
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
#' @seealso \code{\link{bind_to_matrix}},
#'          \code{\link{szKendall.diss}},
#'          \code{\link{szKendall.diss2}}
#'
#' @examples
#' foreach::registerDoSEQ()
#' 
#' 
#' @export
szKendall <- function(counts,
                      dist.method = "szKendall") {

  # Normalize input types into matrix form
  mat_count <- bind_to_matrix(counts)

  # Map method names to computation functions
  method_dispatch <- list(
    szKendall  = function(mat_count) szKendall.diss(mat_count),
    szKendall1 = function(mat_count) szKendall2.diss(mat_count)
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

  # TODO: add Euclidean and Kendall
  res <- method_dispatch[[dist.method]](sim, truth)
  return(res)
}