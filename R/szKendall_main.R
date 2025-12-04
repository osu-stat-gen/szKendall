utils::globalVariables(c(
  "foreach", "%dopar%", "para",
  "szkendall.dist"
))

#' Calculate szKendall-based distance between two sets of contact maps.
#'
#' @description
#' This wrapper function harmonizes the input formats of \code{sim.data}
#' and \code{true.data}, and computes the pairwise dissimilarity using
#' the selected szKendall method.
#' 
#' @useDynLib szKendall, .registration = TRUE
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @importFrom Rcpp sourceCpp 
#' 
#' @param sim.data Input data corresponding to simulated contact maps.
#'   It can be a list of \eqn{\#LP \times \#LP} matrices (one per cell),
#'   a 3D array of dimension \eqn{\#LP \times \#LP \times \#cells},
#'   or a \eqn{(\#LP) \times (\#cells)} matrix. 
#'   See \code{\link{bind_to_matrix}} for details.
#' @param true.data Input data corresponding to true contact maps, in the same format as \code{sim.data}.
#' @param dist.method Character string specifying which szKendall dissimilarity measure to use.
#'   Supported options: \code{"szKendall"} (default), \code{"szKendall1"}. See @detail below.
#'
#' @return A square symmetric szKendall dissimilarity matrix, where the dimension is the number of single cells.
#'
#' @details
#' Internally, this function first converts the given data structures into consistent \eqn{\#LP \times \#cells} matrix form via
#' \code{bind_to_matrix()}, and then dispatches to the corresponding szKendall dissimilarity function.
#' 
#' #' When \code{dist.method = "szKendall"}, pairwise contributions are weighted
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
#' szKendall.diss(sim1.data, true1.data)
#' 
#' 
#' @export
szKendall <- function(sim.data, true.data,
                      dist.method = "szKendall") {

  # Normalize input types into matrix form
  mat_data <- bind_to_matrix(sim.data, true.data)
  sim      <- mat_data$sim
  truth    <- mat_data$true
  rm(mat_data)

  # Map method names to computation functions
  method_dispatch <- list(
    szKendall  = function(sim, truth) szKendall.diss(sim, truth),
    szKendall1 = function(sim, truth) szKendall.diss2(sim, truth)
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

  # Execute selected method
  res <- method_dispatch[[dist.method]](sim, truth)
  return(res)
}