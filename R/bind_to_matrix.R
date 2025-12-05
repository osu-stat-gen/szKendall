#' Bind / Bind array or list
#'
#' @description 
#' This function constructs a \eqn{(\#LP) \times (\#cells)} matrix, where each column corresponds to the upper‐triangular entries (excluding the diagonal) of a contact map for a single cell.
#'
#' @param counts  An array or a list or a matrix object. See detail.
#'
#' @return  A \eqn{(\#LP) \times (\#cells)} matrix. 
#' 
#' @examples
#' ## ---------------------------------------------------------
#' ## Example 1: sim.data and true.data as lists
#' ## ---------------------------------------------------------
#' ## generate an example of list type counts
#' nLP   <- 10   # Number of LPs
#' nCell <- 5    # Number of cells
#'
#' counts_list  <- vector("list", nCell)
#'
#' for (i in seq_len(nCell)) {
#'   counts_list[[i]]  <- matrix(0, nLP, nLP)
#'   idx <- upper.tri(counts_list[[i]], diag = FALSE)
#'   counts_list[[i]][idx]  <- seq_len(sum(idx))
#' }
#' 
#' ## run 
#' bind_to_matrix(counts_list)
#'
#'
#' ## ---------------------------------------------------------
#' ## Example 2: sim.data and true.data as 3D arrays
#' ## ---------------------------------------------------------
#' ## generate an example of array type counts
#' sim.mat  <- matrix(0, nLP, nLP)
#' idx      <- upper.tri(sim.mat, diag = FALSE)
#' sim.mat[idx]  <- seq_len(sum(idx))
#' counts_array  <- array(sim.mat,  dim = c(nLP, nLP, nCell))  # 3D array
#'
#' ## run 
#' bind_to_matrix(counts_array)
#'
#'
#' ## ---------------------------------------------------------
#' ## Example 3: sim.data and true.data already in matrix form
#' ##  Dimension = (#LP*(#LP-1)/2) x (#cells)
#' ## ---------------------------------------------------------
#' ## generate an example of matrix type counts
#' counts_matrix  <- matrix(seq_len(sum(upper.tri(matrix(0, nLP, nLP)))),
#'                     nrow = sum(upper.tri(matrix(0, nLP, nLP))),
#'                     ncol = nCell)
#' 
#' ## run 
#' bind_to_matrix(counts_matrix)
#' 
#' @details
#' When `counts` is a list, it must satisfy:
#' \itemize{
#'   \item \code{length(counts)} equals the number of cells.
#'   \item Each element \code{x[[i]]} is a square numeric matrix of size
#'         (\# of LP) \eqn{\times} (\# of LP), representing the LP-by-LP
#'         interaction matrix for the \eqn{i}-th cell.
#' }
#' 
#' When `counts` is a 3D array, it must satisfy:
#' \itemize{
#'   \item Dimension: \eqn{(\#LP) \times (\#LP) \times (\#cells)}.
#'   \item For every slice \code{[ , , i]}, the matrix is square of size
#'         \eqn{\#LP \times \#LP}, representing the LP-by-LP contact map
#'         of the \eqn{i}-th cell.
#' }
#' 
#' When `counts` is a already matrix of dimension \eqn{\#LP \times \#cells}, no further conversion is applied.
#' 
#' @export
bind_to_matrix <- function(counts) {
  allowed_classes <- c("array", "list", "matrix")
  counts_class  <- class(counts)[1]

  if (!counts_class %in% allowed_classes) {
    stop(
      sprintf(
        "Unsupported input class '%s'.\nAllowed classes: %s",
        counts_class,
        paste(allowed_classes, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  upper_tri_cbind_array <- function(arr,diag=FALSE){
    id_tmp <- upper.tri(arr[,,1], diag = FALSE)
    if (length(dim(arr)) != 3) {
        stop("`arr` must be a 3D array (n x n x k).", call. = FALSE)
    }
    if(dim(arr)[1] != dim(arr)[2]){
        stop("First two dimensions of `arr` must be equal (n x n).", call. = FALSE)
    }
    cols <- lapply(seq_len(dim(arr)[3]),function(x){arr[,,x][id_tmp]})
    res <- do.call(cbind,cols)
    return(res)       
  }

  upper_tri_cbind_list <- function(lst,diag=FALSE){
    if (!all(vapply(lst, is.matrix, logical(1)))) {
        stop("All elements of `mat_list` must be matrices.", call. = FALSE)
        }
    # check dimension
    dims <- lapply(lst, dim)
    if (!all(vapply(dims, function(x) x[1] == x[2], logical(1)))) {
        stop("All matrices must be square (n x n).", call. = FALSE)
        } 
    n <- dims[[1]][1]
    if (!all(vapply(dims, function(x) x[1] == n, logical(1)))) {
        stop("All matrices must have the same dimension.", call. = FALSE)
        }
    
    id_tmp <- upper.tri(lst[[1]], diag = diag)
    cols <- lapply(lst, function(x) x[id_tmp])
    res <- do.call(cbind, cols)
    return(res)      
  }
  
  res <- switch(counts_class,
    "matrix" = counts,
    "array"  = upper_tri_cbind_array(counts),
    "list"   = upper_tri_cbind_list(counts) 
  )
  return(res)
}