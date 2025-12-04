#' Bind / Bind array or list
#'
#' @description 
#' This function constructs a \eqn{(\#LP) \times (\#cells)} matrix, where each column corresponds to the upper‐triangular entries (excluding the diagonal) of a contact map for a single cell.
#'
#' @param sim.data  An array or a list or a matrix object. See detail.
#' @param true.data An array or a list or a matrix object. See detail.
#'
#' @return A list with two elements: \code{sim.data} and \code{true.data}, matrices. 
#' 
#' @examples
#' ## ---------------------------------------------------------
#' ## Example 1: sim.data and true.data as lists
#' ## ---------------------------------------------------------
#' ## generate sim.data and true.data
#' nLP   <- 10   # Number of LPs
#' nCell <- 5    # Number of cells
#'
#' sim.data  <- vector("list", nCell)
#' true.data <- vector("list", nCell)
#'
#' for (i in seq_len(nCell)) {
#'   sim.data[[i]]  <- matrix(0, nLP, nLP)
#'   idx <- upper.tri(sim.data[[i]], diag = FALSE)
#'   sim.data[[i]][idx]  <- seq_len(sum(idx))
#'
#'   true.data[[i]] <- sim.data[[i]] * 10
#' }
#' 
#' ## run 
#' bind_to_matrix(sim.data, true.data)
#'
#'
#' ## ---------------------------------------------------------
#' ## Example 2: sim.data and true.data as 3D arrays
#' ## ---------------------------------------------------------
#' ## generate sim.data and true.data
#' sim.mat  <- matrix(0, nLP, nLP)
#' idx      <- upper.tri(sim.mat, diag = FALSE)
#' sim.mat[idx]  <- seq_len(sum(idx))
#'
#' true.mat <- sim.mat * 10
#'
#' sim.data  <- array(sim.mat,  dim = c(nLP, nLP, nCell))  # 3D array
#' true.data <- array(true.mat, dim = c(nLP, nLP, nCell))
#'
#' ## run 
#' bind_to_matrix(sim.data, true.data)
#'
#'
#' ## ---------------------------------------------------------
#' ## Example 3: sim.data and true.data already in matrix form
#' ##  Dimension = (#LP*(#LP-1)/2) x (#cells)
#' ## ---------------------------------------------------------
#' ## generate sim.data and true.data
#' sim.data  <- matrix(seq_len(sum(upper.tri(matrix(0, nLP, nLP)))),
#'                     nrow = sum(upper.tri(matrix(0, nLP, nLP))),
#'                     ncol = nCell)
#' true.data <- sim.data * 10
#' 
#' ## run 
#' bind_to_matrix(sim.data, true.data)
#' 
#' @details
#' When `sim.data` and `true.data` are lists, they must satisfy:
#' \itemize{
#'   \item \code{length(sim.data)} equals the number of cells.
#'   \item Each element \code{x[[i]]} is a square numeric matrix of size
#'         (\# of LP) \eqn{\times} (\# of LP), representing the LP-by-LP
#'         interaction matrix for the \eqn{i}-th cell.
#' }
#' 
#' When `sim.data` and `true.data` are 3D arrays, they must satisfy:
#' \itemize{
#'   \item Dimension: \eqn{(\#LP) \times (\#LP) \times (\#cells)}.
#'   \item For every slice \code{[ , , i]}, the matrix is square of size
#'         \eqn{\#LP \times \#LP}, representing the LP-by-LP contact map
#'         of the \eqn{i}-th cell.
#' }
#' 
#' When `sim.data` and `true.data` are already matrices of dimension \eqn{\#LP \times \#cells}, no further conversion is applied.
#' 
#' @export
bind_to_matrix <- function(sim.data, true.data) {
  allowed_classes <- c("array", "list", "matrix")
  sim_class  <- class(sim.data)[1]
  true_class <- class(true.data)[1]

  if (!identical(sim_class, true_class)) {
    stop(
      sprintf(
        "`sim.data` and `true.data` must have the same class.\n  sim.data:  '%s'\n  true.data: '%s'",
        sim_class, true_class
      ),
      call. = FALSE
    )
  } 

  if (!sim_class %in% allowed_classes) {
    stop(
      sprintf(
        "Unsupported input class '%s'.\nAllowed classes: %s",
        sim_class,
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
  
  
  
  sim.data <- switch(sim_class,
    "matrix" = sim.data,
    "array"  = upper_tri_cbind_array(sim.data),
    "list"   = upper_tri_cbind_list(sim.data) 
  )
  true.data <- switch(true_class,
    "matrix" = true.data,
    "array"  = upper_tri_cbind_array(true.data),
    "list"   = upper_tri_cbind_list(true.data) 
  )

  out <- list(
    sim  = sim.data,
    truth = true.data
  )
  return(out)
}