#' Heatmap visualization of a dissimilarity matrix
#'
#' This function visualizes a dissimilarity matrix as a heatmap.
#'
#' @importFrom plsgenomics matrix.heatmap
#' @param data  An `n_single`x`n_single` symmetric dissimilarity matrix.
#' @param title The title given to the output plot. Default is empty ("").
#' @return A plot that visualizes the dissimilarity matrix as heatmap. Each element in the heatmap is scaled to be between 0 and 1 (by subtracting the minimum and dividing the range).
#' @examples
#' scHiC_hm_dist(euclid.diss, title="Euclid")
#' @export
scHiC_hm_dist <- function(data, title = ""){

  normmatr <- data
  normmatr[lower.tri(normmatr, diag=TRUE)] <- 0
  maxvalue <- max(data[upper.tri(data, diag=FALSE)])
  minvalue <- min(data[upper.tri(data, diag=FALSE)])
  normmatr[upper.tri(normmatr)] <- (data[upper.tri(data)] - minvalue)/(maxvalue - minvalue)
  normmatr <- normmatr+t(normmatr)

  return(matrix.heatmap(normmatr, main = title, legend.width=0.8, legend.mar=3.1))
}

