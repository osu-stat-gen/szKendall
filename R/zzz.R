# R/zzz.R

#' szKendall: Structural-zero-aware Kendall's tau for scHi-C data
#' 
#' @docType package
#' @name szKendall
#' @useDynLib szKendall, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @keywords internal
"_PACKAGE"