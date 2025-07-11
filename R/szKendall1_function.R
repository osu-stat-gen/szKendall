
utils::globalVariables(c(
  "foreach", "%dopar%", "para",
  "szkendall1.dist"
))

#' Calculate szKendall1 dissimilarity
#'
#' This function computes the szKendall1 dissimilarity matrix given an "observed" locus-pair by cell single-cell Hi-C matrix and the "true" expected contact count matrix (where only structural zero positions have the value zero).
#'
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @useDynLib szKendall, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param sim.data A simulated or observed single-cell Hi-C matrix, where the rows represent locus pairs and columns represent cells.
#' @param true.data The expected or imputed single-cell Hi-C matrix, which has the same dimension as sim.data.
#' @return A square szKendall1 dissimilarity matrix, where the dimension is the number of single cells.
#' @examples
#' foreach::registerDoSEQ()
#' szKendall1.diss(sim1.data, true1.data)
#' @export
szKendall1.diss <- function(sim.data, true.data){

  # Register a parallel backend using the following lines if it is not done first:
  if (!foreach::getDoParRegistered()) {
    numCores <- max(1, parallel::detectCores() - 1)
    doParallel::registerDoParallel(cores = numCores)
    on.exit(doParallel::stopImplicitCluster())  # Clean up
  }

  n.cells <- ncol(sim.data)
  n1 <- ceiling(sqrt(2*nrow(sim.data)))
  weight_vec1 <- cal_weight_vec1(n=n1)
  weight_sz_vec1 <- cal_weight_sz_vec1(n=n1)

  szkendall1.dist <- matrix(0, nrow=n.cells, ncol=n.cells)

  fail.index <- c()
  k <- 0
  for(i in 1:(n.cells-1)){
    for(j in (i+1):(n.cells)){
      if(szkendall1.dist[i,j]<1e-5){
        k <- k+1
        fail.index <- rbind(fail.index, c(i,j))
      }
    }
  }

  score1 <- foreach(para=1:nrow(fail.index), .combine = 'rbind', .noexport="szkendall1") %dopar%{
    value <- szkendall1(sim.data[,fail.index[para,1]], sim.data[,fail.index[para,2]], which(true.data[,fail.index[para,1]]==0), which(true.data[,fail.index[para,2]]==0), weight_vec1, weight_sz_vec1, type="Nodiag")
    c(fail.index[para,1], fail.index[para,2], value)
  }

  for(m in 1:nrow(fail.index)){
    szkendall1.dist[score1[m,1], score1[m,2]] <- score1[m,3]
  }

  fail.index <- c()
  k <- 0
  for(i in 1:(n.cells-1)){
    for(j in (i+1):(n.cells)){
      if(szkendall1.dist[i,j]<1e-5){
        k <- k+1
        fail.index <- rbind(fail.index, c(i,j))
      }
    }
  }

  if(!is.null(dim(fail.index))) {
    warning("Parallel computing of szKendall1 among the cell pairs is not done properly. Please re-run the function and/or re-set the parallel backend.")
  }

  szkendall1.dist <- (szkendall1.dist + t(szkendall1.dist))

  rownames(szkendall1.dist) <- NULL
  colnames(szkendall1.dist) <- NULL

  return(szkendall1.dist)
}


