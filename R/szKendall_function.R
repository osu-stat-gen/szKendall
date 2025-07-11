
utils::globalVariables(c(
  "foreach", "%dopar%", "para",
  "szkendall.dist"
))

#' Calculate szKendall dissimilarity
#'
#' This function computes the szKendall dissimilarity matrix given an "observed" locus-pair by cell single-cell Hi-C matrix and the "true" expected contact count matrix (where only structural zero positions have the value zero).
#'
#' @importFrom foreach foreach %dopar% registerDoSEQ
#' @useDynLib szKendall, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param sim.data A simulated or observed single-cell Hi-C matrix, where the rows represent locus pairs and columns represent cells.
#' @param true.data The expected or imputed single-cell Hi-C matrix, which has the same dimension as sim.data.
#' @return A square szKendall dissimilarity matrix, where the dimension is the number of single cells.
#' @examples
#' foreach::registerDoSEQ()
#' szKendall.diss(sim1.data, true1.data)
#' @export
szKendall.diss <- function(sim.data, true.data){

  # Register a parallel backend using the following lines if it is not done first:
  if (!foreach::getDoParRegistered()) {
    numCores <- max(1, parallel::detectCores() - 1)
    doParallel::registerDoParallel(cores = numCores)
    on.exit(doParallel::stopImplicitCluster())  # Clean up
  }

  n.cells <- ncol(sim.data)
  n1 <- ceiling(sqrt(2*nrow(sim.data)))
  weight_vec3 <- cal_weight_vec3(n=n1)
  weight_sz_vec3 <- cal_weight_sz_vec3(n=n1)

  szkendall.dist <- matrix(0, nrow=n.cells, ncol=n.cells)

  fail.index <- c()
  k <- 0
  for(i in 1:(n.cells-1)){
    for(j in (i+1):(n.cells)){
      if(szkendall.dist[i,j]<1e-5){
        k <- k+1
        fail.index <- rbind(fail.index, c(i,j))
      }
    }
  }

  score3 <- foreach(para=1:nrow(fail.index), .combine = 'rbind', .noexport="szkendall") %dopar%{
    value <- szkendall(sim.data[,fail.index[para,1]], sim.data[,fail.index[para,2]], which(true.data[,fail.index[para,1]]==0), which(true.data[,fail.index[para,2]]==0), weight_vec3, weight_sz_vec3, type="Nodiag")
    c(fail.index[para,1], fail.index[para,2], value)
  }

  for(m in 1:nrow(fail.index)){
    szkendall.dist[score3[m,1], score3[m,2]] <- score3[m,3]
  }

  fail.index <- c()
  k <- 0
  for(i in 1:(n.cells-1)){
    for(j in (i+1):(n.cells)){
      if(szkendall.dist[i,j]<1e-5){
        k <- k+1
        fail.index <- rbind(fail.index, c(i,j))
      }
    }
  }

  if(!is.null(dim(fail.index))) {
    warning("Parallel computing of szKendall among the cell pairs is not done properly. Please re-run the function and/or re-set the parallel backend.")
  }

  szkendall.dist <- (szkendall.dist + t(szkendall.dist))

  rownames(szkendall.dist) <- NULL
  colnames(szkendall.dist) <- NULL

  return(szkendall.dist)
}


