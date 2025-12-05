#' Calculate szKendall dissimilarity
#'
#' This function computes the szKendall dissimilarity matrix given a  (locus-pair by cell single-cell Hi-C) matrix
#' @param mat_count A \eqn{(\#LP) \times (\#cells)} matrix.  
#' @return A square szKendall dissimilarity matrix, where the dimension is the number of single cells.
#' @examples
#' foreach::registerDoSEQ()
#' szKendall.diss(mat_count)
#' @export
szKendall.diss <- function(mat_count,...){

  # Register a parallel backend using the following lines if it is not done first:
  if (!foreach::getDoParRegistered()) {
    numCores <- max(1, parallel::detectCores() - 1)
    doParallel::registerDoParallel(cores = numCores)
    on.exit(doParallel::stopImplicitCluster())  # Clean up
  }
  
  n.cells <- ncol(mat_count)
  n1 <- ceiling(sqrt(2*nrow(mat_count)))
  weight_vec <- cal_weight_vec(n=n1)
  weight_sz_vec <- cal_weight_sz_vec(n=n1)

  szkendall.dist <- matrix(0, nrow=n.cells, ncol=n.cells)

  fail.index <-  matrix(ncol = 2, nrow = 0)
  for(i in 1:(n.cells-1)){
    for(j in (i+1):(n.cells)){
      if(szkendall.dist[i,j]<1e-5){
        fail.index <- rbind(fail.index, c(i,j))
      }
    }
  }

  scoreS <- foreach(para = 1:nrow(fail.index), .combine = 'rbind', .noexport = "szkendall_cpp") %dopar% {
                        i <- fail.index[para, 1]; j <- fail.index[para, 2]
                        szi <- which(mat_count[, i] == 0)
                        szj <- which(mat_count[, j] == 0)
                        val <- szkendall_cpp(mat_count[, i],mat_count[, j], szi, szj, weight_vec, weight_sz_vec, type = "Nodiag")
                        c(i, j, val)
                      }

  for(m in 1:nrow(fail.index)){
    szkendall.dist[scoreS[m, 1], scoreS[m, 2]] <- scoreS[m, 3]
  }

  
  if(any(szkendall.dist[upper.tri(szkendall.dist,diag=F)] < 1e-5)) {
    warning("Parallel computing of szKendall among the cell pairs is not done properly. Please re-run the function and/or re-set the parallel backend.")
  }

  szkendall.dist <- (szkendall.dist + t(szkendall.dist))

  rownames(szkendall.dist) <- NULL
  colnames(szkendall.dist) <- NULL

  return(szkendall.dist)
}


