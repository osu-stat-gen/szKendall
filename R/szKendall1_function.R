#' Calculate szKendall1 dissimilarity
#'
#' This function computes the szKendall1 dissimilarity matrix given a  (locus-pair by cell single-cell Hi-C) matrix
#' @param mat_count A \eqn{(\#LP) \times (\#cells)} matrix.  
#' @return A square szKendall1 dissimilarity matrix, where the dimension is the number of single cells.
#' @examples
#' foreach::registerDoSEQ()
#' szKendall.diss(mat_count)
#' @export
szKendall1.diss <- function(mat_count,kendall.dist,...){

  # Register a parallel backend using the following lines if it is not done first:
  if (!foreach::getDoParRegistered()) {
    numCores <- max(1, parallel::detectCores() - 1)
    doParallel::registerDoParallel(cores = numCores)
    on.exit(doParallel::stopImplicitCluster())  # Clean up
  }

  n.cells <- ncol(mat_count)
  n1 <- ceiling(sqrt(2*nrow(mat_count)))
  weight_sz_vec1 <- cal_weight_sz_vec1(n=n1)

  szkendall1.dist <- matrix(0, nrow=n.cells, ncol=n.cells)

  fail.index <- matrix(ncol = 2, nrow = 0)
  for(i in 1:(n.cells-1)){
    for(j in (i+1):(n.cells)){
      if(szkendall1.dist[i,j]<1e-5){
        fail.index <- rbind(fail.index, c(i,j))
      }
    }
  }

  scoreS1 <- foreach(para = 1:nrow(fail.index), .combine = 'rbind',.noexport = "szkendall1_cpp") %dopar% {
                        i <- fail.index[para, 1]; j <- fail.index[para, 2]
                        szi <- which(mat_count[, i] == 0)
                        szj <- which(mat_count[, j] == 0)
                        val <- szkendall1_cpp_SZpart(mat_count[, i],
                                              mat_count[, j],
                                              szi, szj, weight_sz_vec1, type = "Nodiag")
                        c(i, j, val)
                      }

  for(m in 1:nrow(fail.index)){
    szkendall1.dist[scoreS1[m,1], scoreS1[m,2]] <- scoreS1[m,3]
  }

  

  szkendall1.dist <- (szkendall1.dist + t(szkendall1.dist)) # symmetrize

  szkendall1.dist <- (szkendall1.dist + kendall.dist)      # add Kendall part
  
  if(any(szkendall1.dist[upper.tri(szkendall1.dist,diag=F)] < 1e-5)) {
    warning("Parallel computing of szKendall among the cell pairs is not done properly. Please re-run the function and/or re-set the parallel backend.")
  }

  rownames(szkendall1.dist) <- NULL
  colnames(szkendall1.dist) <- NULL

  return(szkendall1.dist)
}


