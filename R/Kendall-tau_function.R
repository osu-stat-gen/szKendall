#' kendall
#' Kendall's tau matrix for loci pair (i,j) in single cell 1 and loci pair (u,v) in single cell 2 for all (i,j) and (u,v)
#' @param mat_count A \eqn{(\#LP) \times (\#cells)} matrix.  
#' @return A square Kendall dissimilarity matrix, where the dimension is the number of single cells.
#' @importFrom pcaPP cor.fk
#' @export
kendall.diss <- function(mat_count){
  # Register a parallel backend using the following lines if it is not done first:
  if (!foreach::getDoParRegistered()) {
    numCores <- max(1, parallel::detectCores() - 1)
    doParallel::registerDoParallel(cores = numCores)
    on.exit(doParallel::stopImplicitCluster())  # Clean up
  }

  n.cells <- ncol(mat_count)
  n1 <- ceiling(sqrt(2*nrow(mat_count)))

  kendall.dist <- matrix(0, nrow=n.cells, ncol=n.cells)

  fail.index <- c()
  for(i in 1:(n.cells-1)){
    for(j in (i+1):(n.cells)){
      if(kendall.dist[i,j]<1e-5){
        fail.index <- rbind(fail.index, c(i,j))
      }
    }
  }
  scoreK <- foreach(para = 1:nrow(fail.index), .combine = 'rbind',.noexport = c("szkendall1_cpp_SZpart")) %dopar% {
                        i <- fail.index[para, 1]; j <- fail.index[para, 2]
                        val <- kendall_distance2(mat_count[, i],
                                                 mat_count[, j])
                        c(i, j, val)
                        }

  
  for (m in 1:nrow(scoreK)) {
      kendall.dist[scoreK[m, 1], scoreK[m, 2]] <- scoreK[m, 3]
  }
  kendall.dist <- (kendall.dist + t(kendall.dist))

  rownames(kendall.dist) <- NULL
  colnames(kendall.dist) <- NULL

  return(kendall.dist)
}




# Traditional kendall's tau for loci pair (i,j) in single-cell 1 and loci pair (u,v) in single-cell 2 for all (i,j) and (u,v)
kendall <- function(Y1, Y2){
  if(length(Y1) != length(Y2)){
    stop("Error: The two single cell contact count vectors do not have the same length!")
  }else{
    len <- length(Y1)  # len = n*(n-1)/2 
    n <- ceiling(sqrt(2*len)) # number of bins, at least 2 
    K <- matrix(0, nrow=len, ncol=len)
    
    row.idx <- c()
    col.idx <- c()
    for(i in 2:n){
      row.idx <- c(row.idx, 1:(i-1))
      col.idx <- c(col.idx, rep(i, i-1))
    }
    
    for(r in 1:(len-1)){
      for(s in (r+1):len){
        K[r, s] <- (as.numeric((Y1[r]-Y1[s])*(Y2[r]-Y2[s])<0) + 0.5*as.numeric((Y1[r]-Y1[s])*(Y2[r]-Y2[s])==0))   
      }
    }
    
    return(K)
  } 
}




#------------------------------------------------------------------------------------------------------------------#

# A faster function to calculate the Kendall's tau distance 

# Given a vector x, this function calculate the sum of t_i*(t_i-1)/2, 
# where t_i is the number of tied values in the i^th group of ties in x
tie_x_count <- function(x){
  
  x <- x[order(x)]
  tieCount <- 0
  m1 <- 0
  
  if(length(x) > 1){
    for(k in 2:length(x)){
      if(x[k-1] == x[k]){
        tieCount <- tieCount+1
      }else if(tieCount > 0){
        m1 <- m1+tieCount*(tieCount+1)/2
        tieCount <- 0
      }
    }
  }
  
  if(tieCount > 0){
    m1 <- m1+tieCount*(tieCount+1)/2
  }
  
  return(m1)
}


# Given two vectors x and y, this function calculate the sum of t_i*(t_i-1)/2, 
# where t_i is the number of tied values in the i^th group of ties in (x,y)
tie_xy_count <- function(x, y){
  
  ord <- order(x,y)
  x <- x[ord]
  y <- y[ord]
  tieCount <- 0
  m2 <- 0
  
  for(k in 2:length(x)){
    if(x[k-1] == x[k] & y[k-1] == y[k]){
      tieCount <- tieCount+1
    }else if(tieCount > 0){
      m2 <- m2+tieCount*(tieCount+1)/2
      tieCount <- 0
    }
  }
  
  if(tieCount > 0){
    m2 <- m2+tieCount*(tieCount+1)/2
  }
  
  return(m2)
}



# Use Kendall's tau corrlation (K_c, tau_b) to calculate Kendall's tau distance (K_d) -- sum of all the 0,1/2,1 values in the n*n K matrix
kendall_distance <- function(Y1, Y2){
  if(length(Y1) != length(Y2)){
    stop("Error: The two single cell contact count vectors do not have the same length!")
  }else{
    len <- length(Y1)  # len = n*(n-1)/2 

    if(len <= 1){
      K_c <- 0
    } else {
      K_c <- pcaPP::cor.fk(Y1, Y2)
    }
    m1x <- tie_x_count(Y1)
    m1y <- tie_x_count(Y2)
    m2 <- tie_xy_count(Y1, Y2)
    
    K_d <- 1/2*(len*(len-1)/2 - m1x - m1y + m2 - K_c*sqrt(len*(len-1)/2-m1x)*sqrt(len*(len-1)/2-m1y))
    
    return(K_d+0.5*(m1x+m1y-m2))
  } 
}


kendall_distance2 <- function(Y1, Y2){
  if(length(Y1) != length(Y2)){
    stop("Error: The two single cell contact count vectors do not have the same length!")
  }else{
    len <- length(Y1)  # len = n*(n-1)/2 
    
    if(len <= 1){
      K_c <- 0
    } else {
      K_c <- pcaPP::cor.fk(Y1, Y2)
    }
    m1x <- tie_x_count(Y1)
    m1y <- tie_x_count(Y2)
    
    total <- len*(len-1)/4-0.5*K_c*sqrt(len*(len-1)/2-m1x)*sqrt(len*(len-1)/2-m1y)
    
    return(total)
  } 
}


# scaled to be between 0 and 1! 
kendall_distance3 <- function(Y1, Y2){
  if(length(Y1) != length(Y2)){
    stop("Error: The two single cell contact count vectors do not have the same length!")
  }else{
    len <- length(Y1)  # len = n*(n-1)/2 
    
    if(len <= 1){
      K_c <- 0
    } else {
      K_c <- pcaPP::cor.fk(Y1, Y2)
    }
    m1x <- tie_x_count(Y1)
    m1y <- tie_x_count(Y2)
    
    total <- len*(len-1)/4-0.5*K_c*sqrt(len*(len-1)/2-m1x)*sqrt(len*(len-1)/2-m1y)
    
    return(total*2/(len*(len-1)))
  } 
}


# # Test # 
# 
# tic <- Sys.time()
# kendall_distance(c(5,3,0), c(0,2,1))
# Sys.time() - tic 
# 
# tic <- Sys.time()
# kendall_distance2(c(5,3,0), c(0,2,1))
# Sys.time() - tic 
# 
# tic <- Sys.time()
# sum(kendall(c(5,3,0), c(0,2,1)))
# Sys.time() - tic 
# 
# tic <- Sys.time()
# kendall_distance3(c(5,3,0), c(0,2,1))
# Sys.time() - tic 
# 
# 
# tic <- Sys.time()
# kendall_distance(c(rep(2,500), rep(1,500)), c(1:1000))
# Sys.time() - tic 
# 
# tic <- Sys.time()
# kendall_distance2(c(rep(2,500), rep(1,500)), c(1:1000))
# Sys.time() - tic 
# 
# tic <- Sys.time()
# sum(kendall(c(rep(2,500), rep(1,500)), c(1:1000)))
# Sys.time() - tic 
# 
# tic <- Sys.time()
# kendall_distance3(c(rep(2,500), rep(1,500)), c(1:1000))
# Sys.time() - tic 
# 
# 
# tic <- Sys.time()
# kendall_distance(c(1000:1), c(1:1000))
# Sys.time() - tic 
# 
# tic <- Sys.time()
# kendall_distance2(c(1000:1), c(1:1000))
# Sys.time() - tic 
# 
# tic <- Sys.time()
# sum(kendall(c(1000:1), c(1:1000)))
# Sys.time() - tic 
# 
# tic <- Sys.time()
# kendall_distance3(c(1000:1), c(1:1000))
# Sys.time() - tic 
# 
# 
# a <- sample(1:1830)
# b <- sample(1:1830)
# 
# tic <- Sys.time()
# kendall_distance(a, b)
# Sys.time() - tic 
# 
# tic <- Sys.time()
# kendall_distance2(a, b)
# Sys.time() - tic 
# 
# tic <- Sys.time()
# sum(kendall(a, b))
# Sys.time() - tic 
# 
# tic <- Sys.time()
# kendall_distance3(a, b)
# Sys.time() - tic 


