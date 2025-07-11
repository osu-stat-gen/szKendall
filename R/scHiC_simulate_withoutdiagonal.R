
#' Simulate single-cell Hi-C data
#' Given the 3D coordinates of N loci, this function simulates the single-cell Hi-C data for n_single cells, along with the true structural zero positions
#' @importFrom stats runif rpois rbinom quantile
#' @param data An nx3 matrix representing the 3D coordinates of n loci.
#' @param alpha_0 Parameter that controls the sequencing depth of the simulated single-cell Hi-C data.
#' @param alpha_1 Parameter that reflects the biophysical law between the expected contact counts and the 3D Euclidean distances among the loci. Default value is -1.
#' @param beta_l Parameter that mimics the effect of fragment length in generating single-cell Hi-C data. Default value is 0.9.
#' @param beta_g Parameter that mimics the effect of GC content in generating single-cell Hi-C data. Default value is 0.9.
#' @param beta_m Parameter that mimics the effect of mappability in generating single-cell Hi-C data. Default value is 0.9.
#' @param gamma Parameter that defines the threshold for expected contact counts; locus pairs with expected counts below this percentile could become structural zero candidates.
#' @param tau1 Parameter that sets the probability that a locus pair with expected count below gamma percentile turning into a structural zero candidate.
#' @param eta Parameter that decides the proportion of common structural zeros among the structural zero candidates.
#' @param tau2 Parameter that decides the proportion of cell-specific structural zeros among the remaining (1 - eta) proportion of the structural zero candidates.
#' @param n_single The number of single cells to be generated.
#' @return A list with 2 elements:
#' \describe{
#'   \item{singledat}{The simulated single-cell Hi-C contact counts. The dimension is N x n_single, where N = n*(n-1)/2 is the number of unique locus pairs given n loci. }
#'   \item{truecount}{The expected ("true") single-cell Hi-C contact counts. Same dimension as singledat. }
#' }
#' @examples
#' n_single <- 50
#' data("cooord3D_T1")
#' set.seed(1234)
#' ST1 <- scHiC_simulate_withoutdiagonal(
#'   data = coord3D_T1,
#'   alpha_0 = 1.8,
#'   alpha_1 = -1,
#'   beta_l = 0.9,
#'   beta_g = 0.9,
#'   beta_m = 0.9,
#'   gamma = 0.6,
#'   tau1 = 0.7,
#'   eta = 0.1,
#'   tau2=0.6,
#'   n_single=n_single
#' )
#' dim(ST1$singledat)
#' dim(ST1$truecount)
#' @export
scHiC_simulate_withoutdiagonal <- function(data, alpha_0, alpha_1=-1, beta_l=0.9, beta_g=0.9, beta_m=0.9, gamma, tau1, eta, tau2, n_single){
  nposi <- dim(data)[1]
  position <- cbind(data, runif(nposi, min = 0.2, max = 0.3),
                    runif(nposi, min = 0.4, max = 0.5), runif(nposi, min = 0.9, max = 1))

  distance <- matrix(0, nrow = nposi, ncol = nposi)
  for (i in 1:(nposi - 1)) {
    for (j in (i + 1):nposi) {
      distance[i, j] = sqrt(sum((position[i, 1:3] - position[j, 1:3])^2))
    }
  }

  lambda <- matrix(0, nrow = nposi, ncol = nposi)
  for (i in 1:(nposi - 1)) {
    for (j in (i + 1):nposi) {
      l = position[i, 4] * position[j, 4]
      g = position[i, 5] * position[j, 5]
      m = position[i, 6] * position[j, 6]
      lambda[i, j] <- exp(alpha_0 + alpha_1 * log(distance[i, j]) + beta_l * log(l) + beta_g * log(g) + beta_m * log(m))
    }
  }

  seqdepth <- sum(lambda[upper.tri(lambda)])
  thresh <- quantile(lambda[upper.tri(lambda)], gamma)    # here should be "gamma", the original code was "eta"
  posi0 <- which(lambda < thresh & upper.tri(lambda), TRUE)
  true0 <- posi0[sample(nrow(posi0), size = tau1 * nrow(posi0), replace = FALSE), ]

  rpoi <- function(x) {
    r <- rpois(1, x)
    return(r)
  }

  matrow <- function(x, y) {
    r <- x + (y - 1) * (y - 2)/2
    return(r)
  }

  downsamplemat <- function(mat, samplerate = 0.5) {
    new <- matrix(0, nrow(mat), ncol(mat))
    for (i in 1:nrow(mat)) {
      for (j in 1:ncol(mat)) {
        new[i, j] <- sum(runif(mat[i, j], 0, 1) < samplerate)
      }
    }
    return(new)
  }

  subsampling <- function(v, eta) {
    out = rep(0, length(v))
    v1 = v
    v1[v == 0] = 1
    v2 = cumsum(v1)
    v3 = cumsum(v1) - v + 1
    total = sum(v1)
    index = sample(1:total, total * eta, replace = F)
    for (i in 1:length(index)) {
      intersect.row = (v3 <= index[i] & v2 >= index[i])
      out[intersect.row] = out[intersect.row] + 1
    }
    out[v == 0] = 0
    return(out)
  }

  l <- dim(true0)[1]
  if (eta == 0) {
    underline <- lambda
    underline <- underline[upper.tri(underline, diag = FALSE)]
    single <- matrix(rep(underline, each = n_single), ncol = n_single, byrow = TRUE)
    random0 <- NULL
    my_list <- list(single, random0)
    names(my_list) <- c("single", "random0")
    return(my_list)
  }
  else {
    true0rows <- NULL
    for (i in 1:(dim(true0)[1])) {
      true0rows <- c(true0rows, matrow(true0[i, 1], true0[i, 2]))
    }
    random0 <- true0rows[sample(1:l, floor(l * (1-eta)))]
    underline <- lambda
    underline <- underline[upper.tri(underline, diag = FALSE)]
    randlam <- underline[random0]
    underline[true0rows] <- 0
    truecount <- NULL
    for (j in 1:n_single) {
      s <- underline
      s[random0] <- randlam * rbinom(length(random0), 1, 1-tau2)
      truecount <- cbind(truecount, s)
    }
    singledat = apply(truecount, c(1, 2), rpoi)
    # my_list <- list(truecount, random0, singledat)
    # names(my_list) <- c("truecount", "random0", "singledat")
    my_list <- list(singledat, truecount)
    names(my_list) <- c("singledat", "truecount")
    return(my_list)
  }
}


