% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scHiC_simulate_withoutdiagonal.R
\name{scHiC_simulate_withoutdiagonal}
\alias{scHiC_simulate_withoutdiagonal}
\title{Simulate single-cell Hi-C data
Given the 3D coordinates of n loci, this function simulates the single-cell Hi-C data for \code{n_single} cells, along with the true structural zero positions.}
\usage{
scHiC_simulate_withoutdiagonal(
  data,
  alpha_0,
  alpha_1 = -1,
  beta_l = 0.9,
  beta_g = 0.9,
  beta_m = 0.9,
  gamma,
  tau1,
  eta,
  tau2,
  n_single
)
}
\arguments{
\item{data}{An \code{n}x3 matrix representing the 3D coordinates of \code{n} loci.}

\item{alpha_0}{Parameter that controls the sequencing depth of the simulated single-cell Hi-C data.}

\item{alpha_1}{Parameter that reflects the biophysical law between the expected contact counts and the 3D Euclidean distances among the loci. Default value is -1.}

\item{beta_l}{Parameter that mimics the effect of fragment length in generating single-cell Hi-C data. Default value is 0.9.}

\item{beta_g}{Parameter that mimics the effect of GC content in generating single-cell Hi-C data. Default value is 0.9.}

\item{beta_m}{Parameter that mimics the effect of mappability in generating single-cell Hi-C data. Default value is 0.9.}

\item{gamma}{Parameter that defines the threshold for expected contact counts; locus pairs with expected counts below this percentile could become structural zero candidates.}

\item{tau1}{Parameter that sets the probability that a locus pair with expected count below \code{gamma} percentile turning into a structural zero candidate.}

\item{eta}{Parameter that decides the proportion of common structural zeros among the structural zero candidates.}

\item{tau2}{Parameter that decides the proportion of cell-specific structural zeros among the remaining \code{1 - eta} proportion of the structural zero candidates.}

\item{n_single}{The number of single cells to be generated.}
}
\value{
A list with 2 elements:
\describe{
\item{singledat}{The simulated single-cell Hi-C contact counts. The dimension is \code{N}x\code{n_single}, where \code{N = n*(n-1)/2} is the number of unique locus pairs given \code{n} loci. }
\item{truecount}{The expected ("true") single-cell Hi-C contact counts. Same dimension as singledat. }
}
}
\description{
Simulate single-cell Hi-C data
Given the 3D coordinates of n loci, this function simulates the single-cell Hi-C data for \code{n_single} cells, along with the true structural zero positions.
}
\examples{
n_single <- 50
data("coord3D_T1")
set.seed(1234)
ST1 <- scHiC_simulate_withoutdiagonal(
  data = coord3D_T1,
  alpha_0 = 1.8,
  alpha_1 = -1,
  beta_l = 0.9,
  beta_g = 0.9,
  beta_m = 0.9,
  gamma = 0.6,
  tau1 = 0.7,
  eta = 0.1,
  tau2=0.6,
  n_single=n_single
)
dim(ST1$singledat)
dim(ST1$truecount)
}
