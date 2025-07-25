% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/szKendall1_function.R
\name{szKendall1.diss}
\alias{szKendall1.diss}
\title{Calculate szKendall1 dissimilarity}
\usage{
szKendall1.diss(sim.data, true.data)
}
\arguments{
\item{sim.data}{A simulated or observed single-cell Hi-C matrix, where the rows represent locus pairs and columns represent cells.}

\item{true.data}{The expected or imputed single-cell Hi-C matrix, which has the same dimension as sim.data.}
}
\value{
A square szKendall1 dissimilarity matrix, where the dimension is the number of single cells.
}
\description{
This function computes the szKendall1 dissimilarity matrix given an "observed" locus-pair by cell single-cell Hi-C matrix and the "true" expected contact count matrix (where only structural zero positions have the value zero).
}
\examples{
foreach::registerDoSEQ()
szKendall1.diss(sim1.data, true1.data)
}
