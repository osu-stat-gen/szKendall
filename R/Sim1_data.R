
#' Simulated count data from Sim1 setting
#'
#' This dataset contains simulated single-cell Hi-C data used to demonstrate the functionality
#' of the szKendall functions. It consists of 150 single cells of 3 subtypes (each with 50 cells).
#'
#' @format A simulated count matrix with 1830 rows and 150 columns:
#' \describe{
#'   Each row corresponds to one locus pair. The 61 loci result in 61*60/2 = 1830 locus pairs.
#'   Each column corresponds to one single cell. There are 3 subtypes, each contains 50 cells, resulting in 150 cells.
#' }
#'
#' @source These subtypes are derived based on the 3D structure of a segment in chromosome 19 (61 loci) from a K562 cell (accession: GSM2109974) obtained from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80006. The parameter setting is: alpha_0=2.8, gamma=0.5, eta=0.3 (or 0.5 or 0.8), tau1=0.8, tau2=0.6.
"sim1.data"


#' Expected count data from Sim1 setting
#'
#' This dataset contains expected single-cell Hi-C data used to demonstrate the functionality
#' of the szKendall functions. It consists of 150 single cells of 3 subtypes (each with 50 cells).
#'
#' @format The underlying expected count matrix with 1830 rows and 150 columns for sim1.data:
#' \describe{
#'   Each row corresponds to one locus pair. The 61 loci result in 61*60/2 = 1830 locus pairs.
#'   Each column corresponds to one single cell. There are 3 subtypes, each contains 50 cells, resulting in 150 cells.
#' }
#'
#' @source These subtypes are derived based on the 3D structure of a segment in chromosome 19 (61 loci) from a K562 cell (accession: GSM2109974) obtained from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80006. The parameter setting is: alpha_0=2.8, gamma=0.5, eta=0.3 (or 0.5 or 0.8), tau1=0.8, tau2=0.6.
"true1.data"

