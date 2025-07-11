
#' Observed human prefrontal cortex (L4/L5) whole-genome data
#'
#' This dataset contains the observed human prefrontal cortex (L4/L5) single-cell Hi-C data.
#'
#' @format A count matrix with 224456 rows and 311 columns:
#' \describe{
#'   Each row corresponds to one locus pair. The 22 autosomes have in total 224456 intra-chromosomal locus pairs. For the number of locus pairs per chromosome, see L4L5.loci.by.chr.
#'   Each column corresponds to one single cell. The first 131 columns are L4 cells, and the last 180 columns are L5 cells.
#' }
#'
#' @source The sn-m3C-seq dataset containing single-cell Hi-C profiles from human prefrontal cortex cells was downloaded from https://github.com/dixonlab/scm3C-seq
"L4L5.wholegenome"


#' Number of locus pairs per chromosome for the human prefrontal cortex (L4/L5) whole-genome data
#'
#' This vector lists the number of intra-chromosomal locus pairs for each autosome in the human prefrontal cortex (L4/L5) data.
#'
#' @format A numeric vector with 22 values corresponding to the number of intra-chromosomal locus pairs for each of the 22 autosomes.
#'
#' @source The sn-m3C-seq dataset containing single-cell Hi-C profiles from human prefrontal cortex cells was downloaded from https://github.com/dixonlab/scm3C-seq
"L4L5.loci.by.chr"

