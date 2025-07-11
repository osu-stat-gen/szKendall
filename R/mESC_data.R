
#' Observed mouse Embryonic Stem Cell (mESC) chromosome 1 data
#'
#' This dataset contains the observed mESC single-cell Hi-C data for 42 cells. The cells are at three different cell-cycle phases (G1, early-S, late-S/G2) --- see mESC.cell.type.
#'
#' @format A count matrix with 19503 rows and 42 columns:
#' \describe{
#'   Each row corresponds to one locus pair. The 198 loci result in 198*197/2 = 19503 locus pairs.
#'   Each column corresponds to one single cell. There are 42 cells: 30 cells in G1 phase, 4 cells in early-S phase, and 8 cells in late-S/G2 phase.
#' }
#'
#' @source The observed data was processed by scHi-CSim (Fan et al., 2023) (https://github.com/zhanglabtools/scHi-CSim/tree/main/data, with GEO accession: GSE94489).
"mESC.sim.data"


#' Expected mouse Embryonic Stem Cell (mESC) chromosome 1 data
#'
#' This dataset contains the expected count ("true") mESC single-cell Hi-C data for 42 cells. The cells are at three different cell-cycle phases (G1, early-S, late-S/G2) --- see mESC.cell.type.
#'
#' @format A count matrix with 19503 rows and 42 columns:
#' \describe{
#'   Each row corresponds to one locus pair. The 198 loci result in 198*197/2 = 19503 locus pairs.
#'   Each column corresponds to one single cell. There are 42 cells: 30 cells in G1 phase, 4 cells in early-S phase, and 8 cells in late-S/G2 phase.
#' }
#'
#' @source The expected count data was obtained by smoothing the observed count data. See szKendall manuscript for details.
"mESC.true.data"


#' Phases in cell cycle for the mouse Embryonic Stem Cell (mESC) chromosome 1 data
#'
#' This vector contains the cell-cycle phase each mESC single cell resides.
#'
#' @format A factor vector with 3 levels corresponding to 3 different cell-cycle phases (G1, early-S, and late-S/G2).
#'
#' @source The cell-cycle phases for the mESC single cells were extracted from the "group" column of "GSE94489_2i_diploids_features_table.txt" downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94489
"mESC.cell.type"


#' Simulated true structural zero (SZ) positions for the mouse Embryonic Stem Cell (mESC) chromosome 1 data
#'
#' This list contains the simulated true SZ position indexes for the 42 mESC cells.
#'
#' @format A list of 42 vectors:
#' \describe{
#'   Each vector represents the locus pair indexes of the simulated true SZ positions in one mESC cell.
#' }
#'
#' @source The cell-cycle phases for the mESC single cells are extracted from the "group" column of "GSE94489_2i_diploids_features_table.txt" downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94489
"mESC.SZ.index.list"

