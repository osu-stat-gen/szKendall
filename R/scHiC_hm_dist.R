#' Heatmap visualization of a dissimilarity matrix
#'
#' This function visualizes a dissimilarity matrix as a heatmap.
#'
#' @importFrom plsgenomics matrix.heatmap
#' @param data  An `n_single`x`n_single` symmetric dissimilarity matrix.
#' @param title The title given to the output plot. Default is empty ("").
#' @return A plot that visualizes the dissimilarity matrix as heatmap. Each element in the heatmap is scaled to be between 0 and 1 (by subtracting the minimum and dividing the range).
#' @examples
#' scHiC_hm_dist(euclid.diss, title="Euclid")
#' @export
scHiC_hm_dist <- function(data, title = ""){

  normmatr <- data
  normmatr[lower.tri(normmatr, diag=TRUE)] <- 0
  maxvalue <- max(data[upper.tri(data, diag=FALSE)])
  minvalue <- min(data[upper.tri(data, diag=FALSE)])
  normmatr[upper.tri(normmatr)] <- (data[upper.tri(data)] - minvalue)/(maxvalue - minvalue)
  normmatr <- normmatr+t(normmatr)

  col_pal = c('#30123BFF' ,'#321645FF' ,'#341A4EFF' ,'#351E57FF' ,
              '#362160FF' ,'#372569FF' ,'#392971FF' ,'#3A2D79FF' ,
              '#3B3082FF' ,'#3D3489FF' ,'#3E3790FF' ,'#3F3B97FF' ,
              '#3F3F9EFF' ,'#4142A5FF' ,'#4146ACFF' ,'#4249B2FF' ,
              '#424CB7FF' ,'#4450BEFF' ,'#4454C3FF' ,'#4457C8FF' ,
              '#455BCDFF' ,'#455ED3FF' ,'#4662D7FF' ,'#4665DCFF' ,
              '#4668DFFF' ,'#466BE3FF' ,'#476FE7FF' ,'#4772EAFF' ,
              '#4776EEFF' ,'#4778F0FF' ,'#477CF3FF' ,'#467FF5FF' ,
              '#4682F8FF' ,'#4686FAFF' ,'#4589FCFF' ,'#458CFDFF' ,
              '#448FFEFF' ,'#4392FEFF' ,'#4195FFFF' ,'#4099FFFF' ,
              '#3E9CFEFF' ,'#3C9FFEFF' ,'#3AA2FCFF' ,'#38A5FBFF' ,
              '#36A9F9FF' ,'#34ACF7FF' ,'#31AFF5FF' ,'#2FB2F4FF' ,
              '#2DB6F1FF' ,'#2AB9EFFF' ,'#28BCEBFF' ,'#26BFE8FF' ,
              '#24C2E5FF' ,'#22C5E2FF' ,'#20C7DFFF' ,'#1FCADCFF' ,
              '#1CCCD8FF' ,'#1BD0D5FF' ,'#1AD3D1FF' ,'#19D5CEFF' ,
              '#18D7CAFF' ,'#18D9C8FF' ,'#18DCC4FF' ,'#18DEC1FF' ,
              '#18E0BDFF' ,'#19E2BBFF' ,'#19E3B7FF' ,'#1CE6B4FF' ,
              '#1DE7B2FF' ,'#1FEAAEFF' ,'#21EBABFF' ,'#25ECA7FF' ,
              '#28EEA3FF' ,'#2BF09FFF' ,'#2EF19CFF' ,'#32F298FF' ,
              '#36F393FF' ,'#3BF58FFF' ,'#3FF68AFF' ,'#44F786FF' ,
              '#48F882FF' ,'#4DF97EFF' ,'#52FA7AFF' ,'#56FA75FF' ,
              '#5CFC70FF' ,'#61FC6CFF' ,'#66FD68FF' ,'#6BFD64FF' ,
              '#70FE60FF' ,'#75FE5CFF' ,'#7AFE58FF' ,'#7FFF54FF' ,
              '#84FF51FF' ,'#89FF4DFF' ,'#8DFF4AFF' ,'#91FF48FF' ,
              '#96FE44FF' ,'#9AFE41FF' ,'#9EFD3FFF' ,'#A1FD3DFF' ,
              '#A4FC3CFF' ,'#A8FC39FF' ,'#ABFB38FF' ,'#AFFA37FF' ,
              '#B2F936FF' ,'#B6F735FF' ,'#B9F635FF' ,'#BCF534FF' ,
              '#BFF434FF' ,'#C2F234FF' ,'#C6F034FF' ,'#C9EF34FF' ,
              '#CCEC34FF' ,'#CFEA34FF' ,'#D2E935FF' ,'#D5E635FF' ,
              '#D8E436FF' ,'#DBE236FF' ,'#DDE037FF' ,'#E0DE37FF' ,
              '#E3DB38FF' ,'#E5D938FF' ,'#E8D639FF' ,'#EAD439FF' ,
              '#ECD13AFF' ,'#EECF3AFF' ,'#F0CC3AFF' ,'#F2CA3AFF' ,
              '#F4C73AFF' ,'#F5C43AFF' ,'#F7C23AFF' ,'#F8BE39FF' ,
              '#F9BC39FF' ,'#FBB939FF' ,'#FBB737FF' ,'#FCB336FF' ,
              '#FCB036FF' ,'#FDAD34FF' ,'#FEA933FF' ,'#FEA732FF' ,
              '#FEA331FF' ,'#FE9F2FFF' ,'#FE9B2DFF' ,'#FE982CFF' ,
              '#FE942AFF' ,'#FE9129FF' ,'#FD8D27FF' ,'#FD8926FF' ,
              '#FC8524FF' ,'#FB8122FF' ,'#FB7D21FF' ,'#F97A1EFF' ,
              '#F9761DFF' ,'#F8721CFF' ,'#F76E1AFF' ,'#F56A18FF' ,
              '#F46617FF' ,'#F36215FF' ,'#F25F14FF' ,'#F05C12FF' ,
              '#EF5811FF' ,'#ED5410FF' ,'#EB510EFF' ,'#EA4E0DFF' ,
              '#E84B0CFF' ,'#E6480CFF' ,'#E4460AFF' ,'#E2430AFF' ,
              '#E14009FF' ,'#DE3E08FF' ,'#DC3B07FF' ,'#DA3907FF' ,
              '#D73606FF' ,'#D53405FF' ,'#D23105FF' ,'#D02F05FF' ,
              '#CD2C04FF' ,'#CA2A04FF' ,'#C82803FF' ,'#C42603FF' ,
              '#C22402FF' ,'#BE2102FF' ,'#BB2002FF' ,'#B81D02FF' ,
              '#B51B01FF' ,'#B21A01FF' ,'#AE1801FF' ,'#AA1701FF' ,
              '#A71401FF' ,'#A31301FF' ,'#A01101FF' ,'#9C0F01FF' ,
              '#980E01FF' ,'#940C01FF' ,'#900A01FF' ,'#8B0902FF' ,
              '#880802FF' ,'#830702FF' ,'#7F0502FF' ,'#7A0403FF')

  return(matrix.heatmap(normmatr, main = title, legend.width=0.8, legend.mar=3.1, col = col_pal))
}

