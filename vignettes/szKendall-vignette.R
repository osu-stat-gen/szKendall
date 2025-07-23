## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----warning=FALSE, message=FALSE, results='hide'-----------------------------
# If not already installed:
# install.packages("devtools")

# library(devtools)
# devtools::install_github("https://github.com/osu-stat-gen/szKendall")
library(szKendall)

# For visualization 
library(ggplot2)  
library(plsgenomics)
library(Rtsne)
library(umap)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(doParallel)
library(foreach)
registerDoParallel(cores = 4)  
# Run foreach::registerDoSEQ() to override any parallel backend and force sequential calculation 
# foreach::registerDoSEQ()

