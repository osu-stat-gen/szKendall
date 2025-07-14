
<!-- README.md is generated from README.Rmd. Please edit that file -->

# szKendall

<!-- badges: start -->
<!-- badges: end -->

The goal of szKendall is to measure the dissimilarity between two single
cells based on their observed Hi-C data, accounting for the
discrepancies in structural zero positions.

## Installation

You can install the szKendall package from GitHub:

``` r
library(devtools)
<<<<<<< HEAD
devtools::install_github("https://github.com/osu-stat-gen/szKendall")
=======
devtools::install_github("https://github.com/sl-lin/szKendall")
>>>>>>> a5325578bc228be39ed9650f92adea39c6de3c57
```

## Example

This is an example which shows you how to calculate the szKendall
dissimilarity between cells. The simulated data contain 150 single cells
(i.e., 150 columns). Therefore, the resulting szKendall matrix (which is
symmetric with the main-diagonal all being 0) is of dimension
150-by-150.

``` r

library(szKendall)

# Example code
data("sim1.data")
data("true1.data")
# Calculate szKendall dissimilarity 
doParallel::registerDoParallel(cores = 2)  
# Run foreach::registerDoSEQ() to override any default parallel backend and force sequential calculation 
# foreach::registerDoSEQ()
szK.sim1 <- szKendall.diss(sim1.data, true1.data)
szK1.sim1 <- szKendall1.diss(sim1.data, true1.data)
szK2.sim1 <- szKendall2.diss(sim1.data, true1.data)

print(dim(szK.sim1))
#> [1] 150 150
print(dim(szK1.sim1))
#> [1] 150 150
print(dim(szK2.sim1))
#> [1] 150 150
```
