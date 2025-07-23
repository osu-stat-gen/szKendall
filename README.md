
<!-- README.md is generated from README.Rmd. Please edit that file -->

# szKendall

<!-- badges: start -->
<!-- badges: end -->

The goal of szKendall is to measure the dissimilarity between two single
cells based on their observed Hi-C data, accounting for the
discrepancies in structural zero positions.

## System Requirements

- R version: \>= 3.5
- Tested on: macOS 15.5, Windows 11
- Dependencies: None
- No special hardware required

## Installation

You can install the szKendall package from GitHub:

``` r
library(devtools)
devtools::install_github("https://github.com/osu-stat-gen/szKendall")
```

The installation would take a few minutes on a standard desktop
computer.

## Example

This is an example which shows you how to calculate the szKendall
dissimilarity between cells. The simulated data contain 150 single cells
(i.e., 150 columns). Therefore, the resulting szKendall matrix (which is
symmetric with the main-diagonal all being 0) is of dimension
150-by-150. With sequential calculation, the calculation time of one
szKendall dissimilarity in this example should be less than 2 minutes,
and half (or even less) with parallel computing.

``` r

library(szKendall)

# Example code
data("sim1.data")
data("true1.data")
# Calculate szKendall dissimilarity 
doParallel::registerDoParallel(cores = 4)  
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
