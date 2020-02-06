imagine
=======

[![packageversion](https://img.shields.io/badge/Package%20version-1.5.3-orange.svg?style=flat-square)](commits/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/imagine)](https://cran.r-project.org/package=imagine) [![CRAN_time_from_release](https://www.r-pkg.org/badges/ago/imagine)](https://cran.r-project.org/package=imagine) [![metacran downloads](https://cranlogs.r-pkg.org/badges/imagine)](https://cran.r-project.org/package=imagine) [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.1.0-6666ff.svg)](https://cran.r-project.org/)

**[IMAG]ing eng[INE]s, Tools for Application of Image Filters to Data Matrices**

Provides fast application of image filters to data matrices by using C++ algorithms called 'engines'. More details are shown in vignette.

Installation
------------

Get the development version from github:

``` r
# install.packages("devtools")
devtools::install_github("LuisLauM/imagine")
```

Or install the CRAN version

``` r
install.packages("imagine")
```

Input data
----------

For all functions, the main input data must be a `numeric matrix` object. Depending on each funtion, user must indicate some extra arguments for the filter.

Examples
--------

Next, we show the utility of `quantileFilter`, one of the six functions that `imagine` performs.

``` r
# Load imagine
library(imagine)

# Build an example matrix
n <- 1e3
origMatrix <- matrix(seq(n^2), nrow = n)

# Add some NAs
origMatrix_withNA <- origMatrix
origMatrix_withNA[sample(seq(n^2), 0.7*n^2, replace = FALSE)] <- NA

# Apply filter
newMatrix <- quantileFilter(X = origMatrix_withNA, radius = 3, x = 0.1, times = 1)

# Plot matrices for compare
cols <- colorRampPalette(c("green3", "red4"))(n)

par(mar = c(0, 2, 0, 0), mfrow = c(3, 1))

image(origMatrix, col = cols, axes = FALSE)
mtext(text = "Original", side = 2, line = 0.5, font = 2)

image(origMatrix_withNA, col = cols, axes = FALSE)
mtext(text = "Original with NAs", side = 2, line = 0.5, font = 2)

image(newMatrix, col = cols, axes = FALSE)
mtext(text = "Filtered", side = 2, line = 0.5, font = 2)
```
