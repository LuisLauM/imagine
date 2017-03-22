imagine
=======

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/imagine)](http://cran.r-project.org/package=imagine) [![](http://cranlogs.r-pkg.org/badges/imagine)](http://cran.rstudio.com/web/packages/imagine/index.html)

**Provides fast application of image filters to data matrices**

This package uses C++ algorithms called 'engines'. More details are shown in vignette.

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
newMatrix <- quantileFilter(dataMatrix = origMatrix_withNA, radius = 3, x = 0.1, times = 1)

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
