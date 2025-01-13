# imagine 2.1.2

* Minor correction on `is_extreme` C++ function.
* Minor corrections on `contextualMF` documentation.

# imagine 2.1.1

* Adding the `na_only` argument in several functions allowing to apply the filters only in the cells where there is `NA` and replacing with the original value in the rest.

# imagine 2.0.0

* Added a `NEWS.md` file to track changes to the package.
* Important corrections, improvements and changes in engines 2, 4 and 5, so as in functions `convolutionQuantile`, `convolutionMedian`, `quantileFilter` and `contextualMF`.
* `times` argument in `contextualMF` is not longer available.
* `na` argument is removed from previous functions: only `NA` will be considered as a `NA`.

# imagine 2.1.0

* Adding new function (`agenbagFilters`) that performs two methods for gradient calculation.
* Some minor improvements in documentation, vignettes and code.
