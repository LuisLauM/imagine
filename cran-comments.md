## Test environments
* ubuntu 16.04 (on travis-ci), R 3.3.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

## First submission

> Found the following (possibly) invalid URLs:
>  URL: http://cran.r-project.org/package=imagine
>    From: README.md
>    CRAN URL not in canonical form
>  URL: http://cran.rstudio.com/web/packages/imagine/index.html (moved to https://cran.rstudio.com/web/packages/imagine/index.html)
>    From: README.md
>    Status: 404
>    Message: Not Found
>    CRAN URL not in canonical form
>  A canonical CRAN URL starts with https://CRAN.R-project.org/

Done.


> Size of tarball: 15453551 bytes
>
> Pls fix: the CRAN Policy asks for no more than 5MB.

Done.


> Demos with missing or empty index information:
>  single_filters
> Demo index entries without corresponding demo:
>  median_filter
>  mean_filter

Done.


## Second submission

> Non-standard files/directories found at top level:
>  ‘README-unnamed-chunk-4-1.png’ ‘README.Rmd’ ‘README_cache’

Deleted.
