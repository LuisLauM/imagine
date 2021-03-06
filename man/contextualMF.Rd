% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imagine-main.R
\name{contextualMF}
\alias{contextualMF}
\title{Performs Contextual Median Filter}
\usage{
contextualMF(X, times = 1, na = NA)
}
\arguments{
\item{X}{A \code{numeric matrix} object used for apply filters.}

\item{times}{How many times do you want to apply the filter?}

\item{na}{\code{NA} as default. But, if specified, it must be an integer value higher
than the maximum of \code{X}.}
}
\value{
\code{contextualMF} returns a \code{matrix} object with the same dimensions
of \code{X}.
}
\description{
This function performs the Contextual Median Filter proposed by Belkin &
O'Reilly (2009), based on the pseudo-code written on the paper.
}
\details{
The users can change the number of recursive applications by using \code{times} argument
(\code{times = 1} as default).
}
\examples{
# Generate example matrix
nRows <- 50
nCols <- 100

myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)

# Make convolution
myOutput <- contextualMF(X = myMatrix)

# Plot results
image(myOutput, zlim = c(0, 100))
}
\references{
Belkin, I. M., & O'Reilly, J. E. (2009). An algorithm for oceanic front
detection in chlorophyll and SST satellite imagery. Journal of Marine Systems, 78(3),
319-326 (\url{http://dx.doi.org/10.1016/j.jmarsys.2008.11.018}).
}
