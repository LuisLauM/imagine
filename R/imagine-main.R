#' @title IMAGing engINE, Tools for application of image filters to data matrices
#'
#' @author Wencheng Lau-Medrano, \email{luis.laum@gmail.com}
#' @name image-package
#' @description Provides fast application of image filters to data matrices, using R and C++ algorithms.
#' @details This package uses C++ algorithms called 'engines'. More details are shown in the vignette.
#' @aliases imagine-package imagine
#' @docType package
#' @keywords image-matrix, image-filter
NULL

#' @title Make convolution calculations from numeric matrix
#'
#' @rdname convolutions
#'
#' @param X A \code{numeric matrix} object used for apply filters.
#' @param kernel A little matrix used as mask for each cell of \code{X}.
#' @param probs \code{numeric} vector of probabilities with values in [0,1].
#' @param times How many times do you want to apply the filter?
#' @param noNA \code{logical} indicating whether to make convolution only if all values within
#' kernel are not \code{NA}.
#'
#' @description This function takes a \code{matrix} object, and for each cell multiplies its neighborhood by
#' the \code{kernel}. Finally, it returns for each cell the mean of the kernel-weighted sum.
#'
#' @return \code{convolution2D} returns a \code{matrix} object with the same dimensions of \code{X}.
#'
#' @details
#' Convolution is a  mathematical operation which allows the multiplication of two arrays of numbers, in order
#' to produce an array of numbers of the same dimensionality. Functions use C++ algorithms. More details are shown in the vignette.
#'
#' @export
#'
#' @examples
#' # Generate example matrix
#' nRows <- 50
#' nCols <- 100
#'
#' myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)
#' kernel <- diag(3)
#'
#' # Make convolution
#' myOutput1 <- convolution2D(myMatrix, kernel)
#' myOutput2 <- convolutionQuantile(myMatrix, kernel, probs = 0.7)
#'
#' # Plot results
#' par(mfrow = c(2, 2))
#' image(myMatrix, zlim = c(0, 100))
#' image(myOutput1, zlim = c(0, 100))
#' image(myOutput2, zlim = c(0, 100))
convolution2D <- function(X, kernel, times = 1, noNA = FALSE){

  # Check and validation of arguments
  checkedArgs <- list(X = X, kernel = kernel, times = times, noNA = noNA)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = "convolution2D")

  # Apply filters
  output <- checkedArgs$X
  for(i in seq(checkedArgs$times)){
    gc(reset = TRUE)

    output <- with(checkedArgs, engine1(data = output, kernel = kernel, noNA = noNA))
  }

  return(output)
}

#' @rdname convolutions
#' @return \code{convolutionQuantile} uses the kernel but, for each cell, it returns the position
#' of quantile 'probs' (value between 0 and 1).
#' @export
convolutionQuantile <- function(X, kernel, probs, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(X = X, kernel = kernel, probs = probs, times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = "convolutionQuantile")

  checkedArgs$probs <- ceiling(checkedArgs$probs*prod(dim(kernel))) - 1

  # Apply filters
  output <- checkedArgs$X
  for(i in seq(checkedArgs$times)){
    gc(reset = TRUE)

    output <- with(checkedArgs, engine2(data = output, kernel = kernel, probs = probs))
  }

  return(output)
}


#' @rdname convolutions
#' @return \code{convolutionMedian} is a wrapper of \code{convolutionQuantile} with probs = 0.5.
#' @export
convolutionMedian <- function(X, kernel, times = 1){

  output <- convolutionQuantile(X = X, kernel = kernel, probs = 0.5, times = times)

  return(output)
}


#' @title Make a 2D filter calculations from numeric matrix
#'
#' @rdname basic2DFilter
#'
#' @param X A \code{numeric matrix} object used for apply filters.
#' @param radius Size of squared kernel to apply median.
#' @param probs \code{numeric} vector of probabilities with values in [0,1].
#' @param times How many times do you want to apply the filter?
#'
#' @description This functions take a \code{matrix} object, and for each cell multiplies its neighborhood by
#' the squared matrix of dimension \eqn{radius*radius}. Finally and according to the type of function, they
#' return for each cell the mean, median or the quantile of the weighted sum.
#'
#' @return \code{meanFilter} returns a \code{matrix} object with the same dimensions of \code{X}.
#'
#' @details Functions use C++ algorithms. More details are shown in the vignette.
#' @export
#'
#' @examples
#' # Generate example matrix
#' nRows <- 50
#' nCols <- 100
#'
#' myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)
#' radius <- 3
#'
#' # Make convolution
#' myOutput1 <- meanFilter(myMatrix, radius)
#' myOutput2 <- quantileFilter(myMatrix, radius, 0.1)
#' myOutput3 <- medianFilter(myMatrix, radius)
#'
#' # Plot results
#' image(myMatrix, zlim = c(0, 100))
#' image(myOutput1, zlim = c(0, 100))
#' image(myOutput2, zlim = c(0, 100))
#' image(myOutput3, zlim = c(0, 100))
meanFilter <- function(X, radius, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(X = X, radius = radius, times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = "meanFilter")

  # Apply filters
  output <- checkedArgs$X
  for(i in seq(checkedArgs$times)){
    gc(reset = TRUE)

    output <- with(checkedArgs, engine3(data = output, radius = radius))
  }

  return(output)
}


#' @rdname basic2DFilter
#' @return \code{quantileFilter} don't use a kernel but, for each cell, it returns the position
#' of quantile 'probs' (value between 0 and 1).
#' @export
quantileFilter <- function(X, radius, probs, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(X = X, radius = radius, probs = probs, times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = "quantileFilter")

  checkedArgs$probs <- ceiling(checkedArgs$probs*radius^2) - 1

  # Apply filters
  output <- checkedArgs$X
  for(i in seq(checkedArgs$times)){
    gc(reset = TRUE)

    output <- with(checkedArgs, engine4(data = output, radius = radius, probs = probs))
  }

  return(output)
}


#' @rdname basic2DFilter
#' @return \code{medianFilter} is a wrapper of \code{quantileFilter} with probs = 0.5.
#' @export
medianFilter <- function(X, radius, times = 1){

  output <- quantileFilter(X = X, radius = radius, probs = 0.5, times = times)

  return(output)
}

#' Performs Contextual Median Filter
#'
#' @param X A \code{numeric matrix} object used for apply filters.
#' @param inner_radius Size of the Inner squared kernel to apply median. Check Details.
#' @param outer_radius Size of the Outer squared kernel to apply median. Check Details.
#' @param probs \code{numeric} vector of probabilities with values in [0,1].
#' @param times How many times do you want to apply the filter?
#'
#' @description This function performs a generalization of the Contextual Median Filter propose by
#' Belkin & O'Reilly (2009). The default parameters reproduce the algorithm of the paper, but it allows
#' the users to modify values like the inner/outer matrices (check the paper for extra information) as
#' well as the quantile (as default, for median: \code{probs = 0.5}) and the number of applications (\code{times}).
#'
#' @references Belkin, I. M., & O'Reilly, J. E. (2009). An algorithm for oceanic front detection in chlorophyll
#' and SST satellite imagery. Journal of Marine Systems, 78(3), 319-326
#' (\url{http://dx.doi.org/10.1016/j.jmarsys.2008.11.018}).
#'
#' @return \code{contextualMF} returns a \code{matrix} object with the same dimensions of \code{X}.
#' @export
#'
#' @examples
#' # Generate example matrix
#' nRows <- 50
#' nCols <- 100
#'
#' myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)
#'
#' # Make convolution
#' myOutput <- contextualMF(X = myMatrix)
#'
#' # Plot results
#' image(myOutput, zlim = c(0, 100))
contextualMF <- function(X, inner_radius = 3, outer_radius = 5, probs = 0.5, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(X = X, inner_radius = inner_radius, outer_radius = outer_radius,
                      probs = probs, times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = "contextualMF")

  checkedArgs$probs <- ceiling(checkedArgs$probs*inner_radius^2) - 1

  # Apply filters
  output <- checkedArgs$X
  for(i in seq(checkedArgs$times)){
    gc(reset = TRUE)

    output <- with(checkedArgs, engine5(data = output, I_radius = inner_radius, O_radius = outer_radius, probs = probs))
  }

  return(output)
}


#' @title Data matrix to be used as example image.
#' @name wbImage
#' @description \code{matrix} object containig numeric data to plot a image. The photo was taken
#' by the author at 2016.
#' @aliases wbImage
#' @docType data
#' @usage wbImage
#' @format A \code{matrix} with dimnensions 1280x720.
NULL
