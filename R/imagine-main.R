# imagine: Imaging Engine, Tools for Fast application of image filters ----
#'
#' @title Imaging Engine, Tools for application of image filters to data matrices
#'
#' @author Wencheng Lau-Medrano, \email{luis.laum@gmail.com}
#' @name image-package
#' @description Package for faster application of image filters to data matrices, using R and C++ algorithms.
#' @details This package uses C++ algorithms called 'engines'. More details are shown in vignette.
#' @aliases imagine-package imagine
#' @docType package
#' @keywords image-matrix, image-filter
NULL

#' @title Make convolution calculations from numeric matrix
#'
#' @rdname convolutions
#'
#' @param dataMatrix A \code{numeric matrix} object used for apply filters.
#' @param kernel A little matrix used as mask for each cell of \code{dataMatrix}.
#' @param x \code{numeric} vector of probabilities with values in [0,1].
#' @param times How many times do you want to apply the filter?
#'
#' @description This function takes a \code{matrix} object, and for each cell multiplies its neighborhood by
#' the \code{kernel}. Finally, it returns for each cell the mean of the kernel-weighted sum.
#'
#' @return It returns a \code{matrix} object with the same dimensions of \code{dataMatrix}.
#'
#' @details This function uses the \code{engine1} C++ algorithm. More details are shown in vignette.
#' @export
#'
#' @examples
#' # Generate example matrix
#' nRows <- 10
#' nCols <- 20
#'
#' myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)
#' kernel <- diag(3)
#'
#' # Make convolution
#' myOutput1 <- convolution2D(myMatrix, kernel)
#' myOutput2 <- convolutionMean(myMatrix, kernel)
#' myOutput3 <- convolutionQuantile(myMatrix, kernel, x = 0.7)
#'
#' # Plot results
#' par(mfrow = c(2, 2))
#' image(myMatrix, zlim = c(0, 100))
#' image(myOutput1, zlim = c(0, 100))
#' image(myOutput2, zlim = c(0, 100))
#' image(myOutput3, zlim = c(0, 100))
convolution2D <- function(dataMatrix, kernel, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(dataMatrix = dataMatrix, kernel = kernel, times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = as.character(match.call())[1])

  # Apply filters
  output <- checkedArgs$dataMatrix
  for(i in seq(times)){
    output <- with(checkedArgs,
                   convolution2D_internal(dataMatrix = output, kernel = kernel))
    gc()
  }

  return(output)
}

#' @rdname convolutions
#' @return \code{convolutionMean} uses the kernel but, for each cell, it returns the average.
#' @export
convolutionMean <- function(dataMatrix, kernel, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(dataMatrix = dataMatrix, kernel = kernel, times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = as.character(match.call())[1])

  # Apply filters
  output <- checkedArgs$dataMatrix
  for(i in seq(times)){
    output <- with(checkedArgs,
                   convolutionMean_internal(dataMatrix = output, kernel = kernel))
    gc()
  }

  return(output)
}

#' @rdname convolutions
#' @return \code{convolutionQuantile} uses the kernel but, for each cell, it returns the position
#' of quantile 'x' (value between 0 and 1).
#' @export
convolutionQuantile <- function(dataMatrix, kernel, x, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(dataMatrix = dataMatrix, kernel = kernel, x = x, times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = as.character(match.call())[1])

  # Apply filters
  output <- checkedArgs$dataMatrix
  for(i in seq(times)){
    output <- with(checkedArgs,
                   convolutionQuantile_internal(dataMatrix = output, kernel = kernel, x = x))
    gc()
  }

  return(output)
}


#' @rdname convolutions
#' @return \code{convolutionMedian} is a wrapper of \code{convolutionQuantile} with x = 0.5
#' @export
convolutionMedian <- function(dataMatrix, kernel, times = 1){

  output <- convolutionQuantile(dataMatrix = dataMatrix, kernel = kernel, x = 0.5, times = times)

  return(output)
}


#' @title Make a 2D median filter calculation from numeric matrix
#'
#' @rdname medianFilter2D
#'
#' @param dataMatrix A \code{numeric matrix} object used for apply filters.
#' @param radius Size of squared kernel to apply median.
#' @param x \code{numeric} vector of probabilities with values in [0,1].
#' @param times How many times do you want to apply the filter?
#'
#' @description This function takes a \code{matrix} object, and for each cell multiplies its neighborhood by
#' the squared matrix of dimension \eqn{radius*radius}. Finally, it returns for each cell the median of the
#' weighted sum.
#'
#' @return It returns a \code{matrix} object with the same dimensions of \code{dataMatrix}.
#'
#' @details This function uses the \code{engine2} C++ algorithm. More details are shown in vignette.
#' @export
#'
#' @examples
#' # Generate example matrix
#' nRows <- 10
#' nCols <- 20
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
meanFilter <- function(dataMatrix, radius, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(dataMatrix = dataMatrix, radius = radius, times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = as.character(match.call())[1])

  # Apply filters
  output <- checkedArgs$dataMatrix
  for(i in seq(times)){
    output <- with(checkedArgs,
                   meanFilter_internal(dataMatrix = output, radius = radius))
    gc()
  }

  return(output)
}


#' @rdname medianFilter2D
#' @return \code{quantileFilter} uses the kernel but, for each cell, it returns the position
#' of quantile 'x' (value between 0 and 1).
#' @export
quantileFilter <- function(dataMatrix, radius, x, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(dataMatrix = dataMatrix, radius = radius, x = x, times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = as.character(match.call())[1])

  # Apply filters
  output <- checkedArgs$dataMatrix
  for(i in seq(times)){
    output <- with(checkedArgs,
                   quantileFilter_internal(dataMatrix = output, radius = radius, x = x))
    gc()
  }

  return(output)
}


#' @rdname medianFilter2D
#' @return \code{medianFilter} is a wrapper of \code{quantileFilter} with x = 0.5
#' @export
medianFilter <- function(dataMatrix, radius, times = 1){

  output <- quantileFilter(dataMatrix = dataMatrix, radius = radius, x = 0.5, times = times)

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
