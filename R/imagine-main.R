#' @title Imaging Engine, Tools for application of image filters to data matrices
#'
#' @author Wencheng Lau-Medrano, \email{luis.laum@gmail.com}
#' @name image-package
#' @description Provides fast application of image filters to data matrices, using R and C++ algorithms.
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
#' @return \code{convolution2D} returns a \code{matrix} object with the same dimensions of \code{dataMatrix}.
#'
#' @details Functions use the \code{engine1}, \code{engine2}, \code{engine3} C++ algorithms.
#' More details are shown in vignette.
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
  for(i in seq(checkedArgs$times)){
    gc(reset = TRUE)

    output <- with(checkedArgs, engine1(data = output, kernel = kernel))

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
    gc(reset = TRUE)

    output <- with(checkedArgs, engine2(data = output, kernel = kernel))

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

  checkedArgs$x <- ceiling(checkedArgs$x*prod(dim(kernel))) - 1

  # Apply filters
  output <- checkedArgs$dataMatrix
  for(i in seq(checkedArgs$times)){
    gc(reset = TRUE)

    maxValue <- max(an(output), na.rm = TRUE)*max(an(checkedArgs$kernel), na.rm = TRUE)

    output <- with(checkedArgs, engine3(data = output, kernel = kernel, x = x, maxValue = maxValue))
    output[output > (maxValue - 1)] <- NA

  }

  return(output)
}


#' @rdname convolutions
#' @return \code{convolutionMedian} is a wrapper of \code{convolutionQuantile} with x = 0.5.
#' @export
convolutionMedian <- function(dataMatrix, kernel, times = 1){

  output <- convolutionQuantile(dataMatrix = dataMatrix, kernel = kernel, x = 0.5, times = times)

  return(output)
}


#' @title Make a 2D filter calculations from numeric matrix
#'
#' @rdname basic2DFilter
#'
#' @param dataMatrix A \code{numeric matrix} object used for apply filters.
#' @param radius Size of squared kernel to apply median.
#' @param x \code{numeric} vector of probabilities with values in [0,1].
#' @param times How many times do you want to apply the filter?
#'
#' @description This functions take a \code{matrix} object, and for each cell multiplies its neighborhood by
#' the squared matrix of dimension \eqn{radius*radius}. Finally and according to the type of function, they
#' return for each cell the mean, median or the quantile of the weighted sum.
#'
#' @return \code{meanFilter} returns a \code{matrix} object with the same dimensions of \code{dataMatrix}.
#'
#' @details Functions use the \code{engine4} and \code{engine5} C++ algorithms. More details are shown in vignette.
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
meanFilter <- function(dataMatrix, radius, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(dataMatrix = dataMatrix, radius = radius, times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = as.character(match.call())[1])

  # Apply filters
  output <- checkedArgs$dataMatrix
  for(i in seq(checkedArgs$times)){
    gc(reset = TRUE)

    output <- with(checkedArgs, engine4(data = output, radius = radius))

  }

  return(output)
}


#' @rdname basic2DFilter
#' @return \code{quantileFilter} don't use a kernel but, for each cell, it returns the position
#' of quantile 'x' (value between 0 and 1).
#' @export
quantileFilter <- function(dataMatrix, radius, x, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(dataMatrix = dataMatrix, radius = radius, x = x, times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = as.character(match.call())[1])

  checkedArgs$x <- ceiling(checkedArgs$x*radius^2) - 1

  # Apply filters
  output <- checkedArgs$dataMatrix
  for(i in seq(checkedArgs$times)){

    maxValue <- max(an(output), na.rm = TRUE)*99

    output <- with(checkedArgs, engine5(data = output, radius = radius, x = x, maxValue = maxValue))
    output[output > (maxValue - 1)] <- NA

    gc(reset = TRUE)
  }

  return(output)
}


#' @rdname basic2DFilter
#' @return \code{medianFilter} is a wrapper of \code{quantileFilter} with x = 0.5.
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
