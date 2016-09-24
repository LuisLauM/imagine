# imagine: Imaging Engine, Tools for Fast application of image filters ----
#'
#' @title Imaging Engine, Tools for application of image filters to data matrices
#'
#' @author Wencheng Lau-Medrano, \email{luis.laum@gmail.com}
#' @name image-package
#' @description Package for faster application of image filters to data matrices, using R and C++ algorithms.
#' @aliases imagine-package imagine
#' @docType package
#' @keywords image-matrix, image-filter
NULL

#' @title Make convolution calculation from numeric matrix
#'
#' @param dataMatrix A \code{numeric matrix} object used for apply filters.
#' @param kernel A little matrix used as mask for each cell of \code{dataMatrix}.
#' @param times How many times do you want to apply the filter?
#'
#' @description This function takes a \code{matrix} object, and for each cell multiplies its neighborhood by
#' the \code{kernel}. Finally, it returns for each cell the mean of the kernel-weighted sum.
#'
#' @return It returns a \code{matrix} object with the same dimensions of \code{dataMatrix}.
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
#' myOutput <- convolution2D(myMatrix, kernel)
#'
#' # Plot results
#' image(myMatrix, zlim = c(0, 100))
#' image(myOutput, zlim = c(0, 100))
convolution2D <- function(dataMatrix, kernel, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(dataMatrix, kernel, times)
  checkedArgs <- checkArgs(checkedArgs, type = "convolution2D")

  # Apply filters
  output <- dataMatrix
  for(i in seq(times)){
    output <- with(checkedArgs,
                   convolution2D_internal(dataMatrix = output, kernel = kernel))
  }

  return(output)
}

#' @title Make a 2D median filter calculation from numeric matrix
#'
#' @param dataMatrix A \code{numeric matrix} object used for apply filters.
#' @param radius Size of squared kernel to apply median.
#' @param times How many times do you want to apply the filter?
#'
#' @description This function takes a \code{matrix} object, and for each cell multiplies its neighborhood by
#' the squared matrix of dimension \eqn{radius*radius}. Finally, it returns for each cell the median of the
#' weighted sum.
#'
#' @return It returns a \code{matrix} object with the same dimensions of \code{dataMatrix}.
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
#' myOutput <- medianFilter2D(myMatrix, radius)
#'
#' # Plot results
#' image(myMatrix, zlim = c(0, 100))
#' image(myOutput, zlim = c(0, 100))
medianFilter2D <- function(dataMatrix, radius, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(dataMatrix, radius, times)
  checkedArgs <- checkArgs(checkedArgs, type = "medianFilter2D")

  # Apply filters
  output <- dataMatrix
  for(i in seq(times)){
    output <- with(checkedArgs,
                   medianFilter2D_internal(dataMatrix = output, radius = radius))
  }

  return(output)
}

#' @title Make 2D mean filter calculation from numeric matrix
#'
#' @param dataMatrix A \code{numeric matrix} object used for apply filters.
#' @param radius Size of squared kernel to apply mean.
#' @param times How many times do you want to apply the filter?
#'
#' @description This function takes a \code{matrix} object, and for each cell multiplies its neighborhood by
#' the squared matrix of dimension \eqn{radius*radius}. Finally, it returns for each cell the mean of the
#' weighted sum.
#'
#' @return It returns a \code{matrix} object with the same dimensions of \code{dataMatrix}.
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
#' myOutput <- meanFilter2D(myMatrix, radius)
#'
#' # Plot results
#' image(myMatrix, zlim = c(0, 100))
#' image(myOutput, zlim = c(0, 100))
meanFilter2D <- function(dataMatrix, radius, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(dataMatrix, radius, times)
  checkedArgs <- checkArgs(checkedArgs, type = "meanFilter2D")

  # Apply filters
  output <- dataMatrix
  for(i in seq(times)){
    output <- with(checkedArgs,
                   meanFilter2D_internal(dataMatrix = output, radius = radius))
  }

  return(output)
}
