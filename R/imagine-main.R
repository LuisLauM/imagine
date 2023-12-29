#' @title Data matrix to be used as example image.
#' @description \code{matrix} object containing numeric data to plot a image.
#' The photo was taken by the author on 2016.
#' @docType data
#' @usage wbImage
#' @format A \code{matrix} with dimensions 1280x720.
"wbImage"


#' @title Make convolution calculations from numeric matrix
#'
#' @rdname convolutions
#'
#' @param X A numeric \code{matrix} object used for apply filters.
#' @param kernel A little matrix used as mask for each cell of \code{X}.
#' @param probs \code{numeric} vector of probabilities with values in [0,1].
#' @param times How many times do you want to apply the filter?
#' @param normalize \code{logical} indicating if results will (or not) be
#' normalized. See details.
#'
#' @description This function takes a \code{matrix} object, and for each cell
#' multiplies its neighborhood by the \code{kernel}. Finally, it returns for
#' each cell the mean of the kernel-weighted sum.
#'
#' @return \code{convolution2D} returns a \code{matrix} object with the same
#' dimensions of \code{X}.
#'
#' @details
#' Convolution is a mathematical operation that combines two arrays of numbers
#' to produce an array of the same structure. The output will consist of only
#' valid values, meaning those where both arrays have non-missing data.
#' Consequently, any missing values (NAs) in the input matrix will propagate
#' outwards to the extent of the convolution kernel.
#'
#' Through normalization, the output of each convolution window is scaled by
#' dividing it by the sum of the absolute values of the kernel
#' (\code{sum(abs(as.numeric(kernel)))}, disabled by default).
#'
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
convolution2D <- function(X, kernel, times = 1, normalize = FALSE){

  armadillo_version(FALSE)

  # Check and validation of arguments
  checkedArgs <- list(X = X, kernel = kernel, times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = "convolution2D")

  # Apply filters
  output <- checkedArgs$X
  for(i in seq(checkedArgs$times)){
    gc(reset = TRUE)

    output <- with(checkedArgs, engine1_2dConv(data = output, kernel = kernel))

    if(normalize){
      output <- output/sum(abs(as.numeric(kernel)), na.rm = TRUE)
    }
  }

  output
}

#' @rdname convolutions
#' @return \code{convolutionQuantile} uses the kernel but, for each cell, it
#' returns the position of quantile 'probs' (value between 0 and 1).
#' @export
convolutionQuantile <- function(X, kernel, probs, times = 1, normalize = FALSE){

  # Check and validation of arguments
  checkedArgs <- list(X = X, kernel = kernel, probs = probs, times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = "convolutionQuantile")

  # Apply filters
  output <- checkedArgs$X
  for(i in seq(checkedArgs$times)){
    gc(reset = TRUE)

    output <- with(checkedArgs,
                   engine2_convWithQuantiles(data = output, kernel = kernel, probs = probs))

    if(normalize){
      output <- output/sum(abs(as.numeric(kernel)), na.rm = TRUE)
    }
  }

  output
}


#' @rdname convolutions
#' @return \code{convolutionMedian} is a wrapper of \code{convolutionQuantile}
#' with probs = 0.5.
#' @export
convolutionMedian <- function(X, kernel, times = 1){

  convolutionQuantile(X = X, kernel = kernel, probs = 0.5, times = times)
}


#' @title Make a 2D filter calculations from numeric matrix
#'
#' @rdname basic2DFilter
#'
#' @param X A numeric \code{matrix} object used for apply filters.
#' @param radius Size of squared or rectangular kernel to apply median. See
#' Details.
#' @param probs \code{numeric} vector of probabilities with values in [0,1].
#' @param times How many times do you want to apply the filter?
#'
#' @description This functions take a \code{matrix} object, and for each cell
#' calculate mean, median or certain quantile around a squared/rectangular
#' neighborhood.
#'
#' @return A \code{matrix} object with the same dimensions of \code{X}.
#'
#' @details \code{radius} must be defined as a 2-length numeric vector
#' specifying the number of rows and columns of the window which will be used to
#' make calculations. If the length of radius is 1, the window will be a square.
#'
#' Functions use C++ algorithms for running some statistical calculations. The
#' mean is far obvious, however, there are several ways to perform quantiles.
#' \code{quantileFilter} function uses
#' \href{https://arma.sourceforge.net/docs.html#quantile}{arma::quantile}: a
#' RcppArmadillo function, which is equivalent to use R \link[stats]{quantile}
#' funtion with \code{type = 5}.
#'
#' \code{medianFilter} is a wraper of \code{quantileFilter}, so the same
#' observations are applied to it.
#'
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
#' myOutput1 <- meanFilter(X = myMatrix, radius = radius)
#' myOutput2 <- quantileFilter(X = myMatrix, radius = radius, probs = 0.1)
#' myOutput3 <- medianFilter(X = myMatrix, radius = radius)
#'
#' # Plot results
#' par(mfrow = c(2, 2))
#' image(myMatrix, zlim = c(0, 100), title = "Original")
#' image(myOutput1, zlim = c(0, 100), title = "meanFilter")
#' image(myOutput2, zlim = c(0, 100), title = "quantileFilter")
#' image(myOutput3, zlim = c(0, 100), title = "medianFilter")
meanFilter <- function(X, radius, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(X = X,
                      radius = rep(x = radius, length.out = 2),
                      times = times)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = "meanFilter")

  # Apply filters
  output <- checkedArgs$X
  for(i in seq(checkedArgs$times)){
    gc(reset = TRUE)

    output <- with(checkedArgs, engine3_meanFilter(data = output, radius = radius))
  }

  output
}


#' @rdname basic2DFilter
#' @return \code{quantileFilter} don't use a kernel but, for each cell, it
#' returns the position of quantile 'probs' (value between 0 and 1).
#' @export
quantileFilter <- function(X, radius, probs, times = 1){

  # Check and validation of arguments
  checkedArgs <- list(X = X,
                      radius = rep(x = radius, length.out = 2),
                      probs = probs,
                      times = times,
                      na = NA)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = "quantileFilter")

  # Apply filters
  output <- checkedArgs$X
  for(i in seq(checkedArgs$times)){
    gc(reset = TRUE)

    output <- with(checkedArgs,
                   engine4_quantileFilter(data = output, radius = radius, probs = probs))
  }

  output
}


#' @rdname basic2DFilter
#' @return \code{medianFilter} is a wrapper of \code{quantileFilter} with
#' \code{probs = 0.5}.
#' @export
medianFilter <- function(X, radius, times = 1){

  quantileFilter(X = X, radius = radius, probs = 0.5, times = times)
}


#' @title Performs Contextual Median Filter
#'
#' @param X A numeric \code{matrix} object used for apply filters.
#'
#' @description This function implements the Contextual Median Filter (CMF)
#' algorithm, which was first described by Belkin & O'Reilly (2009), following
#' the pseudocode provided in their paper.
#'
#' @details
#' Following the definition of CMF, since \strong{imagine} v.2.0.0, \code{times}
#' argument will not be available anymore.
#'
#' \strong{imagine} offers the CMF algorithm but for the using to find out
#' oceanographic fronts, it is recommended to see and use the functions of the
#' \href{https://CRAN.R-project.org/package=grec}{\strong{grec}} package.
#'
#' @references Belkin, I. M., & O'Reilly, J. E. (2009). An algorithm for oceanic
#' front detection in chlorophyll and SST satellite imagery. Journal of Marine
#' Systems, 78(3), 319-326 (\doi{http://dx.doi.org/10.1016/j.jmarsys.2008.11.018}).
#'
#' @export
#'
#' @return \code{contextualMF} returns a \code{matrix} object with the same
#' dimensions of \code{X}.
#'
#' @examples
#' data(wbImage)
#'
#' # Agenbag, gradient algorithm 1
#' cmdOut <- agenbagFilters(X = wbImage, algorithm = 1)
#'
#' # image(cmdOut)
contextualMF <- function(X){

  # Check and validation of arguments
  checkedArgs <- list(X = X)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = "contextualMF")

  # Apply filters
  with(checkedArgs, engine5_CMF(data = X))
}


#' @title Performs algorithms from Agenbag et al. (2003)
#'
#' @param X A numeric \code{matrix} used as main input.
#' @param algorithm \code{integer} specifying the type of method that will be
#' used. See Details.
#' @param ... Not used.
#'
#' @description This function performs two (gradient) calculation approaches for
#' SST, as outlined in the paper by Agenbag et al. (2003).
#'
#' @details
#' Section 2.2.4 of the paper by Agenbag et al. (2003) introduces the following
#' two methods:
#' \itemize{
#'  \item{\strong{Method 1: }}{Based on the equation
#'  \deqn{Y_{i,j}=\sqrt{(X_{i+1,j}-X_{i-1,j})^2 +(X_{i,j+1}-X_{i,j-1})^2}}}
#'  where \eqn{Y_{i,j}} represents the output value for each \eqn{X_{i,j}} pixel value
#'  of a given \eqn{X} matrix.
#'  \item{\strong{Method 2: }}{the standard deviation in a 3x3 pixel area centered on
#'  position \eqn{(i,j)}.}
#' }
#'
#' As outlined in the original study, this method conducts searches within a
#' 1-pixel vicinity of each point. For method 1, it only returns a value for
#' points where none of the four involved values are NA. Conversely, for method
#' 2, the standard deviation calculation is performed only for points where at
#' least 3 non-NA values are found in the 3x3 neighborhood.
#'
#'
#' @references Agenbag, J.J., A.J. Richardson, H. Demarcq, P. Freon, S. Weeks,
#' and F.A. Shillington. "Estimating Environmental Preferences of South African
#' Pelagic Fish Species Using Catch Size- and Remote Sensing Data". Progress in
#' Oceanography 59, No 2-3 (October 2003): 275-300.
#' (\doi{https://doi.org/10.1016/j.pocean.2003.07.004}).
#'
#' @return \code{agenbagFilters} returns a \code{matrix} object with the same
#' dimensions of \code{X}.
#'
#' @export
#'
#' @examples
#' data(wbImage)
#'
#' # Agenbag, method 1
#' agenbag1 <- agenbagFilters(X = wbImage, algorithm = 1)
#'
#' # Agenbag, method 2
#' agenbag2 <- agenbagFilters(X = wbImage, algorithm = 2)
#'
#' # Plotting results
#' par(mfrow = c(3, 1), mar = rep(0, 4))
#'
#' # Original
#' image(wbImage, axes = FALSE, col = gray.colors(n = 1e3))
#'
#' # Calculated
#' cols <- hcl.colors(n = 1e3, palette = "YlOrRd", rev = TRUE)
#' image(agenbag1, axes = FALSE, col = cols)
#' image(agenbag2, axes = FALSE, col = cols)
agenbagFilters <- function(X, algorithm = c(1, 2), ...){

  # Check and validation of arguments
  checkedArgs <- list(X = X,
                      algorithm = algorithm)
  checkedArgs <- checkArgs(imagineArgs = checkedArgs, type = "agenbagFilters")

  switch(checkedArgs$algorithm,
         "1" = engine6_agenbag1(data = checkedArgs$X),
         "2" = engine7_agenbag2(data = checkedArgs$X, ...))
}
