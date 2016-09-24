convolution2D_internal <- function(dataMatrix, kernel){


  sumKernel <- sum(an(kernel))

  newData <- engine1(data = dataMatrix, kernel = kernel)

  return(newData)
}

medianFilter2D_internal <- function(dataMatrix, radius){

  kernel <- matrix(data = 1, nrow = radius, ncol = radius)

  x <- ceiling(0.5*prod(dim(kernel)))

  newData <- engine2(data = dataMatrix, kernel = kernel, x = x)

  return(newData)
}

meanFilter2D_internal <- function(dataMatrix, radius){

  kernel <- matrix(data = 1, nrow = radius, ncol = radius)

  newData <- engine1(data = dataMatrix, kernel = kernel)

  return(newData)
}

checkArgs <- function(imagineArgs, type){

  output <- switch(tolower(type),
                   convolution2D = checkArgs_convolution2D(imagineArgs),
                   medianFilter2D = checkArgs_medianFilter2D(imagineArgs),
                   meanFilter2D = checkArgs_medianFilter2D(imagineArgs))

  return(output)
}

checkArgs_convolution2D <- function(allArgs){

  dataMatrix <- allArgs$dataMatrix
  kernel <- allArgs$kernel
  times <- allArgs$times

  if(class(dataMatrix) != "matrix" || mode(dataMatrix) != "numeric"){
    stop("'datMatrix' must be a numeric matrix.")
  }

  if(class(kernel) != "matrix" || mode(kernel) != "numeric"){
    stop("'kernel' must be a numeric matrix.")
  }

  if((nrow(kernel) >= nrow(dataMatrix)) || ncol(kernel) >= ncol(dataMatrix)){
    stop("Dimensions of 'kernel' matrix must be less than the 'dataMatrix'.")
  }

  kernel <- matrix(data = round(an(kernel), 0), nrow = nrow(kernel))

  if(length(times) != 1 || mode(times) != "numeric"){
    stop("'time' must be a numeric vector with length 1.")
  }

  times <- as.integer(times)

  return(allArgs)
}

checkArgs_medianFilter2D <- function(allArgs){

  dataMatrix <- allArgs$dataMatrix
  radius <- allArgs$radius

  if(class(dataMatrix) != "matrix" || mode(dataMatrix) != "numeric"){
    stop("'datMatrix' must be a numeric matrix.")
  }

  if(length(radius) != 1 || mode(radius) != "numeric"){
    stop("'radius' must be a numeric vector of length 1.")
  }

  if(radius > nrow(dataMatrix) || radius > ncol(dataMatrix)){
    stop("'radius' must be less than the dimensions (number of row and columns) of 'dataMatrix'.")
  }

  radius <- as.integer(radius)

  if(length(times) != 1 || mode(times) != "numeric"){
    stop("'time' must be a numeric vector with length 1.")
  }

  times <- as.integer(times)

  return(allArgs)
}
