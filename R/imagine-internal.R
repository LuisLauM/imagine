convolution2D_internal <- function(dataMatrix, kernel){

  newData <- engine1(data = dataMatrix, kernel = kernel)

  return(newData)
}

convolutionMean_internal <- function(dataMatrix, kernel){

  newData <- engine2(data = dataMatrix, kernel = kernel)

  return(newData)
}

convolutionQuantile_internal <- function(dataMatrix, kernel, x){

  x <- ceiling(x*prod(dim(kernel))) - 1
  maxValue <- max(an(dataMatrix), na.rm = TRUE)*max(an(kernel), na.rm = TRUE)

  newData <- engine3(data = dataMatrix, kernel = kernel, x = x, maxValue = maxValue)
  newData[newData > maxValue - 1] <- NA

  return(newData)
}

meanFilter_internal <- function(dataMatrix, radius){

  newData <- engine4(data = dataMatrix, radius = radius)

  return(newData)
}

quantileFilter_internal <- function(dataMatrix, radius, x){

  x <- ceiling(x*radius^2) - 1
  maxValue <- max(an(dataMatrix), na.rm = TRUE)

  newData <- engine5(data = dataMatrix, radius = radius, x = as.integer(x), maxValue = maxValue)
  newData[newData > maxValue - 1] <- NA

  return(newData)
}

checkArgs <- function(imagineArgs, type){

  output <- switch(type,
                   convolution2D       = checkArgs_convolution2D(imagineArgs),
                   convolutionMean     = checkArgs_convolutionMean(imagineArgs),
                   convolutionQuantile = checkArgs_convolutionQuantile(imagineArgs),
                   meanFilter          = checkArgs_meanFilter(imagineArgs),
                   quantileFilter      = checkArgs_quantileFilter(imagineArgs),
                   "Invalid value for 'type'.")

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

checkArgs_convolutionMean <- function(allArgs){

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

  if(sum(an(kernel)) <= 0){
    stop("For 'convolutionMean', kernel values sum must be higher than zero.")
  }

  if(length(times) != 1 || mode(times) != "numeric"){
    stop("'time' must be a numeric vector with length 1.")
  }

  times <- as.integer(times)

  return(allArgs)
}

checkArgs_convolutionQuantile <- function(allArgs){

  dataMatrix <- allArgs$dataMatrix
  kernel <- allArgs$kernel
  x <- allArgs$x
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

  if(class(x) != "numeric" || length(x) != 1 || x < 0 || x > 1){
    stop("'x' must be a numeric between 0 and 1 with length = 1.")
  }

  if(length(times) != 1 || mode(times) != "numeric"){
    stop("'time' must be a numeric vector with length 1.")
  }

  times <- as.integer(times)

  return(allArgs)
}

checkArgs_meanFilter <- function(allArgs){

  dataMatrix <- allArgs$dataMatrix
  radius <- allArgs$radius
  times <- allArgs$times

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

checkArgs_quantileFilter <- function(allArgs){

  dataMatrix <- allArgs$dataMatrix
  radius <- allArgs$radius
  x <- allArgs$x
  times <- allArgs$times

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

  if(class(x) != "numeric" || length(x) != 1 || x < 0 || x > 1){
    stop("'x' must be a numeric between 0 and 1 with length = 1.")
  }

  if(length(times) != 1 || mode(times) != "numeric"){
    stop("'time' must be a numeric vector with length 1.")
  }

  times <- as.integer(times)

  return(allArgs)
}
