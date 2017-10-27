checkArgs <- function(imagineArgs, type){

  output <- switch(type,
                   convolution2D       = checkArgs_convolution2D(imagineArgs),
                   convolutionQuantile = checkArgs_convolutionQuantile(imagineArgs),
                   meanFilter          = checkArgs_meanFilter(imagineArgs),
                   quantileFilter      = checkArgs_quantileFilter(imagineArgs),
                   contextualMF        = checkArgs_contextualMF(imagineArgs),
                   "Invalid value for 'type'.")

  return(output)
}


checkArgs_convolution2D <- function(allArgs){

  dataMatrix <- allArgs$dataMatrix
  kernel <- allArgs$kernel
  times <- allArgs$times
  noNA <- allArgs$noNA

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

  if(!is.logical(noNA)){
    noNA <- isTRUE(noNA)

    warning("Given 'noNA' was not logical. Default value (FALSE) was taken.")
  }

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

checkArgs_contextualMF <- function(allArgs){
  dataMatrix <- allArgs$dataMatrix
  inner_radius <- allArgs$inner_radius
  outer_radius <- allArgs$outer_radius
  x <- allArgs$x
  times <- allArgs$times

  if(class(dataMatrix) != "matrix" || mode(dataMatrix) != "numeric"){
    stop("'datMatrix' must be a numeric matrix.")
  }

  if(length(inner_radius) != 1 || mode(inner_radius) != "numeric" ||
     length(outer_radius) != 1 || mode(outer_radius) != "numeric"){
    stop("'outer_radius' and 'inner_radius' must be numeric vectors of length 1.")
  }

  if(inner_radius > nrow(dataMatrix) || inner_radius > ncol(dataMatrix) ||
     outer_radius > nrow(dataMatrix) || outer_radius > ncol(dataMatrix)){
    stop("'outer_radius' and 'inner_radius' must be less than the dimensions (number of row and columns) of 'dataMatrix'.")
  }

  if(inner_radius >= outer_radius){
    stop("'outer_radius' must be greater than 'inner_radius'.")
  }

  inner_radius <- as.integer(inner_radius)
  outer_radius <- as.integer(outer_radius)

  if(class(x) != "numeric" || length(x) != 1 || x < 0 || x > 1){
    stop("'x' must be a numeric between 0 and 1 with length = 1.")
  }

  if(length(times) != 1 || mode(times) != "numeric"){
    stop("'time' must be a numeric vector with length 1.")
  }

  times <- as.integer(times)

  return(allArgs)
}
