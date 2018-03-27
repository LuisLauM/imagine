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

  X <- allArgs$X
  kernel <- allArgs$kernel
  times <- allArgs$times
  noNA <- allArgs$noNA

  if(class(X) != "matrix" || mode(X) != "numeric"){
    stop("'datMatrix' must be a numeric matrix.")
  }

  if(class(kernel) != "matrix" || mode(kernel) != "numeric"){
    stop("'kernel' must be a numeric matrix.")
  }

  if((nrow(kernel) >= nrow(X)) || ncol(kernel) >= ncol(X)){
    stop("Dimensions of 'kernel' matrix must be less than the 'X'.")
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

  X <- allArgs$X
  kernel <- allArgs$kernel
  probs <- allArgs$probs
  times <- allArgs$times

  if(class(X) != "matrix" || mode(X) != "numeric"){
    stop("'datMatrix' must be a numeric matrix.")
  }

  if(class(kernel) != "matrix" || mode(kernel) != "numeric"){
    stop("'kernel' must be a numeric matrix.")
  }

  if((nrow(kernel) >= nrow(X)) || ncol(kernel) >= ncol(X)){
    stop("Dimensions of 'kernel' matrix must be less than the 'X'.")
  }

  kernel <- matrix(data = round(an(kernel), 0), nrow = nrow(kernel))

  if(class(probs) != "numeric" || length(probs) != 1 || probs < 0 || probs > 1){
    stop("'probs' must be a numeric between 0 and 1 with length = 1.")
  }

  if(length(times) != 1 || mode(times) != "numeric"){
    stop("'time' must be a numeric vector with length 1.")
  }

  times <- as.integer(times)

  return(allArgs)
}

checkArgs_meanFilter <- function(allArgs){

  X <- allArgs$X
  radius <- allArgs$radius
  times <- allArgs$times

  if(class(X) != "matrix" || mode(X) != "numeric"){
    stop("'datMatrix' must be a numeric matrix.")
  }

  if(length(radius) != 1 || mode(radius) != "numeric"){
    stop("'radius' must be a numeric vector of length 1.")
  }

  if(radius > nrow(X) || radius > ncol(X)){
    stop("'radius' must be less than the dimensions (number of row and columns) of 'X'.")
  }

  radius <- as.integer(radius)

  if(length(times) != 1 || mode(times) != "numeric"){
    stop("'time' must be a numeric vector with length 1.")
  }

  times <- as.integer(times)

  return(allArgs)
}

checkArgs_quantileFilter <- function(allArgs){

  X <- allArgs$X
  radius <- allArgs$radius
  probs <- allArgs$probs
  times <- allArgs$times

  if(class(X) != "matrix" || mode(X) != "numeric"){
    stop("'datMatrix' must be a numeric matrix.")
  }

  if(length(radius) != 1 || mode(radius) != "numeric"){
    stop("'radius' must be a numeric vector of length 1.")
  }

  if(radius > nrow(X) || radius > ncol(X)){
    stop("'radius' must be less than the dimensions (number of row and columns) of 'X'.")
  }

  radius <- as.integer(radius)

  if(class(probs) != "numeric" || length(probs) != 1 || probs < 0 || probs > 1){
    stop("'probs' must be a numeric between 0 and 1 with length = 1.")
  }

  if(length(times) != 1 || mode(times) != "numeric"){
    stop("'time' must be a numeric vector with length 1.")
  }

  times <- as.integer(times)

  return(allArgs)
}

checkArgs_contextualMF <- function(allArgs){
  X <- allArgs$X
  inner_radius <- allArgs$inner_radius
  outer_radius <- allArgs$outer_radius
  probs <- allArgs$probs
  times <- allArgs$times

  if(class(X) != "matrix" || mode(X) != "numeric"){
    stop("'datMatrix' must be a numeric matrix.")
  }

  if(length(inner_radius) != 1 || mode(inner_radius) != "numeric" ||
     length(outer_radius) != 1 || mode(outer_radius) != "numeric"){
    stop("'outer_radius' and 'inner_radius' must be numeric vectors of length 1.")
  }

  if(inner_radius > nrow(X) || inner_radius > ncol(X) ||
     outer_radius > nrow(X) || outer_radius > ncol(X)){
    stop("'outer_radius' and 'inner_radius' must be less than the dimensions (number of row and columns) of 'X'.")
  }

  if(inner_radius >= outer_radius){
    stop("'outer_radius' must be greater than 'inner_radius'.")
  }

  inner_radius <- as.integer(inner_radius)
  outer_radius <- as.integer(outer_radius)

  if(class(probs) != "numeric" || length(probs) != 1 || probs < 0 || probs > 1){
    stop("'probs' must be a numeric between 0 and 1 with length = 1.")
  }

  if(length(times) != 1 || mode(times) != "numeric"){
    stop("'time' must be a numeric vector with length 1.")
  }

  times <- as.integer(times)

  return(allArgs)
}
