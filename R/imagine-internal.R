checkArgs <- function(imagineArgs, type){

  output <- switch(type,
                   convolution2D       = checkArgs_convolution2D(imagineArgs),
                   convolutionQuantile = checkArgs_convolutionQuantile(imagineArgs),
                   meanFilter          = checkArgs_meanFilter(imagineArgs),
                   quantileFilter      = checkArgs_quantileFilter(imagineArgs),
                   contextualMF        = checkArgs_contextualMF(imagineArgs),
                   agenbagFilters      = checkArgs_agenbagFilters(imagineArgs),
                   "Invalid value for 'type'.")

  return(output)
}


checkArgs_convolution2D <- function(allArgs){

  X <- allArgs$X
  kernel <- allArgs$kernel
  times <- allArgs$times

  if(!is.matrix(X) || !is.numeric(X)){
    stop("'X' must be a numeric matrix.")
  }

  if(!is.matrix(kernel) || !is.numeric(kernel)){
    stop("'kernel' must be a numeric matrix.")
  }

  if((nrow(kernel) >= nrow(X)) || ncol(kernel) >= ncol(X)){
    stop("Dimensions of 'kernel' matrix must be less than the 'X'.")
  }

  if(length(times) != 1 || !is.numeric(times)){
    stop("'time' must be a numeric vector with length 1.")
  }

  if(!isTRUE(all.equal(times, as.integer(times)))){
    warning("'times' must be an integer number. It will be coerced as integer.")
    allArgs$times <- as.integer(allArgs$times)
  }

  return(allArgs)
}

checkArgs_convolutionQuantile <- function(allArgs){

  X <- allArgs$X
  kernel <- allArgs$kernel
  probs <- allArgs$probs
  times <- allArgs$times

  if(!is.matrix(X) || !is.numeric(X)){
    stop("'X' must be a numeric matrix.")
  }

  if(any(is.na(kernel))){
    stop("'kernel' must not have NA inside.")
  }

  if(!is.matrix(kernel) || !is.numeric(kernel)){
    stop("'kernel' must be a numeric matrix.")
  }

  if((nrow(kernel) >= nrow(X)) || ncol(kernel) >= ncol(X)){
    stop("Dimensions of 'kernel' matrix must be less than the 'X'.")
  }

  if(!is.numeric(probs) || length(probs) != 1 || probs < 0 || probs > 1){
    stop("'probs' must be a numeric between 0 and 1 with length = 1.")
  }

  if(length(times) != 1 || !is.numeric(times)){
    stop("'time' must be a numeric vector with length 1.")
  }

  if(!isTRUE(all.equal(times, as.integer(times)))){
    warning("'times' must be an integer number. It will be coerced as integer.")
    allArgs$times <- as.integer(allArgs$times)
  }

  return(allArgs)
}

checkArgs_meanFilter <- function(allArgs){

  X <- allArgs$X
  radius <- allArgs$radius
  times <- allArgs$times

  if(!is.matrix(X) || !is.numeric(X)){
    stop("'X' must be a numeric matrix.")
  }

  if(!is.atomic(radius) || !is.numeric(radius)){
    stop("'radius' must be a numeric vector of length 1 or 2.")
  }

  if(radius[1] > nrow(X) || radius[2] > ncol(X)){
    stop("'radius' must be less than the dimensions (number of row and columns) of 'X'.")
  }

  if(!isTRUE(all.equal(radius, as.integer(radius)))){
    warning("'radius' must be an integer number. It will be coerced as integer.")
    allArgs$radius <- as.integer(allArgs$radius)
  }

  if(length(times) != 1 || !is.numeric(times)){
    stop("'time' must be a numeric vector with length 1.")
  }

  if(!isTRUE(all.equal(times, as.integer(times)))){
    warning("'times' must be an integer number. It will be coerced as integer.")
    allArgs$times <- as.integer(allArgs$times)
  }

  times <- as.integer(times)

  return(allArgs)
}

checkArgs_quantileFilter <- function(allArgs){

  X <- allArgs$X
  radius <- allArgs$radius
  probs <- allArgs$probs
  times <- allArgs$times

  if(!is.matrix(X) || !is.numeric(X)){
    stop("'X' must be a numeric matrix.")
  }

  if(!is.atomic(radius) || !is.numeric(radius)){
    stop("'radius' must be a numeric vector of length 1 or 2.")
  }

  if(!isTRUE(all.equal(radius, as.integer(radius)))){
    warning("'radius' must be an integer vector. It will be coerced as integer.")
    allArgs$radius <- as.integer(allArgs$radius)
  }

  if(any(radius < 1) || radius[1] > nrow(X) || radius[2] > ncol(X)){
    stop("Some value in 'radius' is out of acceptable bounds in regards to 'X'.")
  }

  if(!is.numeric(probs) || length(probs) != 1 || probs < 0 || probs > 1){
    stop("'probs' must be a numeric between 0 and 1 with length = 1.")
  }

  if(length(times) != 1 || !is.numeric(times)){
    stop("'time' must be a numeric vector with length 1.")
  }

  if(!isTRUE(all.equal(times, as.integer(times)))){
    warning("'times' must be an integer number. It will be coerced as integer.")
    allArgs$times <- as.integer(allArgs$times)
  }

  return(allArgs)
}

checkArgs_contextualMF <- function(allArgs){
  X <- allArgs$X

  if(!is.matrix(X) || !is.numeric(X)){
    stop("'X' must be a numeric matrix.")
  }

  return(allArgs)
}

checkArgs_agenbagFilters <- function(allArgs){

  X <- allArgs$X
  algorithm <- as.numeric(allArgs$algorithm[1])

  if(!is.matrix(X) || !is.numeric(X)){
    stop("'X' must be a numeric matrix.")
  }

  if(length(algorithm) != 1 || is.na(algorithm) || !algorithm %in% 1:2){
    stop("'algorithm' argument only accepts next numeric values: 1 or 2.")
  }

  return(allArgs)
}
