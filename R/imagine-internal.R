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
  na <- allArgs$na

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

  if(!is.na(na)){
    if(sum(is.na(X)) > 0){
      stop("If 'na' has been set as a value, 'X' must not contain any NA.")
    }

    if(length(na) != 1 || !is.numeric(na)){
      stop("'na' must be set whether as NA or as a numeric value of length 1.")
    }

    if(max(X, na.rm = TRUE) >= as.integer(na)){
      stop("'na' value must be an integer value higher than max(X).")
    }

    if(!isTRUE(all.equal(na, as.integer(na)))){
      warning("'na' must be an integer number. It will be coerced as integer.")
    }

    allArgs$naVal <- as.integer(allArgs$na)
  }else{
    allArgs$naVal <- as.integer(ceiling(max(X, na.rm = TRUE)) + 1)
    allArgs$X[is.na(allArgs$X)] <- allArgs$naVal
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

  if(length(radius) != 1 || !is.numeric(radius)){
    stop("'radius' must be a numeric vector of length 1.")
  }

  if(radius > nrow(X) || radius > ncol(X)){
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
  na <- allArgs$na

  if(!is.matrix(X) || !is.numeric(X)){
    stop("'X' must be a numeric matrix.")
  }

  if(length(radius) != 1 || !is.numeric(radius)){
    stop("'radius' must be a numeric vector of length 1.")
  }

  if(!isTRUE(all.equal(radius, as.integer(radius)))){
    warning("'radius' must be an integer number. It will be coerced as integer.")
    allArgs$radius <- as.integer(allArgs$radius)
  }

  if(radius > nrow(X) || radius > ncol(X)){
    stop("'radius' must be less than the dimensions (number of row and columns) of 'X'.")
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

  if(!is.na(na)){
    if(sum(is.na(X)) > 0){
      stop("If 'na' has been set as a value, 'X' must not contain any NA.")
    }

    if(length(na) != 1 || !is.numeric(na)){
      stop("'na' must be set whether as NA or as a numeric value of length 1.")
    }

    if(max(X, na.rm = TRUE) >= as.integer(na)){
      stop("'na' value must be an integer value higher than max(X).")
    }

    if(!isTRUE(all.equal(na, as.integer(na)))){
      warning("'na' must be an integer number. It will be coerced as integer.")
    }

    allArgs$naVal <- as.integer(allArgs$na)
  }else{
    allArgs$naVal <- as.integer(ceiling(max(X, na.rm = TRUE)) + 1)
    allArgs$X[is.na(allArgs$X)] <- allArgs$naVal
  }

  return(allArgs)
}

checkArgs_contextualMF <- function(allArgs){
  X <- allArgs$X
  # inner_radius <- allArgs$inner_radius
  # outer_radius <- allArgs$outer_radius
  # probs <- allArgs$probs
  times <- allArgs$times
  na <- allArgs$na

  if(!is.matrix(X) || !is.numeric(X)){
    stop("'X' must be a numeric matrix.")
  }

  # if(length(inner_radius) != 1 || !is.numeric(inner_radius) ||
  #    length(outer_radius) != 1 || !is.numeric(outer_radius)){
  #   stop("'outer_radius' and 'inner_radius' must be numeric vectors of length 1.")
  # }
  #
  # if(inner_radius > nrow(X) || inner_radius > ncol(X) ||
  #    outer_radius > nrow(X) || outer_radius > ncol(X)){
  #   stop("'outer_radius' and 'inner_radius' must be less than the dimensions (number of row and columns) of 'X'.")
  # }
  #
  # if(inner_radius >= outer_radius){
  #   stop("'outer_radius' must be greater than 'inner_radius'.")
  # }
  #
  # if(!isTRUE(all.equal(inner_radius, as.integer(inner_radius))) || !isTRUE(all.equal(outer_radius, as.integer(outer_radius)))){
  #   warning("'inner_radius' and 'outer_radius' must be integer numbers. They will be coerced as integers.")
  #   inner_radius <- as.integer(inner_radius)
  #   outer_radius <- as.integer(outer_radius)
  # }
  #
  # if(!is.numeric(probs) || length(probs) != 1 || probs < 0 || probs > 1){
  #   stop("'probs' must be a numeric vector with values between 0 and 1 with length = 1.")
  # }

  if(length(times) != 1 || !is.numeric(times)){
    stop("'time' must be a numeric vector with length 1.")
  }

  if(!isTRUE(all.equal(times, as.integer(times)))){
    warning("'times' must be an integer number. It will be coerced as integer.")
    allArgs$times <- as.integer(allArgs$times)
  }

  if(!is.na(na)){
    if(sum(is.na(X)) > 0){
      stop("If 'na' has been set as a value, 'X' must not contain any NA.")
    }

    if(length(na) != 1 || !is.numeric(na)){
      stop("'na' must be set whether as NA or as a numeric value of length 1.")
    }

    if(max(X, na.rm = TRUE) >= as.integer(na)){
      stop("'na' value must be an integer value higher than max(X).")
    }

    if(!isTRUE(all.equal(na, as.integer(na)))){
      warning("'na' must be an integer number. It will be coerced as integer.")
    }

    allArgs$naVal <- as.integer(allArgs$na)
  }else{
    allArgs$naVal <- as.integer(ceiling(max(X, na.rm = TRUE)) + 1)
    allArgs$X[is.na(allArgs$X)] <- allArgs$naVal
  }

  return(allArgs)
}
