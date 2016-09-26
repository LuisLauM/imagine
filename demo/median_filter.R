# Median filter --------------------------------------------------------------------------

# Generate example matrix
nRows <- 1000
nCols <- 2000

myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)
radius <- 3

# Make convolution
myOutput <- medianFilter2D(myMatrix, radius)

# Plot results
image(myMatrix, zlim = c(0, 100))
image(myOutput, zlim = c(0, 100))
