# Mean filter ----------------------------------------------------------------------------

# Generate example matrix
nRows <- 100
nCols <- 200

myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)
radius <- 3

# Make convolution
myOutput <- meanFilter2D(myMatrix, radius)

# Plot results
image(myMatrix, zlim = c(0, 100))
image(myOutput, zlim = c(0, 100))
