# 2D convolution -------------------------------------------------------------------------
####### First example: A complete filled data matrix #######
# Generate example matrix
nRows <- 1000
nCols <- 2000

myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)
radius <- 3
kernel <- matrix(c(2, 1, 2,
                   1, 1, 1,
                   2, 1, 2), nrow = 3)

# Make convolution
myOutput <- convolution2D(myMatrix, kernel)

# Plot results
image(myMatrix, zlim = c(0, 100), col = colorRampPalette(c("red", "blue"))(1e3))
image(myOutput, zlim = c(0, 100), col = colorRampPalette(c("red", "blue"))(1e3))

####### Second example: A data matrix with NA values #######
# Generate example matrix
nRows <- 100
nCols <- 200

myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)

# Add some NA random values
index <- sample(x = seq(nRows*nCols), size = as.integer(nRows*nCols*0.2), replace = FALSE)
myMatrix[index] <- NA

# Build kernel
radius <- 3
kernel <- matrix(c(2, 1, 2,
                   1, 1, 1,
                   2, 1, 2), nrow = 3)

# Make convolution
myOutput <- convolution2D(myMatrix, kernel, times = 3)

# Plot results
image(myMatrix, zlim = c(0, 100), col = colorRampPalette(c("red", "blue"))(1e3))
image(myOutput, zlim = c(0, 100), col = colorRampPalette(c("red", "blue"))(1e3))
