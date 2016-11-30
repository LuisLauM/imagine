# 2D convolution -------------------------------------------------------------------------
####### First example: A complete filled data matrix #######
# Generate example matrix
nRows <- 50
nCols <- 100

myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)
radius <- 3
kernel <- matrix(c(2, 1, 2,
                   1, 1, 1,
                   2, 1, 2), nrow = 3)

# Make convolution
myOutput1 <- convolution2D(myMatrix, kernel)
myOutput2 <- convolutionMean(myMatrix, kernel)
myOutput3 <- convolutionQuantile(myMatrix, kernel, x = 0.7)

# Plot results
cols <- colorRampPalette(colors = c("black", "white"))(1e4)
par(mfrow = c(2, 2))
image(myMatrix, zlim = c(0, 100), col = cols)
image(myOutput1, zlim = c(0, 100), col = cols)
image(myOutput2, zlim = c(0, 100), col = cols)
image(myOutput3, zlim = c(0, 100), col = cols)

####### Second example: A data matrix with NA values #######
# Generate example matrix
nRows <- 50
nCols <- 100

myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)

# Add some NA random values
index <- sample(x = seq(nRows*nCols), size = as.integer(nRows*nCols*0.2), replace = FALSE)
myMatrix[index] <- NA

# Build kernel
radius <- 7
kernel <- matrix(c(2, 1, 2, 1, 2,
                   1, 1, 1, 1, 1,
                   2, 1, 2, 2, 1), nrow = 3)

# Make convolution
myOutput1 <- convolution2D(myMatrix, kernel)
myOutput2 <- convolutionMean(myMatrix, kernel)
myOutput3 <- convolutionQuantile(myMatrix, kernel, x = 0.7)

# Plot results
par(mfrow = c(2, 2))
image(myMatrix, zlim = c(0, 100), col = cols)
image(myOutput1, zlim = c(0, 100), col = cols)
image(myOutput2, zlim = c(0, 100), col = cols)
image(myOutput3, zlim = c(0, 100), col = cols)
