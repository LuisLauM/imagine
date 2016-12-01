# 2D convolution -------------------------------------------------------------------------
####### First example: The effect of convolution filters #######
myMatrix <- wbImage
kernel1 <- matrix(c(-1, -2, -1,
                    0,  0,  0,
                    1,  2,  1), nrow = 3)

# Make convolution
myOutput1 <- convolution2D(myMatrix, kernel)
myOutput2 <- convolutionMean(myMatrix, kernel)
myOutput3 <- convolutionQuantile(myMatrix, kernel, x = 0.7)

# Plot results
cols <- colorRampPalette(colors = c("black", "white"))(1e4)
par(mfrow = c(2, 2), mar = rep(0, 4))
image(myMatrix, zlim = c(0, 1), axes = FALSE, col = cols)
image(myOutput1, zlim = c(0, 1), axes = FALSE, col = cols)
image(myOutput2, zlim = c(0, 1), axes = FALSE, col = cols)
image(myOutput3, zlim = c(0, 1), axes = FALSE, col = cols)

####### Second example: A data matrix with NA values #######
# Generate example matrix
nRows <- 300
nCols <- 300

myMatrix <- matrix(seq(from = 0, to = 100, length.out = nRows*nCols),
                   nrow = nRows, ncol = nCols)

# Add some NA random values
index <- sample(x = seq(nRows*nCols), size = as.integer(nRows*nCols*0.7), replace = FALSE)
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
cols <- colorRampPalette(colors = c("red", "green"))(1e4)
par(mfrow = c(2, 2), mar = rep(0, 4))
image(myMatrix, zlim = c(0, 100), axes = FALSE, col = cols)
image(myOutput1, zlim = c(0, 100), axes = FALSE, col = cols)
image(myOutput2, zlim = c(0, 100), axes = FALSE, col = cols)
image(myOutput3, zlim = c(0, 100), axes = FALSE, col = cols)

####### Third example: Influence of number of applications (times) #######
# Generate example matrix
nRows <- 300
nCols <- 600

myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)

# Add some NA random values
index <- sample(x = seq(nRows*nCols), size = as.integer(nRows*nCols*0.7), replace = FALSE)
myMatrix[index] <- NA

# Build kernel
radius <- 3
kernel <- matrix(c(2, 1, 2, 1, 2,
                   1, 1, 1, 1, 1,
                   2, 1, 2, 2, 1), nrow = 3)

# Make convolution
myOutput1 <- convolutionMean(myMatrix, kernel, times = 1)
myOutput2 <- convolutionMean(myMatrix, kernel, times = 3)
myOutput3 <- convolutionMean(myMatrix, kernel, times = 5)

# Plot results
cols <- colorRampPalette(colors = c("red", "green"))(1e4)
par(mfrow = c(2, 2), mar = rep(0, 4))
image(myMatrix, zlim = c(0, 100), axes = FALSE, col = cols)
image(myOutput1, zlim = c(0, 100), axes = FALSE, col = cols)
image(myOutput2, zlim = c(0, 100), axes = FALSE, col = cols)
image(myOutput3, zlim = c(0, 100), axes = FALSE, col = cols)


####### Fourth example: Reconstruct image #######
# Build an example matrix
origMatrix <- wbImage

# Add some NAs
origMatrix_withNA <- origMatrix
origMatrix_withNA[sample(seq(prod(dim(origMatrix))), 0.7*prod(dim(origMatrix)), replace = FALSE)] <- NA

# Apply filter
newMatrix <- quantileFilter(dataMatrix = origMatrix_withNA, radius = 3, x = 0.2, times = 4)

# Plot matrices for compare
cols <- colorRampPalette(c("white", "black"))(n)

par(mar = c(0, 2, 0, 0), mfrow = c(3, 1))

image(origMatrix, col = cols, axes = FALSE)
mtext(text = "Original", side = 2, line = 0.5, font = 2)

image(origMatrix_withNA, col = cols, axes = FALSE)
mtext(text = "Original with NAs", side = 2, line = 0.5, font = 2)

image(newMatrix, col = cols, axes = FALSE)
mtext(text = "Filtered", side = 2, line = 0.5, font = 2)
