
# Single filters ----------------------------------------------------------

# Generate example matrix
nRows <- 100
nCols <- 200

myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)
radius <- 3

# Mean filter
myOutput1 <- meanFilter(myMatrix, radius)

# Quantile filter
myOutput2 <- quantileFilter(myMatrix, radius, 0.1)

# Median filter
myOutput3 <- medianFilter(myMatrix, radius)

# Plot results
cols <- colorRampPalette(colors = c("black", "white"))(1e4)
par(mfrow = c(2, 2))

image(myMatrix, zlim = c(0, 100))
image(myOutput1, zlim = c(0, 100))
image(myOutput2, zlim = c(0, 100))
image(myOutput3, zlim = c(0, 100))
