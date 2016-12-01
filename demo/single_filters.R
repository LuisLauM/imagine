
# Single filters ----------------------------------------------------------

# Generate example matrix
nRows <- 1000
nCols <- 1000

myMatrix <- matrix(runif(nRows*nCols, 0, 100), nrow = nRows, ncol = nCols)
radius <- 3

index <- sample(x = seq(prod(dim(myMatrix))), size = as.integer(0.5*prod(dim(myMatrix))), replace = FALSE)
myMatrix[index] <- NA

# Mean filter
myOutput1 <- meanFilter(myMatrix, radius, times = 3)

# Quantile filter
myOutput2 <- quantileFilter(myMatrix, radius, 0.9)

# Median filter
myOutput3 <- medianFilter(myMatrix, radius)

# Plot results
cols <- colorRampPalette(colors = c("red", "green"))(1e4)
par(mfrow = c(2, 2), mar = rep(0, 4))

image(myMatrix, zlim = c(0, 100), axes = FALSE, col = cols)
image(myOutput1, zlim = c(0, 100), axes = FALSE, col = cols)
image(myOutput2, zlim = c(0, 100), axes = FALSE, col = cols)
image(myOutput3, zlim = c(0, 100), axes = FALSE, col = cols)
