## ---- echo = FALSE, message = FALSE-------------------------------------------
library(imagine)

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("imagine")

## ---- eval=FALSE--------------------------------------------------------------
#  # Build kernels
#  # Kernel 1: For bottom edge recognition
#  kernel1 <- matrix(c(-1, -2, -1,
#                       0,  0,  0,
#                       1,  2,  1),
#                    nrow = 3)
#  
#  # Kernel 2: Diagonal weighting
#  kernel2 <- matrix(c(-2, 0, 0,
#                       0, 1, 0,
#                       0, 0, 2),
#                    nrow = 3)
#  
#  # Apply filters
#  convolutionExample  <- convolution2D(X = wbImage, kernel = kernel1)
#  convQuantileExample <- convolutionQuantile(X = wbImage, kernel = kernel2, probs = 0.1)
#  

## ---- message=FALSE, fig.height=3, fig.width=5.33, fig.cap = "Figure 2: Original matrix", results='hide', fig.pos="h", echo=FALSE----
par(mar = rep(0, 4), mfrow = c(1, 1))
cols <- colorRampPalette(colors = c("black", "white"))(1e4)

image(wbImage, col = cols)

## ---- message=FALSE, fig.height=6, fig.width=5.33, fig.cap = "Figure 1: Filtered matrices", results='hide', fig.pos="h", echo=FALSE----

myMatrix <- wbImage

# Build kernels
# Kernel 1: For bottom edge recognition
kernel1 <- matrix(c(-1, -2, -1,
                     0,  0,  0,
                     1,  2,  1), 
                  nrow = 3)

# Kernel 2: Diagonal weighting
kernel2 <- matrix(c(-2, 0, 0,
                    0, 1, 0,
                    0, 0, 2), 
                  nrow = 3)

# Apply filters
convolutionExample  <- convolution2D(X = myMatrix, kernel = kernel1)
convQuantileExample <- convolutionQuantile(X = myMatrix, kernel = kernel2, probs = 0.1)

# Make plots
par(mar = c(0, 0.5, 0, 0.5), oma = c(0, 0, 2, 0), mfrow = c(2, 1))

image(convolutionExample, col = cols, axes = FALSE)
mtext(text = "2D convolution", side = 1, line = -1.5, col = "white", font = 2, adj = 0.99)

image(convQuantileExample, col = cols, axes = FALSE)
mtext(text = "2D quantile convolution", side = 1, line = -1.5, col = "black", font = 2, adj = 0.99)

## ---- eval=FALSE--------------------------------------------------------------
#  # Add some noise (NA) to the image (matrix)
#  set.seed(7)
#  naIndex <- sample(x = seq(prod(dim(myMatrix))), size = as.integer(0.4*prod(dim(myMatrix))), replace = FALSE)
#  myMatrix[naIndex] <- NA
#  
#  # Build kernel
#  radius <- 3
#  
#  # Apply filters
#  meanfilterExample     <- meanFilter(X = myMatrix, radius = radius)
#  quantilefilterExample <- quantileFilter(X = myMatrix, radius = radius, probs = 0.1)
#  medianfilterExample   <- medianFilter(X = myMatrix, radius = radius, times = 10)
#  

## ---- message=FALSE, fig.height=3, fig.width=5.33, fig.cap = "Figure 2: Original matrix", results='hide', fig.pos="h", echo=FALSE----
set.seed(7)
naIndex <- sample(x = seq(prod(dim(myMatrix))), size = as.integer(0.4*prod(dim(myMatrix))), replace = FALSE)
myMatrix[naIndex] <- NA

par(mar = rep(0, 4), mfrow = c(1, 1))
image(myMatrix, col = cols)

## ---- message=FALSE, fig.height=9, fig.width=5.33, fig.cap = "Figure 1: Filtered matrices", results='hide', fig.pos="h", echo=FALSE----
# Build kernel
radius <- 3

# Add some noise (NA) to the image (matrix)
set.seed(7)
naIndex <- sample(x = seq(prod(dim(myMatrix))), size = as.integer(0.4*prod(dim(myMatrix))), replace = FALSE)
myMatrix[naIndex] <- NA

# Build kernel
radius <- 3

# Apply filters
meanfilterExample     <- meanFilter(X = myMatrix, radius = radius)
quantilefilterExample <- quantileFilter(X = myMatrix, radius = radius, probs = 0.1)
medianfilterExample   <- medianFilter(X = myMatrix, radius = radius, times = 10)

# Make plots
par(mar = c(0, 0.5, 0, 0.5), oma = c(0, 0, 2, 0), mfrow = c(3, 1))

# meanfilterExample[meanfilterExample < 0] <- 0
image(meanfilterExample, col = cols, axes = FALSE)
mtext(text = "Mean filter", side = 1, line = -1.5, col = "white", font = 2, adj = 0.99)

# quantilefilterExample[quantilefilterExample < 0] <- 0
image(quantilefilterExample, col = cols, axes = FALSE)
mtext(text = "Quantile filter (probs = 0.1)", side = 1, line = -1.5, col = "white", font = 2, adj = 0.99)

# medianfilterExample[medianfilterExample < 0] <- 0
image(medianfilterExample, col = cols, axes = FALSE)
mtext(text = "2D median filter", side = 1, line = -1.5, col = "white", font = 2, adj = 0.99)

