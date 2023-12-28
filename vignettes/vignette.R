## ----echo=FALSE---------------------------------------------------------------
library(imagine)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("imagine")

## ----eval=FALSE---------------------------------------------------------------
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

## ----message=FALSE, fig.height=6, fig.width=5.33, fig.cap = "Figure 1: 2D vs 2D quantile convolutions", results='hide', fig.pos="h", echo=FALSE----

# Defining a copy of wbImage
myMatrix <- wbImage

# Defining color palette
cols <- gray.colors(n = 1e3, start = 0, end = 1)

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

## ----eval=FALSE---------------------------------------------------------------
#  # Add some noise (NA) to the image (matrix)
#  set.seed(7)
#  naIndex <- sample(x       = seq(prod(dim(myMatrix))),
#                    size    = as.integer(0.4*prod(dim(myMatrix))),
#                    replace = FALSE)
#  myMatrix[naIndex] <- NA
#  
#  # Build kernel
#  radius <- 3
#  
#  # Apply filters
#  meanfilterExample     <- meanFilter(X = myMatrix, radius = radius)
#  quantilefilterExample <- quantileFilter(X = myMatrix, radius = radius, probs = 0.1)
#  medianfilterExample   <- medianFilter(X = myMatrix, radius = radius)
#  

## ----message=FALSE, fig.height=5, fig.width=7.5, fig.cap = "Figure 2: Basic filters comparison", results='hide', fig.pos="h", echo=FALSE----
# Defining a copy of wbImage
myMatrix <- wbImage

# Defining color palette
cols <- gray.colors(n = 1e3, start = 0, end = 1)

# Add some noise (NA) to the image (matrix)
set.seed(7)
naIndex <- sample(x = seq(prod(dim(myMatrix))), size = as.integer(0.4*prod(dim(myMatrix))))
myMatrix[naIndex] <- NA

# Apply filters
meanfilterExample     <- meanFilter(X = myMatrix, radius = 3)
medianfilterExample   <- medianFilter(X = myMatrix, radius = 3)

# Make plots
par(mar = rep(0, 4), oma = rep(0.5, 4), mfrow = c(2, 2))

image(wbImage, col = cols, axes = FALSE)
mtext(text = "Original", side = 3, line = -1.5, font = 2, adj = 0.99)

image(myMatrix, col = cols, axes = FALSE)
mtext(text = "Original with noise (NA)", side = 3, line = -1.5, font = 2, adj = 0.99)

# meanfilterExample[meanfilterExample < 0] <- 0
image(meanfilterExample, col = cols, axes = FALSE)
mtext(text = "Mean filter", side = 3, line = -1.5, font = 2, adj = 0.99)

# medianfilterExample[medianfilterExample < 0] <- 0
image(medianfilterExample, col = cols, axes = FALSE)
mtext(text = "2D median filter", side = 3, line = -1.5, font = 2, adj = 0.99)

## ----eval=FALSE---------------------------------------------------------------
#  medianFilter(X = wbImage, radius = 5, times = 50)

## ----message=FALSE, fig.height=5, fig.width=7.5, fig.cap = "Figure 3: Filters with several time settings", results='hide', fig.pos="h", echo=FALSE----
# Defining color palette
cols <- gray.colors(n = 1e3, start = 0, end = 1)

# Apply filters
median_times01   <- medianFilter(X = wbImage, radius = 5)
median_times10   <- medianFilter(X = wbImage, radius = 5, times = 10)
median_times50   <- medianFilter(X = wbImage, radius = 5, times = 50)

# Make plots
par(mar = rep(0, 4), oma = rep(0.5, 4), mfrow = c(2, 2))

image(wbImage, col = cols, axes = FALSE)
mtext(text = "Original", side = 3, line = -1.5, font = 2, adj = 0.99)

image(median_times01, col = cols, axes = FALSE)
mtext(text = "2D median filter\ntimes = 01", side = 3, line = -2.5, 
      font = 2, adj = 0.99)

image(median_times10, col = cols, axes = FALSE)
mtext(text = "2D median filter\ntimes = 10", side = 3, line = -2.5, 
      font = 2, adj = 0.99)

image(median_times50, col = cols, axes = FALSE)
mtext(text = "2D median filter\ntimes = 50", side = 3, line = -2.5, 
      font = 2, adj = 0.99)

