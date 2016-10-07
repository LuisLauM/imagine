## ---- echo = FALSE, message = FALSE--------------------------------------
library(imagine)

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("imagine")

## ---- message=FALSE, fig.height=3, fig.width=5.33, fig.cap = "Figure 2: Original matrix", results='hide', fig.pos="h", echo=FALSE----
par(mar = rep(0, 4), mfrow = c(1, 1))
cols <- colorRampPalette(colors = c("black", "white"))(1e4)

image(myMatrix, col = cols)

## ---- message=FALSE, fig.height=9, fig.width=5.33, fig.cap = "Figure 1: Filtered matrices", results='hide', fig.pos="h", echo=FALSE----
par(mar = c(0, 0.5, 0, 0.5), oma = c(0, 0, 2, 0), mfrow = c(3, 1))

convolutionExample[convolutionExample < 0] <- 0
image(convolutionExample, col = cols, axes = FALSE)
mtext(text = "2D convolution", side = 1, line = -1.5, col = "white", font = 2, adj = 0.99)

convMeanExample[convMeanExample < 0] <- 0
image(convMeanExample, col = cols, axes = FALSE)
mtext(text = "2D mean convolution", side = 1, line = -1.5, col = "white", font = 2, adj = 0.99)

convQuantileExample[convQuantileExample < 0] <- 0
image(convQuantileExample, col = cols, axes = FALSE)
mtext(text = "2D median convolution", side = 1, line = -1.5, col = "white", font = 2, adj = 0.99)

## ---- cache=TRUE---------------------------------------------------------
myMatrix <- wbImage

# # Add some noise (NA) to the image (matrix)
# naIndex <- sample(x = seq(prod(dim(myMatrix))), size = as.integer(0.4*prod(dim(myMatrix))), replace = FALSE)
# myMatrix[naIndex] <- NA

# Build kernel
radius <- 3

# Apply filters
meanfilterExample     <- meanFilter(dataMatrix = myMatrix, radius = radius)
quantilefilterExample <- quantileFilter(dataMatrix = myMatrix, radius = radius, x = 0.1)
medianfilterExample   <- medianFilter(dataMatrix = myMatrix, radius = radius)


## ---- message=FALSE, fig.height=3, fig.width=5.33, fig.cap = "Figure 2: Original matrix", results='hide', fig.pos="h", echo=FALSE----
par(mar = rep(0, 4), mfrow = c(1, 1))
image(myMatrix, col = cols)

## ---- message=FALSE, fig.height=9, fig.width=5.33, fig.cap = "Figure 1: Filtered matrices", results='hide', fig.pos="h", echo=FALSE----
par(mar = c(0, 0.5, 0, 0.5), oma = c(0, 0, 2, 0), mfrow = c(3, 1))

meanfilterExample[meanfilterExample < 0] <- 0
image(meanfilterExample, col = cols, axes = FALSE)
mtext(text = "Mean filter", side = 1, line = -1.5, col = "white", font = 2, adj = 0.99)

quantilefilterExample[quantilefilterExample < 0] <- 0
image(quantilefilterExample, col = cols, axes = FALSE)
mtext(text = "Quantile filter", side = 1, line = -1.5, col = "white", font = 2, adj = 0.99)

medianfilterExample[medianfilterExample < 0] <- 0
image(medianfilterExample, col = cols, axes = FALSE)
mtext(text = "2D median filter", side = 1, line = -1.5, col = "white", font = 2, adj = 0.99)

## ---- message=FALSE, fig.height=2, fig.width=5.7, fig.cap = "Figure 3: Neighborhood kernel application for different kernel dimensions. Black dot indicates the position of the pixel over the filter will be applied. Arrows indicates the direction of filter application", results='hide', fig.pos="h", echo=FALSE----
par(mar = rep(0, 4), mfrow = c(1, 3), xaxs = "i", yaxs = "i")

xlim <- c(-0.15, 1.1)
ylim <- c(-0.01, 1.3)

plotKernels <- function(dim1 = 3, dim2 = 3, showArrows = FALSE){
  plot(1, 1, pch = NA, axes = FALSE, xlim = xlim, ylim = ylim, xlab = NA, ylab = NA)
  polygon(x = c(0, 1, 1, 0), y = c(0, 0, 1, 1))
  
  delay1 <- ifelse(dim1 %% 2 != 0, 0, 1/(2*dim1))
  delay2 <- ifelse(dim2 %% 2 != 0, 0, 1/(2*dim2))
  
  mtext(text = paste0(dim1, "x", dim2, " kernel"), side = 3, line = -2)
  
  for(i in seq(dim1 - 1)){
    lines(c(0, 1), 1 - rep(i/dim1, 2))
  }
  
  for(i in seq(dim2 - 1)){
    lines(rep(i/dim2, 2), c(0, 1))
  }
  
  # abline(h = seq(0, 1, length.out = dim1 + 1))
  # abline(v = seq(0, 1, length.out = dim2 + 1))
  
  if(isTRUE(showArrows)){
    arrows(x0 = -0.05, y0 = 1.05, x1 = -0.05, y1 = 0.8, length = 0.1)
    arrows(x0 = -0.05, y0 = 1.05, x1 = 0.25, y1 = 1.05, length = 0.1)
  }
  
  points(0.5 - delay2, 0.5 + delay1, pch = 16)
}

# 3x3 kernel
plotKernels(3, 3, showArrows = TRUE)

# 3x4 kernel
plotKernels(3, 4)

# 4x5 kernel
plotKernels(4, 5)

