#include <Rcpp.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine, .registration = TRUE

// ENGINE 1: 2D convolution
// [[Rcpp::export]]
NumericMatrix engine1(NumericMatrix data, NumericMatrix kernel, bool noNA){

  // Define dimension of matrix
  int nrows = data.nrow();
  int ncols = data.ncol();

  int knlrows = kernel.nrow();
  int knlcols = kernel.ncol();

  double knlRowHalfDouble = std::floor(knlrows/2);
  int knlRowHalf = (int)round(knlRowHalfDouble);

  double knlColHalfDouble = std::floor(knlcols/2);
  int knlColHalf = (int)round(knlColHalfDouble);

  int threshold = (knlrows*knlcols - 1);

  // If noNA is TRUE, define threshold as prod of dims of kernel
  if(noNA){
    threshold = 0;
  }

  // Define output matrix, same dims of input
  NumericMatrix emptyData(nrows, ncols);

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      if((i > (knlRowHalf - 1)) && (i < (nrows - knlRowHalf)) && (j > (knlColHalf - 1)) && (j < (ncols - knlColHalf))){
        double cumSum = 0;
        int naSum = 0;

        // Multiply the value of each cell by the corresponding value of the kernel.
        for(int n = 0; n < knlcols; n++){
          for(int m = 0; m < knlrows; m++){
            int a = i + m - knlRowHalf;
            int b = j + n - knlColHalf;

            // If the value is a NA do not consider for sum and increase the naSum index
            if(std::isnan(data(a, b))){
              naSum++;
            }else{
              cumSum += data(a, b)*kernel(m, n);
            }
          }
        }

        // Assign sum of values to corresponding cell, if all the values were NA, result will be NA
        if(naSum > threshold){
          emptyData(i, j) = NA_REAL;
        }else{
          emptyData(i, j) = cumSum;
        }
      }else{
        emptyData(i, j) = NA_REAL;
      }
    }
  }

  return emptyData;
}

// ENGINE 2: Convolution with quantiles
// [[Rcpp::export]]
NumericMatrix engine2(NumericMatrix data, NumericMatrix kernel, double probs){

  int nrows = data.nrow();
  int ncols = data.ncol();

  int knlrows = kernel.nrow();
  int knlcols = kernel.ncol();

  double knlRowHalfDouble = std::floor(knlrows/2);
  int knlRowHalf = (int)round(knlRowHalfDouble);

  double knlColHalfDouble = std::floor(knlcols/2);
  int knlColHalf = (int)round(knlColHalfDouble);

  NumericMatrix emptyData(nrows, ncols);
  NumericVector miniMatrix(knlrows*knlcols);

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      if((i > knlRowHalf) && (i < (nrows - knlRowHalf)) && (j > knlColHalf) && (j < (ncols - knlColHalf))){
        for(int n = 0; n < knlcols; n++){
          for(int m = 0; m < knlrows; m++){
            int index = m*knlcols + n;
            int a = i + m - knlRowHalf;
            int b = j + n - knlColHalf;

            if(std::isnan(data(a, b))){
              miniMatrix[index] = NA_REAL;
            }else{
              miniMatrix[index] = data(a, b)*kernel(m, n);
            }
          }
        }

        // Sort values
        miniMatrix.sort();

        // Get value for position indicated by 'probs'
        emptyData(i, j) = miniMatrix[probs];
      }else{
        emptyData(i, j) = NA_REAL;
      }
    }
  }

  return emptyData;
}

// ENGINE 3: Mean filter
// [[Rcpp::export]]
NumericMatrix engine3(NumericMatrix data, int radius){
  int nrows = data.nrow();
  int ncols = data.ncol();

  NumericMatrix emptyData(nrows, ncols);

  double halfRadiusDouble = std::floor(radius/2);
  int halfRadius = (int)round(halfRadiusDouble);

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      if((i > halfRadius) && (i < (nrows - halfRadius)) && (j > halfRadius) && (j < (ncols - halfRadius))){
        double cumSum = 0;
        int k = 1;

        for(int n = 0; n < radius; n++){
          for(int m = 0; m < radius; m++){

            int a = i + m - halfRadius;
            int b = j + n - halfRadius;

            if(!std::isnan(data(a, b))){
              cumSum += data(a, b);
              k++;
            }
          }
        }

        if(k < 1){
          emptyData(i, j) = NA_REAL;
        }else{
          emptyData(i, j) = cumSum/k;
        }
      }else{
        emptyData(i, j) = NA_REAL;
      }

    }
  }

  return emptyData;
}

// ENGINE 4: Quantile filter
// [[Rcpp::export]]
NumericMatrix engine4(NumericMatrix data, int radius, double probs){
  int nrows = data.nrow();
  int ncols = data.ncol();

  NumericMatrix emptyData(nrows, ncols);
  NumericVector miniMatrix(radius*radius);

  double halfRadiusDouble = std::floor(radius/2);
  int halfRadius = (int)round(halfRadiusDouble);

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      // Only if i and j is within limits, apply filter
      if((i > halfRadius) && (i < (nrows - halfRadius)) && (j > halfRadius) && (j < (ncols - halfRadius))){
        for(int n = 0; n < radius; n++){
          for(int m = 0; m < radius; m++){
            int index = m*radius + n;
            int a = i + m - halfRadius;
            int b = j + n - halfRadius;

            if(std::isnan(data(a, b))){
              miniMatrix[index] = NA_REAL;
            }else{
              miniMatrix[index] = data(a, b);
            }
          }
        }

        // Sort values
        miniMatrix.sort();

        // Get value for position indicated by 'probs'
        emptyData(i, j) = miniMatrix[probs];
      }else{
        emptyData(i, j) = NA_REAL;
      }
    }
  }

  return emptyData;
}

// ENGINE 5: Contextual Median Filter
// Proposed by Belkin et al. (2009), doi:10.1016/j.jmarsys.2008.11.018
// [[Rcpp::export]]
NumericMatrix engine5(NumericMatrix data, double probs, int I_radius, int O_radius){
  int nrows = data.nrow();
  int ncols = data.ncol();

  NumericMatrix emptyData(nrows, ncols);
  NumericVector O_miniMatrix(O_radius*O_radius);
  NumericVector I_miniMatrix(I_radius*I_radius);

  double I_halfRadiusDouble = std::floor(I_radius/2);
  int I_halfRadius = (int)round(I_halfRadiusDouble);

  double O_halfRadiusDouble = std::floor(O_radius/2);
  int O_halfRadius = (int)round(O_halfRadiusDouble);

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      // Only if i and j is within limits, apply filter
      if((i > O_halfRadius) && (i < (nrows - O_halfRadius)) && (j > O_halfRadius) && (j < (ncols - O_halfRadius))){

        // Check Outer filter
        int naCounter = 0;
        for(int n = 0; n < O_radius; n++){
          for(int m = 0; m < O_radius; m++){
            int index = m*O_radius + n;
            int a = i + m - O_halfRadius;
            int b = j + n - O_halfRadius;

            if(std::isnan(data(a, b))){
              O_miniMatrix[index] = NA_REAL;
              naCounter++;
            }else{
              O_miniMatrix[index] = data(a, b);
            }
          }
        }

        if(naCounter == (O_radius*O_radius)){
          emptyData(i, j) = NA_REAL;
        }else{
          double maxMini = max(O_miniMatrix);
          double minMini = min(O_miniMatrix);

          // Check the Outer filter
          // If TRUE, keep the grid value
          if((data(i, j) == minMini) | (data(i, j) == maxMini)){
            emptyData(i, j) = data(i, j);
          }else{
            // If FALSE, apply the Inner filter
            // Check Inner filter
            naCounter = 0;
            for(int n = 0; n < I_radius; n++){
              for(int m = 0; m < I_radius; m++){
                int index = m*I_radius + n;
                int a = i + m - I_halfRadius;
                int b = j + n - I_halfRadius;

                if(std::isnan(data(a, b))){
                  I_miniMatrix[index] = NA_REAL;
                  naCounter++;
                }else{
                  I_miniMatrix[index] = data(a, b);
                }
              }
            }

            if(naCounter == (I_radius*I_radius)){
              emptyData(i, j) = NA_REAL;
            }else{
              double maxMini = max(I_miniMatrix);
              double minMini = min(I_miniMatrix);

              if((data(i, j) == minMini) | (data(i, j) == maxMini)){
                // Sort values
                I_miniMatrix.sort();

                // Put the median value of the Inner miniMatrix
                emptyData(i, j) = I_miniMatrix[probs];
              }else{
                emptyData(i, j) = data(i, j);
              }
            }
          }
        }
      }else{
        emptyData(i, j) = NA_REAL;
      }
    }
  }

  return emptyData;
}
