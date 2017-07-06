#include <Rcpp.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine, .registration = TRUE

// ENGINE 1: 2D convolution
// [[Rcpp::export]]
NumericMatrix engine1(NumericMatrix data, NumericMatrix kernel, bool noNA = false){

  // Define dimension of matrix
  int nrows = data.nrow();
  int ncols = data.ncol();

  int knlrows = kernel.nrow();
  int knlcols = kernel.ncol();

  double knlRowHalfDouble = std::floor(knlrows/2);
  int knlRowHalf = (int)round(knlRowHalfDouble);

  double knlColHalfDouble = std::floor(knlcols/2);
  int knlColHalf = (int)round(knlColHalfDouble);

  bool threshold = 1;

  // If noNA is TRUE, define threshold as prod of dims of kernel
  if(!noNA){
    threshold = (knlrows*knlcols - 1);
  }

  // Define output matrix, same dims of input
  NumericMatrix emptyData(nrows, ncols);

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      if((i > knlRowHalf) && (i < (nrows - knlRowHalf)) && (j > knlColHalf) && (j < (ncols - knlColHalf))){
        double cumSum = 0;
        double naSum = 0;

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
NumericMatrix engine2(NumericMatrix data, NumericMatrix kernel, double x){

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

        // Get value for position indicated by 'x'
        emptyData(i, j) = miniMatrix[x];
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
NumericMatrix engine4(NumericMatrix data, int radius, double x){
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

        // Get value for position indicated by 'x'
        emptyData(i, j) = miniMatrix[x];
      }else{
        emptyData(i, j) = NA_REAL;
      }
    }
  }

  return emptyData;
}
