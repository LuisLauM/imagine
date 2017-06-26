#include <Rcpp.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine, .registration = TRUE

// ENGINE 1: 2D convolution
// [[Rcpp::export]]
NumericMatrix engine1(NumericMatrix data, NumericMatrix kernel){

  // Define dimension of matrix
  int nrows = data.nrow();
  int ncols = data.ncol();

  int knlrows = kernel.nrow();
  int knlcols = kernel.ncol();

  // Define output matrix, same dims of input
  NumericMatrix emptyData(nrows, ncols);

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      double cumSum = 0;
      double naSum = 0;

      // Multiply the value of each cell by the corresponding value of the kernel.
      for(int n = 0; n < knlcols; n++){
        for(int m = 0; m < knlrows; m++){
          int a = i + m - 1;
          int b = j + n - 1;

          // If the value is a NA do not consider for sum and increase the naSum index
          if((a > 1) && (a < (nrows - 1)) && (b > 1) && (b < (ncols - 1))){
            if(std::isnan(data(a, b))){
              naSum++;
            }else{
              cumSum += data(a, b)*kernel(m, n);
            }
          }
        }
      }

      // Assign sum of values to corresponding cell, if all the values were NA, result will be NA
      if(naSum < (knlrows*knlcols)){
        emptyData(i, j) = cumSum;
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

  NumericMatrix emptyData(nrows, ncols);
  NumericVector miniMatrix(knlrows*knlcols);

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      for(int n = 0; n < knlcols; n++){
        for(int m = 0; m < knlrows; m++){
          int index = m*knlcols + n;
          int a = i + m - 1;
          int b = j + n - 1;

          if((a > 1) && (a < (nrows - 1)) && (b > 1) && (b < (ncols - 1))){
            if(std::isnan(data(a, b))){
              miniMatrix[index] = NA_REAL;
            }else{
              miniMatrix[index] = data(a, b)*kernel(m, n);
            }
          }
        }
      }

      // Sort values
      miniMatrix.sort();

      // Get value for position indicated by 'x'
      emptyData(i, j) = miniMatrix[x];

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

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      double cumSum = 0;
      int k = 1;

      for(int n = 0; n < radius; n++){
        for(int m = 0; m < radius; m++){

          int a = i + m - 1;
          int b = j + n - 1;

          if((a > 1) && (a < (nrows - 1)) && (b > 1) && (b < (ncols - 1)) && (!std::isnan(data(a, b)))){
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

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      for(int n = 0; n < radius; n++){
        for(int m = 0; m < radius; m++){
          int index = m*radius + n;
          int a = i + m - 1;
          int b = j + n - 1;

          if((a > 1) && (a < (nrows - 1)) && (b > 1) && (b < (ncols - 1))){
            if(std::isnan(data(a, b))){
              miniMatrix[index] = NA_REAL;
            }else{
              miniMatrix[index] = data(a, b);
            }
          }
        }
      }


      // Sort values
      miniMatrix.sort();

      // Get value for position indicated by 'x'
      emptyData(i, j) = miniMatrix[x];
    }
  }

  return emptyData;
}
