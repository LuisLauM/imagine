#include <Rcpp.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine, .registration = TRUE

// ENGINE 1
// [[Rcpp::export]]
NumericMatrix engine1(NumericMatrix data, NumericMatrix kernel){

  // Engine 1: Basic convolution

  int nrows = data.nrow();
  int ncols = data.ncol();

  int knlrows = kernel.nrow();
  int knlcols = kernel.ncol();

  NumericMatrix emptyData(nrows, ncols);

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      double cumSum = 0;

      for(int n = 0; n < knlcols; n++){
        for(int m = 0; m < knlrows; m++){
          int a = i + m - 1;
          int b = j + n - 1;

          // Multiply the value of each cell by the corresponding value of the kernel.
          if((a > 1) && (a < (nrows - 1)) && (b > 1) && (b < (ncols - 1)) && (!std::isnan(data(a, b)))){
            cumSum += data(a, b)*kernel(m, n);
          }
        }
      }

      // Assign sum of values to corresponding cell
      emptyData(i, j) = cumSum;
    }
  }

  return emptyData;
}

// ENGINE 2
// [[Rcpp::export]]
NumericMatrix engine2(NumericMatrix data, NumericMatrix kernel){

  // Engine 2: Convolution mean
  // For now, this function is not available because the output is equivalent to
  // multiplying the matrices by a constant.

  int nrows = data.nrow();
  int ncols = data.ncol();

  int knlrows = kernel.nrow();
  int knlcols = kernel.ncol();

  NumericMatrix emptyData(nrows, ncols);

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      double cumSum = 0;
      int k = 1;

      for(int n = 0; n < knlcols; n++){
        for(int m = 0; m < knlrows; m++){
          int a = i + m - 1;
          int b = j + n - 1;

          // Multiply the value of each cell by the corresponding value of the kernel.
          if((a > 1) && (a < (nrows - 1)) && (b > 1) && (b < (ncols - 1)) && (!std::isnan(data(a, b)))){
            cumSum += data(a, b)*kernel(m, n);
            k = k + kernel(m, n);
          }

        }
      }

      // If all values were NA, returns NA for this cell
      if(k < 1){
        emptyData(i, j) = NA_REAL;
      }else{
        emptyData(i, j) = cumSum/k;
      }

    }
  }

  return emptyData;
}

// ENGINE 3
// [[Rcpp::export]]
NumericMatrix engine3(NumericMatrix data, NumericMatrix kernel, double x, double maxValue){

  // Engine 3: Convolution with quantiles

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

          if((a < 1 ) || (a > (nrows - 1)) || (b < 1) || (b > (ncols - 1)) || (std::isnan(data(a, b)))){
            miniMatrix[index] = maxValue;
          }else{
            miniMatrix[index] = data(a, b)*kernel(m, n);
          }
        }
      }

      // Sort values
      std::sort(miniMatrix.begin(), miniMatrix.end());

      // Get value for position indicated by 'x'
      emptyData(i, j) = miniMatrix[x];

    }
  }

  return emptyData;
}

// ENGINE 4
// [[Rcpp::export]]
NumericMatrix engine4(NumericMatrix data, int radius){
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

// ENGINE 5
// [[Rcpp::export]]
NumericMatrix engine5(NumericMatrix data, int radius, double x, double maxValue){
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

          if((a < 1) || (a > (nrows - 1)) || (b < 1) || (b > (ncols - 1)) || (std::isnan(data(a, b)))){
            miniMatrix[index] = maxValue;
          }else{
            miniMatrix[index] = data(a, b);
          }
        }
      }

      // Sort values
      std::sort(miniMatrix.begin(), miniMatrix.end());

      // Get value for position indicated by 'x'
      emptyData(i, j) = miniMatrix[x];
    }
  }

  return emptyData;
}
