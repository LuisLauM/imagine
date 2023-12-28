#include <Rcpp.h>
#include <algorithm>
#include <math.h>
#include <fstream>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine, .registration = TRUE

// ENGINE 1: 2D convolution
// [[Rcpp::export]]
NumericMatrix engine1_2dConv(NumericMatrix data, NumericMatrix kernel){

  // Get dimension of input matrix
  int nrows = data.nrow();
  int ncols = data.ncol();

  // Get dimension of input kernel
  int knlrows = kernel.nrow();
  int knlcols = kernel.ncol();

  // Get half row size for kernel
  double knlRowHalfDouble = std::floor(knlrows/2);
  int knlRowHalf = (int)round(knlRowHalfDouble);

  // Get half column size for kernel
  double knlColHalfDouble = std::floor(knlcols/2);
  int knlColHalf = (int)round(knlColHalfDouble);

  // int threshold = knlrows*knlcols;
  int threshold = 1;

  // Define output matrix, same dims of input
  NumericMatrix emptyData(nrows, ncols);
  std::fill(emptyData.begin(), emptyData.end(), NA_REAL);

  for(int j = knlColHalf; j < (ncols - knlColHalf); j++){
    for(int i = knlRowHalf; i < (nrows - knlRowHalf); i++){

      double cumSum = 0;
      int naCounter = 0;

      // Multiply the value of each cell by the corresponding value of the kernel.
      for(int n = 0; n < knlcols; n++){
        for(int m = 0; m < knlrows; m++){
          int a = i + m - knlRowHalf;
          int b = j + n - knlColHalf;

          // If the value is a NA do not consider for sum and increase the naCounter
          if(std::isnan(data(a, b))){
            naCounter++;
          }else{
            cumSum += data(a, b)*kernel(m, n);
          }
        }
      }

      // Assign sum of values to corresponding cell, if all the values were NA,
      // result will be NA
      if(naCounter < threshold){
        emptyData(i, j) = cumSum;
      }
    }
  }

  return emptyData;
}
