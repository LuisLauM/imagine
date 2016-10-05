#include <Rcpp.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine
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
          cumSum += data(a, b)*kernel(m, n);
        }
      }

      // Assign sum of values to corresponding cell
      emptyData(i, j) = cumSum;
    }
  }

  return emptyData;
}
