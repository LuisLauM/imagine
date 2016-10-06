#include <Rcpp.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine
// [[Rcpp::export]]
NumericMatrix engine2(NumericMatrix data, NumericMatrix kernel){

  // Engine 2: Convolution mean

  int nrows = data.nrow();
  int ncols = data.ncol();

  int knlrows = kernel.nrow();
  int knlcols = kernel.ncol();

  NumericMatrix emptyData(nrows, ncols);

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      double cumSum = 0;
      int k = 0;

      for(int n = 0; n < knlcols; n++){
        for(int m = 0; m < knlrows; m++){
          int a = i + m - 1;
          int b = j + n - 1;

          // Multiply the value of each cell by the corresponding value of the kernel.
          if((a > 0) & (a < (nrows - 1)) & (b > 0) & (b < (ncols - 1)) & (!std::isnan(data(a, b)))){
            cumSum += data(a, b)*kernel(m, n);
            k += kernel(m, n);
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
