#include <Rcpp.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine
// [[Rcpp::export]]
NumericMatrix engine4(NumericMatrix data, int radius){
  int nrows = data.nrow();
  int ncols = data.ncol();

  NumericMatrix emptyData(nrows, ncols);

  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){

      double cumSum = 0;
      int k = 0;

      for(int n = 0; n < radius; n++){
        for(int m = 0; m < radius; m++){

          int a = i + m - 1;
          int b = j + n - 1;

          if((a > 0) & (a < (nrows - 1)) & (b > 0) & (b < (ncols - 1)) & (!std::isnan(data(a, b)))){
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
