#include <Rcpp.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine
// [[Rcpp::export]]
NumericMatrix engine1(NumericMatrix data, NumericMatrix kernel){
  double nrows = data.nrow();
  double ncols = data.ncol();

  double knlrows = kernel.nrow();
  double knlcols = kernel.ncol();

  // double rowY = (knlrows - 1)/2;
  // double colY = (knlcols - 1)/2;
  //
  // double rowZ = std::floor(rowY);
  // double colZ = std::ceil(colY);

  NumericMatrix emptyData(nrows, ncols);

  // for(double j = colZ - 1; j < ncols - colZ; j++){
  for(double j = 0; j < ncols; j++){
    for(double i = 0; i < nrows; i++){

      double cumSum = 0;
      double k = 0;

      for(double n = 0; n < knlcols; n++){
        for(double m = 0; m < knlrows; m++){
          double a = i + m - 1;
          double b = j + n - 1;

          if(!std::isnan(data(a, b))){
            cumSum += data(a, b)*kernel(m, n);
            k += kernel(m, n);
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
