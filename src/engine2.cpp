#include <Rcpp.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine
// [[Rcpp::export]]
NumericMatrix engine2(NumericMatrix data, NumericMatrix kernel, double x){
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
  NumericVector miniMatrix(knlrows*knlcols);

  // for(double j = colZ - 1; j < ncols - colZ; j++){
  for(double j = 0; j < ncols; j++){
    for(double i = 0; i < nrows; i++){
      for(double n = 0; n < knlcols; n++){
        for(double m = 0; m < knlrows; m++){
          double index = m*knlrows + n;
          double a = i + m - 1;
          double b = j + n - 1;

          miniMatrix[index] = data(a, b);
        }
      }

      std::sort(miniMatrix.begin(), miniMatrix.end());

      emptyData(i, j) = miniMatrix[x];
    }
  }

  return emptyData;
}
