#include <Rcpp.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine
// [[Rcpp::export]]
NumericMatrix engine2(NumericMatrix data, double radius, double x){
  double nrows = data.nrow();
  double ncols = data.ncol();

  // double rowY = (knlrows - 1)/2;
  // double colY = (knlcols - 1)/2;
  //
  // double rowZ = std::floor(rowY);
  // double colZ = std::ceil(colY);

  NumericMatrix emptyData(nrows, ncols);
  NumericVector miniMatrix(radius*radius);

  // for(double j = colZ - 1; j < ncols - colZ; j++){
  for(double j = 0; j < ncols; j++){
    for(double i = 0; i < nrows; i++){
      for(double n = 0; n < radius; n++){
        for(double m = 0; m < radius; m++){
          double index = m*radius + n;
          double a = i + m - 1;
          double b = j + n - 1;

          miniMatrix[index] = data(a, b);
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
