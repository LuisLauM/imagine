#include <Rcpp.h>
#include <algorithm>
#include <math.h>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine
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

          if((a < 0) | (a > (nrows - 1)) | (b < 0) | (b > (ncols - 1)) | (std::isnan(data(a, b)))){
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
