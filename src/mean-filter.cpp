#include <Rcpp.h>
#include <algorithm>
#include <math.h>
#include <fstream>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine, .registration = TRUE

// ENGINE 3: Mean filter
// [[Rcpp::export]]
NumericMatrix engine3_meanFilter(NumericMatrix data, NumericVector radius){

  // Get dimension of input matrix
  int nrows = data.nrow();
  int ncols = data.ncol();

  // Get dimension of input kernel
  int radius_row = radius[0];
  int radius_col = radius[1];

  // Get half row size for kernel
  double halfRadiusDouble = std::floor(radius_row/2);
  int halfRadius_row = (int)round(halfRadiusDouble);

  // Get half column size for kernel
  halfRadiusDouble = std::floor(radius_col/2);
  int halfRadius_col = (int)round(halfRadiusDouble);

  // Define an empty matrix as large as the original where the output values
  // will be storaged
  NumericMatrix emptyData(nrows, ncols);
  std::fill(emptyData.begin(), emptyData.end(), NA_REAL);

  // Loop along every cell of original matrix
  for(int j = halfRadius_col; j < (ncols - halfRadius_col); j++){
    for(int i = halfRadius_row; i < (nrows - halfRadius_row); i++){

      // Initialize cumsum of window values (cumSum) and no-NA values (k)
      double cumSum = 0;
      int k = 0;

      for(int n = 0; n < radius_col; n++){
        for(int m = 0; m < radius_row; m++){

          int a = i + m - halfRadius_row;
          int b = j + n - halfRadius_col;

          if(!std::isnan(data(a, b))){
            cumSum += data(a, b);
            k++;
          }
        }
      }

      // Only if there is one or more of no-NA values
      if(k > 0){
        emptyData(i, j) = cumSum/k;
      }
    }
  }

  return emptyData;
}
