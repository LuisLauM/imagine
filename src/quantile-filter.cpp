#include <algorithm>
#include <math.h>
#include <RcppArmadillo.h>
#include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @importFrom RcppArmadillo armadillo_version
//' @useDynLib imagine, .registration = TRUE

// ENGINE 4: Quantile filter
// [[Rcpp::export]]
NumericMatrix engine4_quantileFilter(arma::mat data, NumericVector radius, arma::vec probs){

  // Get dimension of input matrix
  int nrows = data.n_rows;
  int ncols = data.n_cols;

  int radius_row = radius[0];
  int radius_col = radius[1];

  NumericMatrix emptyData(nrows, ncols);
  std::fill(emptyData.begin(), emptyData.end(), NA_REAL);

  arma::mat miniMatrix(radius_row, radius_col);

  double halfRadiusDouble = std::floor(radius_row/2);
  int halfRadius_row = (int)round(halfRadiusDouble);

  halfRadiusDouble = std::floor(radius_col/2);
  int halfRadius_col = (int)round(halfRadiusDouble);

  for(int j = halfRadius_col; j < (ncols - halfRadius_col); j++){
    for(int i = halfRadius_row; i < (nrows - halfRadius_row); i++){

      int naCounter = 0;
      for(int n = 0; n < radius_col; n++){
        for(int m = 0; m < radius_row; m++){

          int a = i + m - halfRadius_row;
          int b = j + n - halfRadius_col;

          miniMatrix(m, n) = data(a, b);

          if(std::isnan(data(a, b))){
            naCounter++;
          }
        }
      }

      if(naCounter < (radius_row*radius_col)){
        arma::mat miniMatrix2 = arma::reshape(miniMatrix, miniMatrix.n_elem, 1);

        arma::mat quantVal = arma::quantile(miniMatrix2.elem(arma::find_finite(miniMatrix2)), probs);

        emptyData(i, j) = arma::conv_to < double >::from(quantVal);
      }
    }
  }

  return emptyData;
}
