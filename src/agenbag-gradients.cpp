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

// ENGINE 6: Agenbag et al (2003)
// [[Rcpp::export]]
NumericMatrix engine6_agenbag1(NumericMatrix data){

  // Get dimension of input matrix
  int nrows = data.nrow();
  int ncols = data.ncol();

  // Get dimension of input kernel
  int radius_row = 3;
  int radius_col = 3;

  // Get half row size for kernel
  double halfRadiusDouble = std::floor(radius_row/2);
  int halfRadius_row = (int)round(halfRadiusDouble);

  // Get half column size for kernel
  halfRadiusDouble = std::floor(radius_col/2);
  int halfRadius_col = (int)round(halfRadiusDouble);

  // Define an empty matrix as large as the original where the output values
  // will be stored
  NumericMatrix emptyData(nrows, ncols);
  std::fill(emptyData.begin(), emptyData.end(), NA_REAL);

  // Loop along every cell of original matrix
  for(int j = halfRadius_col; j < (ncols - halfRadius_col); j++){
    for(int i = halfRadius_row; i < (nrows - halfRadius_row); i++){

      if(!std::isnan(data(i + 1, j)) && !std::isnan(data(i - 1, j)) && !std::isnan(data(i, j + 1)) && !std::isnan(data(i, j - 1))){
        double cellValue = pow(pow(data(i + 1, j) - data(i - 1, j), 2) + pow(data(i, j + 1) - data(i, j - 1), 2), 0.5);
        emptyData(i, j) = cellValue;
      }
    }
  }

  return emptyData;
}


// ENGINE 7: Agenbag et al (2003)
// [[Rcpp::export]]
NumericMatrix engine7_agenbag2(arma::mat data, int threshold = 3){

  // Get dimension of input matrix
  int nrows = data.n_rows;
  int ncols = data.n_cols;

  int radius_row = 3;
  int radius_col = 3;

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

      if((radius_row*radius_col - naCounter) >= threshold){
        arma::mat miniMatrix2 = arma::reshape(miniMatrix, miniMatrix.n_elem, 1);

        arma::mat quantVal = arma::stddev(miniMatrix2.elem(arma::find_finite(miniMatrix2)), 0, 0);

        emptyData(i, j) = arma::conv_to < double >::from(quantVal);
      }
    }
  }

  return emptyData;
}
