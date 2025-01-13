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

// ENGINE 2: Convolution with quantiles
// [[Rcpp::export]]
NumericMatrix engine2_convWithQuantiles(arma::mat data, arma::mat kernel, arma::vec probs, bool na_only){

  // Get dimension of input matrix
  int nrows = data.n_rows;
  int ncols = data.n_cols;

  // Get dimension of input kernel
  int knlrows = kernel.n_rows;
  int knlcols = kernel.n_cols;

  // Get half row size for kernel
  double knlRowHalfDouble = std::floor(knlrows/2);
  int knlRowHalf = (int)round(knlRowHalfDouble);

  // Get half column size for kernel
  double knlColHalfDouble = std::floor(knlcols/2);
  int knlColHalf = (int)round(knlColHalfDouble);

  // Define an empty matrix as large as the original where the output values
  // will be storaged
  NumericMatrix emptyData(nrows, ncols);
  std::fill(emptyData.begin(), emptyData.end(), NA_REAL);

  // Define a matrix as large as the kernel where the temporal convolution
  // output values will be storaged
  arma::mat miniMatrix(knlrows, knlcols);

  // Loop along every cell of original matrix
  for(int j = knlColHalf; j < (ncols - knlcols); j++){
    for(int i = knlRowHalf; i < (nrows - knlrows); i++){

      if(na_only){
        if(!std::isnan(data(i, j))){
          emptyData(i, j) = data(i, j);
          continue;
        }
      }

      // Initialize NA counter
      int naCounter = 0;

      // Loop along each window (kernel-sized) around each cell
      for(int n = 0; n < knlcols; n++){
        for(int m = 0; m < knlrows; m++){

          // Define relative indexes for the current window
          int a = i + m - knlRowHalf;
          int b = j + n - knlColHalf;

          // Calculate the product of original-matrix and kernel values
          miniMatrix(m, n) = data(a, b)*kernel(m, n);

          // If cell value is a NA, increasing the counter
          if(std::isnan(data(a, b))){
            naCounter++;
          }
        }
      }

      // Only if the ammount of NAs is lower than the # of elements of miniMatrix
      if(naCounter < (knlrows*knlcols)){
        // Reshape miniMatrix as an Armadillo matrix of 1 column
        arma::mat miniMatrix2 = arma::reshape(miniMatrix, miniMatrix.n_elem, 1);

        // Calculate the quantile
        arma::mat quantVal = arma::quantile(miniMatrix2.elem(arma::find_finite(miniMatrix2)), probs);

        // Replace the quantile value in the corresponding cell of the output matrix
        emptyData(i, j) = arma::conv_to < double >::from(quantVal);
      }
    }
  }

  return emptyData;
}
