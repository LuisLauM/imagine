// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <algorithm>
#include <math.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @importFrom RcppArmadillo armadillo_version
//' @useDynLib imagine, .registration = TRUE

// ENGINE 1: 2D convolution
// [[Rcpp::export]]
NumericMatrix engine1(NumericMatrix data, NumericMatrix kernel){

  // Get dimension of input matrix
  int nrows = data.nrow();
  int ncols = data.ncol();

  // Get dimension of input kernel
  int knlrows = kernel.nrow();
  int knlcols = kernel.ncol();

  // Get half row size for kernel
  double knlRowHalfDouble = std::floor(knlrows/2);
  int knlRowHalf = (int)round(knlRowHalfDouble);

  // Get half column size for kernel
  double knlColHalfDouble = std::floor(knlcols/2);
  int knlColHalf = (int)round(knlColHalfDouble);

  // int threshold = knlrows*knlcols;
  int threshold = 1;

  // Define output matrix, same dims of input
  NumericMatrix emptyData(nrows, ncols);
  std::fill(emptyData.begin(), emptyData.end(), NA_REAL);

  for(int j = knlColHalf; j < (ncols - knlColHalf); j++){
    for(int i = knlRowHalf; i < (nrows - knlRowHalf); i++){

      double cumSum = 0;
      int naCounter = 0;

      // Multiply the value of each cell by the corresponding value of the kernel.
      for(int n = 0; n < knlcols; n++){
        for(int m = 0; m < knlrows; m++){
          int a = i + m - knlRowHalf;
          int b = j + n - knlColHalf;

          // If the value is a NA do not consider for sum and increase the naCounter
          if(std::isnan(data(a, b))){
            naCounter++;
          }else{
            cumSum += data(a, b)*kernel(m, n);
          }
        }
      }

      // Assign sum of values to corresponding cell, if all the values were NA,
      // result will be NA
      if(naCounter < threshold){
        emptyData(i, j) = cumSum;
      }
    }
  }

  return emptyData;
}

// ENGINE 2: Convolution with quantiles
// [[Rcpp::export]]
NumericMatrix engine2(arma::mat data, arma::mat kernel, arma::vec probs){

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

// ENGINE 3: Mean filter
// [[Rcpp::export]]
NumericMatrix engine3(NumericMatrix data, NumericVector radius){

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

// ENGINE 4: Quantile filter
// [[Rcpp::export]]
NumericMatrix engine4(arma::mat data, NumericVector radius, arma::vec probs){
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

// ENGINE 5: Contextual Median Filter
// Proposed by Belkin et al. (2009), doi:10.1016/j.jmarsys.2008.11.018
// [[Rcpp::export]]
NumericMatrix engine5(NumericMatrix data){
  int nrows = data.nrow();
  int ncols = data.ncol();

  // Create a NA matrix for output
  NumericMatrix emptyData(nrows, ncols);
  std::fill(emptyData.begin(), emptyData.end(), NA_REAL);

  // OUTER INDEX
  // |  0  1  2  3  4 |
  // |  5  6  7  8  9 |
  // | 10 11 12 13 14 |
  // | 15 16 17 18 19 |
  // | 20 21 22 23 24 |

  // N-S split outer matrix
  IntegerVector ns_outer = IntegerVector::create(2, 7, 12, 17, 22);

  // W-E split outer matrix
  IntegerVector we_outer = IntegerVector::create(10, 11, 12, 13, 14);

  // NW-SE split outer matrix
  IntegerVector nwse_outer = IntegerVector::create(0, 6, 12, 18, 24);

  // NE-SW split outer matrix
  IntegerVector nesw_outer = IntegerVector::create(4, 8, 12, 16, 20);


  // INNER INDEX
  // |  1  2  3 |
  // |  4  5  6 |
  // |  7  8  9 |

  // N-S split outer matrix
  IntegerVector ns_inner = IntegerVector::create(2, 5, 8);

  // W-E split outer matrix
  IntegerVector we_inner = IntegerVector::create(4, 5, 6);

  // 1. Check for peaks and troughs within 1D 5-point slices through a sliding
  // 5×5 window. The window slides east–west (E–W), northsouth (N–S) across the
  // matrix:

  // Loop for the whole matrix
  for(int j = 2; j < (ncols - 2); j++){
    for(int i = 2; i < (nrows - 2); i++){
      // If cell data(i, j) is missing, avoid it
      if(!std::isnan(data(i, j))){

        // Create outputs for Outer and Inner mini matrix
        NumericVector O_miniMatrix(25);
        NumericVector I_miniMatrix(9);

        // Check NA on inner matrix
        int naCounter = 0;
        for(int n = 0; n < 3; n++){
          for(int m = 0; m < 3; m++){
            int index = m*3 + n;
            int a = i + m - 1;
            int b = j + n - 1;

            if(std::isnan(data(a, b))){
              // If cell value is missing (NA), increase the counter
              naCounter++;
            }else{
              // Otherwise, add to the miniMatrix
              I_miniMatrix[index] = data(a, b);
            }
          }
        }

        if(naCounter > 8){
          continue;
        }else{
          // Check NA on outer matrix
          for(int n = 0; n < 5; n++){
            for(int m = 0; m < 5; m++){
              int index = m*5 + n;
              int a = i + m - 1;
              int b = j + n - 1;

              if(std::isnan(data(a, b))){
                // If cell value is missing (NA), increase the counter
                naCounter++;
              }else{
                // Otherwise, add to the miniMatrix
                O_miniMatrix[index] = data(a, b);
              }
            }
          }
        }

        // Defining cell value
        double cellValue = data(i, j);

        // PEAK-5 SEARCH -----------------------------------------------------
        // Set a counter for slices
        int peak5 = 0;

        ////////////////////// N-S slice //////////////////////
        // Subset from outer mini matrix
        NumericVector tempVector = O_miniMatrix[ns_outer];
        tempVector = Rcpp::na_omit(tempVector);

        // Search for min & max
        double maxElement = max(tempVector);
        double minElement = min(tempVector);

        // If some min or max was found, increase the counter
        if((abs(cellValue - minElement) < 1e-6) || (abs(cellValue - maxElement) < 1e-6)){
          peak5++;
        }

        ////////////////////// W-E slice //////////////////////
        // Subset from outer mini matrix
        tempVector = O_miniMatrix[we_outer];
        tempVector = Rcpp::na_omit(tempVector);

        maxElement = max(tempVector);
        minElement = min(tempVector);

        // If some min or max was found, increase the counter
        if((abs(cellValue - minElement) < 1e-6) || (abs(cellValue - maxElement) < 1e-6)){
          peak5++;
        }

        ////////////////////// NW-SE slice //////////////////////
        // Subset from outer mini matrix
        tempVector = O_miniMatrix[nwse_outer];
        tempVector = Rcpp::na_omit(tempVector);

        maxElement = max(tempVector);
        minElement = min(tempVector);

        // If some min or max was found, increase the counter
        if((abs(cellValue - minElement) < 1e-6) || (abs(cellValue - maxElement) < 1e-6)){
          peak5++;
        }

        ////////////////////// NE-SW slice //////////////////////
        // Subset from outer mini matrix
        tempVector = O_miniMatrix[nesw_outer];
        tempVector = Rcpp::na_omit(tempVector);

        maxElement = max(tempVector);
        minElement = min(tempVector);

        // If some min or max was found, increase the counter
        if((abs(cellValue - minElement) < 1e-6) || (abs(cellValue - maxElement) < 1e-6)){
          peak5++;
        }

        // Defining tag of peak-5
        bool peak5_tag = peak5 == 4;


        // 2. Check for peaks and troughs within 1D 3-point slices through
        // sliding 3×3 window. The window slides west–east, north–south across
        // the matrix:

        // PEAK-3 SEARCH -----------------------------------------------------
        // Set a counter for slices
        int peak3 = 0;

        ////////////////////// N-S slice //////////////////////
        // Subset from outer mini matrix
        tempVector = I_miniMatrix[ns_inner];
        tempVector = Rcpp::na_omit(tempVector);

        // Search for min & max
        maxElement = max(tempVector);
        minElement = min(tempVector);

        // If some min or max was found, increase the counter
        if((abs(cellValue - minElement) < 1e-6) || (abs(cellValue - maxElement) < 1e-6)){
          peak3++;
        }

        ////////////////////// W-E slice //////////////////////
        // Subset from outer mini matrix
        tempVector = I_miniMatrix[we_inner];
        tempVector = Rcpp::na_omit(tempVector);

        maxElement = max(tempVector);
        minElement = min(tempVector);

        // If some min or max was found, increase the counter
        if((abs(cellValue - minElement) < 1e-6) || (abs(cellValue - maxElement) < 1e-6)){
          peak3++;
        }

        // Defining tag of peak-3
        bool peak3_tag = peak3 == 2;

        // 3. Apply the selective 2D 3×3 median filter within sliding 3×3 window.
        // If the window center is a significant 5-point extremum (Peak-5),
        // leave it intact (do not blunt it with median filter), otherwise if
        // the window center is a spike (Peak-3) use the 2D 3×3 median filter:

        if(!peak5_tag && peak3_tag){

          I_miniMatrix = Rcpp::na_omit(I_miniMatrix);

          arma::vec I_miniMatrix2 = as<arma::vec>(wrap(I_miniMatrix));

          arma::mat I_miniMatrix3 = arma::reshape(I_miniMatrix2, I_miniMatrix2.n_elem, 1);

          arma::mat quantVal = arma::median(I_miniMatrix3);

          // Put the median value of the Inner miniMatrix
          emptyData(i, j) = arma::conv_to < double >::from(quantVal);
        }else{
          emptyData(i, j) = cellValue;
        }
      }
    }
  }

  return emptyData;
}
