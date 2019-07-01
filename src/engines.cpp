#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <numeric>
#include <math.h>

using namespace Rcpp;
// using namespace std;

//' @importFrom Rcpp evalCpp
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
      int naSum = 0;

      // Multiply the value of each cell by the corresponding value of the kernel.
      for(int n = 0; n < knlcols; n++){
        for(int m = 0; m < knlrows; m++){
          int a = i + m - knlRowHalf;
          int b = j + n - knlColHalf;

          // If the value is a NA do not consider for sum and increase the naSum index
          if(std::isnan(data(a, b))){
            naSum++;
          }else{
            cumSum += data(a, b)*kernel(m, n);
          }
        }
      }

      // Assign sum of values to corresponding cell, if all the values were NA, result will be NA
      if(naSum < threshold){
        emptyData(i, j) = cumSum;
      }
    }
  }

  return emptyData;
}

// ENGINE 2: Convolution with quantiles
// [[Rcpp::export]]
NumericMatrix engine2(NumericMatrix data, NumericMatrix kernel, double probs, int naVal){

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

  int miniMatrixSize = knlrows*knlcols;
  NumericVector miniMatrix(miniMatrixSize);

  // Define output matrix, same dims of input
  NumericMatrix emptyData(nrows, ncols);
  std::fill(emptyData.begin(), emptyData.end(), NA_REAL);

  for(int j = knlColHalf; j < (ncols - knlColHalf); j++){
    for(int i = knlRowHalf; i < (nrows - knlRowHalf); i++){

      int naCounter = 0;

      // Multiply the value of each cell by the corresponding value of the kernel.
      for(int n = 0; n < knlcols; n++){
        for(int m = 0; m < knlrows; m++){
          int index = m*knlcols + n;
          int a = i + m - knlRowHalf;
          int b = j + n - knlColHalf;

          // If the value is a NA do not consider for sum and increase the naSum index
          if(abs(data(a, b) - naVal) > 0){
            miniMatrix[index] = data(a, b)*kernel(m, n);
          }else{
            miniMatrix[index] = naVal;
            naCounter++;
          }
        }
      }

      // Assign sum of values to corresponding cell, if all the values were NA, result will be NA
      if(naCounter < miniMatrixSize){
        // Sort values
        miniMatrix.sort();

        // Convert probs to an index
        int validPositions = (miniMatrixSize - naCounter);

        int index = (int)round(std::floor(validPositions*probs));

        double tempValue = miniMatrix[index];
        if(validPositions % 2 == 0){
          tempValue = (tempValue + miniMatrix[index - 1])/2;
        }

        // Put the median value of the Inner miniMatrix
        emptyData(i, j) = tempValue;
      }
    }
  }

  return emptyData;
}

// ENGINE 3: Mean filter
// [[Rcpp::export]]
NumericMatrix engine3(NumericMatrix data, int radius){
  int nrows = data.nrow();
  int ncols = data.ncol();

  NumericMatrix emptyData(nrows, ncols);
  std::fill(emptyData.begin(), emptyData.end(), NA_REAL);

  double halfRadiusDouble = std::floor(radius/2);
  int halfRadius = (int)round(halfRadiusDouble);

  for(int j = halfRadius; j < (ncols - halfRadius); j++){
    for(int i = halfRadius; i < (nrows - halfRadius); i++){

      double cumSum = 0;
      int k = 0;

      for(int n = 0; n < radius; n++){
        for(int m = 0; m < radius; m++){

          int a = i + m - halfRadius;
          int b = j + n - halfRadius;

          if(!std::isnan(data(a, b))){
            cumSum += data(a, b);
            k++;
          }
        }
      }

      if(k > 0){
        emptyData(i, j) = cumSum/k;
      }
    }
  }

  return emptyData;
}

// ENGINE 4: Quantile filter
// [[Rcpp::export]]
NumericMatrix engine4(NumericMatrix data, int radius, double probs, int naVal){
  int nrows = data.nrow();
  int ncols = data.ncol();

  NumericMatrix emptyData(nrows, ncols);
  std::fill(emptyData.begin(), emptyData.end(), NA_REAL);

  int miniMatrixSize = radius*radius;
  NumericVector miniMatrix(miniMatrixSize);

  double halfRadiusDouble = std::floor(radius/2);
  int halfRadius = (int)round(halfRadiusDouble);

  for(int j = halfRadius; j < (ncols - halfRadius); j++){
    for(int i = halfRadius; i < (nrows - halfRadius); i++){

      int naCounter = 0;

      for(int n = 0; n < radius; n++){
        for(int m = 0; m < radius; m++){
          int index = m*radius + n;
          int a = i + m - halfRadius;
          int b = j + n - halfRadius;

          if(abs(data(a, b) - naVal) > 0){
            miniMatrix[index] = data(a, b);
          }else{
            miniMatrix[index] = naVal;
            naCounter++;
          }
        }
      }

      // Assign sum of values to corresponding cell, if all the values were NA, result will be NA
      if(naCounter < miniMatrixSize){
        // Sort values
        miniMatrix.sort();

        // Convert probs to an index
        int validPositions = (miniMatrixSize - naCounter);

        int index = (int)round(std::floor(validPositions*probs));

        double tempValue = miniMatrix[index];
        if(validPositions % 2 == 0){
          tempValue = (tempValue + miniMatrix[index - 1])/2;
        }

        // Put the median value of the Inner miniMatrix
        emptyData(i, j) = tempValue;
      }
    }
  }

  return emptyData;
  // return miniMatrix;
}

// ENGINE 5: Contextual Median Filter
// Proposed by Belkin et al. (2009), doi:10.1016/j.jmarsys.2008.11.018
// [[Rcpp::export]]
NumericMatrix engine5(NumericMatrix data, double probs, int I_radius, int O_radius, int naVal){
  int nrows = data.nrow();
  int ncols = data.ncol();

  // Create a NA matrix for output
  NumericMatrix emptyData(nrows, ncols);
  std::fill(emptyData.begin(), emptyData.end(), NA_REAL);

  double I_halfRadiusDouble = std::floor(I_radius/2);
  int I_halfRadius = (int)round(I_halfRadiusDouble);

  double O_halfRadiusDouble = std::floor(O_radius/2);
  int O_halfRadius = (int)round(O_halfRadiusDouble);

  int O_minimatrixSize = O_radius*O_radius;
  int I_minimatrixSize = I_radius*I_radius;

  // OUTER INDEX
  // N-S split outer matrix
  IntegerVector ns_outer(O_radius);
  int index = 0;
  for(int i = O_halfRadius; i < O_minimatrixSize; i = i + O_radius){
    ns_outer(index) = i;
    index++;
  }

  // W-E split outer matrix
  IntegerVector we_outer(O_radius);
  index = 0;
  for(int i = O_radius*O_halfRadius; i < O_radius*(O_halfRadius + 1); i++){
    we_outer(index) = i;
    index++;
  }

  // NW-SE split outer matrix
  IntegerVector nwse_outer(O_radius);
  index = 0;
  for(int i = 0; i < O_minimatrixSize; i = i + O_radius + 1){
    nwse_outer(index) = i;
    index++;
  }

  // NE-SW split outer matrix
  IntegerVector nesw_outer(O_radius);
  index = 0;
  for(int i = O_radius - 1; i < (O_minimatrixSize - 1); i = i + O_radius - 1){
    nesw_outer(index) = i;
    index++;
  }


  // INNER INDEX
  // N-S split inner matrix
  IntegerVector ns_inner(I_radius);
  index = 0;
  for(int i = I_halfRadius; i < I_minimatrixSize; i = i + I_radius){
    ns_inner(index) = i;
    index++;
  }

  // W-E split inner matrix
  IntegerVector we_inner(I_radius);
  index = 0;
  for(int i = I_radius*I_halfRadius; i < I_radius*(I_halfRadius + 1); i++){
    we_inner(index) = i;
    index++;
  }

  // NW-SE split inner matrix
  IntegerVector nwse_inner(I_radius);
  index = 0;
  for(int i = 0; i < I_minimatrixSize; i = i + I_radius + 1){
    nwse_inner(index) = i;
    index++;
  }

  // NE-SW split inner matrix
  IntegerVector nesw_inner(I_radius);
  index = 0;
  for(int i = I_radius - 1; i < (I_minimatrixSize - 1); i = i + I_radius - 1){
    nesw_inner(index) = i;
    index++;
  }

  for(int j = O_halfRadius; j < (ncols - O_halfRadius); j++){
    for(int i = O_halfRadius; i < (nrows - O_halfRadius); i++){

      // Create outputs for Outer and Inner mini matrix
      NumericVector O_miniMatrix(O_minimatrixSize);
      NumericVector I_miniMatrix(I_minimatrixSize);

      // Build outer minimatrix (checking NA)
      int O_naCounter = 0;
      for(int n = 0; n < O_radius; n++){
        for(int m = 0; m < O_radius; m++){
          int index = m*O_radius + n;
          int a = i + m - O_halfRadius;
          int b = j + n - O_halfRadius;

          if(abs(data(a, b) - naVal) > 0){
            O_miniMatrix[index] = data(a, b);
          }else{
            O_miniMatrix[index] = naVal;
            O_naCounter++;
          }
        }
      }

      // Build outer minimatrix (checking NA)
      int I_naCounter = 0;
      for(int n = 0; n < I_radius; n++){
        for(int m = 0; m < I_radius; m++){
          int index = m*I_radius + n;
          int a = i + m - I_halfRadius;
          int b = j + n - I_halfRadius;

          if(abs(data(a, b) - naVal) > 0){
            I_miniMatrix[index] = data(a, b);
          }else{
            I_miniMatrix[index] = naVal;
            I_naCounter++;
          }
        }
      }

      if((O_naCounter < O_minimatrixSize) && (I_naCounter < I_minimatrixSize)){

        // Set a counter for outer slices
        int outerTag = 0;

        ////////////////////// N-S slice //////////////////////
        // Subset from outer mini matrix
        NumericVector tempVector = O_miniMatrix[ns_outer];
        tempVector = tempVector[tempVector < naVal];

        // Search for min & max
        int maxElement = which_max(tempVector);
        int minElement = which_min(tempVector);

        // If some min or max was found, increase the counter
        if((maxElement == O_halfRadius) | (minElement == O_halfRadius)){
          outerTag++;
        }

        ////////////////////// W-E slice //////////////////////
        // Subset from outer mini matrix
        tempVector = O_miniMatrix[we_outer];
        tempVector = tempVector[tempVector < naVal];

        // Search for min & max
        maxElement = which_max(tempVector);
        minElement = which_min(tempVector);

        // If some min or max was found, increase the counter
        if((maxElement == O_halfRadius) | (minElement == O_halfRadius)){
          outerTag++;
        }

        ////////////////////// NW-SE slice //////////////////////
        // Subset from outer mini matrix
        tempVector = O_miniMatrix[nwse_outer];
        tempVector = tempVector[tempVector < naVal];

        // Search for min & max
        maxElement = which_max(tempVector);
        minElement = which_min(tempVector);

        // If some min or max was found, increase the counter
        if((maxElement == O_halfRadius) | (minElement == O_halfRadius)){
          outerTag++;
        }

        ////////////////////// NE-SW slice //////////////////////
        // Subset from outer mini matrix
        tempVector = O_miniMatrix[nesw_outer];
        tempVector = tempVector[tempVector < naVal];

        // Search for min & max
        maxElement = which_max(tempVector);
        minElement = which_min(tempVector);

        // If some min or max was found, increase the counter
        if((maxElement == O_halfRadius) | (minElement == O_halfRadius)){
          outerTag++;
        }

        // Check if data(i, j) has been a max or min along 4 slices for outer kernel
        if(outerTag == 4){
          emptyData(i, j) = data(i, j);
        }else{
          // Set a counter for inner slices
          int innerTag = 0;

          ////////////////////// N-S slice //////////////////////
          // Subset from inner mini matrix
          NumericVector tempVector = I_miniMatrix[ns_inner];
          tempVector = tempVector[tempVector < naVal];

          // Search for min & max
          maxElement = which_max(tempVector);
          minElement = which_min(tempVector);

          // If some min or max was found, increase the counter
          if((maxElement == I_halfRadius) | (minElement == I_halfRadius)){
            innerTag++;
          }

          ////////////////////// W-E slice //////////////////////
          // Subset from inner mini matrix
          tempVector = I_miniMatrix[we_inner];
          tempVector = tempVector[tempVector < naVal];

          // Search for min & max
          maxElement = which_max(tempVector);
          minElement = which_min(tempVector);

          // If some min or max was found, increase the counter
          if((maxElement == I_halfRadius) | (minElement == I_halfRadius)){
            innerTag++;
          }

          ////////////////////// NW-SE slice //////////////////////
          // Subset from inner mini matrix
          tempVector = I_miniMatrix[nwse_inner];
          tempVector = tempVector[tempVector < naVal];

          // Search for min & max
          maxElement = which_max(tempVector);
          minElement = which_min(tempVector);

          // If some min or max was found, increase the counter
          if((maxElement == I_halfRadius) | (minElement == I_halfRadius)){
            innerTag++;
          }

          ////////////////////// NE-SW slice //////////////////////
          // Subset from inner mini matrix
          tempVector = I_miniMatrix[nesw_inner];
          tempVector = tempVector[tempVector < naVal];

          // Search for min & max
          maxElement = which_max(tempVector);
          minElement = which_min(tempVector);

          // If some min or max was found, increase the counter
          if((maxElement == I_halfRadius) | (minElement == I_halfRadius)){
            innerTag++;
          }

          // Check if data(i, j) has been a max or min along 4 slices for inner kernel
          if(innerTag == 4){
            // Sort values
            I_miniMatrix.sort();

            // Convert probs to an index
            int validPositions = (I_minimatrixSize - I_naCounter);

            int index = (int)round(std::floor(validPositions*probs));

            double tempValue = I_miniMatrix[index];
            if(validPositions % 2 == 0){
              tempValue = (tempValue + I_miniMatrix[index - 1])/2;
            }

            // Put the median value of the Inner miniMatrix
            emptyData(i, j) = tempValue;
          }else{
            emptyData(i, j) = data(i, j);
          }
        }
      }
    }
  }

  return emptyData;
}
