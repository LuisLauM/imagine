#include <Rcpp.h>
#include <algorithm>
#include <math.h>
#include <fstream>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib imagine, .registration = TRUE

// is_extreme function
// Evaluate if the center of a defined matrix is a maximum/minimum in the whole
// matrix
int is_extreme(NumericMatrix in_mat, int direction){

  int side = in_mat.nrow();

  NumericVector tempVector(side ^ 2);
  tempVector = tempVector * NA_REAL;
  int midPos = std::floor(side/2);

  double cellVal = in_mat(midPos, midPos);

  switch(direction){
  // WE slice
  case 1:
    tempVector = in_mat(midPos, _);
    break;

    // NS slice
  case 2:
    tempVector = in_mat(_, midPos);
    break;

    // NW-SE slice
  case 3:
    for(int i = 0; i < side; i++){
      tempVector[i] = in_mat(i, i);
    }
    break;

    // NE-SW slice
  case 4:
    for(int i = 0; i < side; i++){
      tempVector[i] = in_mat(i, side - i);
    }
    break;
  default:
    break;
  }

  // Remove NA of tempVector
  tempVector = na_omit(tempVector);

  // Search for min & max
  double maxElement = max(tempVector);
  double minElement = min(tempVector);

  // If some min or max was found, increase the counter
  int out = (cellVal <= minElement) || (cellVal >= maxElement);

  return out;
}

// ENGINE 5: Contextual Median Filter
// Proposed by Belkin et al. (2009), doi:10.1016/j.jmarsys.2008.11.018
// [[Rcpp::export]]
NumericMatrix engine5_CMF(NumericMatrix data, int i_size = 3, int o_size = 5){

  // Get dimension of input matrix
  int nrows = data.nrow();
  int ncols = data.ncol();

  // Create a NA matrix for output
  NumericMatrix emptyData(nrows, ncols);
  std::fill(emptyData.begin(), emptyData.end(), NA_REAL);

  int i_halfsize = floor(i_size/2);
  int o_halfsize = floor(o_size/2);

  // ---------------------------------------------------------------------------
  // 1. Check for peaks and troughs within 1D 5-point slices through a sliding
  // 5 x 5 window. The window slides east-west (E-W), northsouth (N-S) across
  // the matrix:

  // Loop for the whole matrix
  for(int j = o_halfsize; j < (ncols - o_halfsize); j++){

    for(int i = o_halfsize; i < (nrows - o_halfsize); i++){

      // Defining cell value
      double cellValue = data(i, j);

      // Create outputs for Outer and Inner mini matrix
      NumericMatrix O_miniMatrix = data(Range(i - o_halfsize, i + o_halfsize),
                                        Range(j - o_halfsize, j + o_halfsize));
      NumericMatrix I_miniMatrix = data(Range(i - i_halfsize, i + i_halfsize),
                                        Range(j - i_halfsize, j + i_halfsize));

      // If some of the inner matrices are full of NA, pass to the next cell
      if((all(is_na(O_miniMatrix)) + all(is_na(I_miniMatrix))) > 0) continue;

      // If the cell value IS NOT A NA, evaluating if it is an extreme value for
      // a 3x3 and 5x5 matrix neightborhood
      if(!NumericVector::is_na(cellValue)){

        // PEAK-5 SEARCH
        // Set a counter for slices
        int peak5 = 0;

        ////////////////////// W-E slice //////////////////////
        peak5 += is_extreme(O_miniMatrix, 1);

        ////////////////////// N-S slice //////////////////////
        peak5 += is_extreme(O_miniMatrix, 2);

        ////////////////////// NW-SE slice //////////////////////
        peak5 += is_extreme(O_miniMatrix, 3);

        ////////////////////// NE-SW slice //////////////////////
        peak5 += is_extreme(O_miniMatrix, 4);

        // Defining tag of peak-5
        bool peak5_tag = (peak5 == 4);


        // ---------------------------------------------------------------------
        // 2. Check for peaks and troughs within 1D 3-point slices through
        // sliding 3x3 window. The window slides west-east, north-south across
        // the matrix:

        // PEAK-3 SEARCH
        // Set a counter for slices
        int peak3 = 0;

        ////////////////////// W-E slice //////////////////////
        peak3 += is_extreme(I_miniMatrix, 1);

        ////////////////////// N-S slice //////////////////////
        peak3 += is_extreme(I_miniMatrix, 2);

        // Defining tag of peak-5
        bool peak3_tag = (peak3 == 2);


        // ---------------------------------------------------------------------
        // 3. Apply the selective 2D 3x3 median filter within sliding 3x3 window.
        // If the window center is a significant 5-point extremum (Peak-5),
        // leave it intact (do not blunt it with median filter), otherwise if
        // the window center is a spike (Peak-3) use the 2D 3x3 median filter:
        if(!peak5_tag && peak3_tag){

          // Calculate median of inner matrix
          NumericVector outVal = na_omit(as<NumericVector>(I_miniMatrix));

          // Replace the corresponding position in the output matrix
          emptyData(i, j) = Rcpp::median(outVal);
        }else{
          emptyData(i, j) = cellValue;
        }
      }
    }
  }

  return emptyData;
}

