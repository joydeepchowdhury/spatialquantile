# include <RcppArmadillo.h>
# include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

// The function `spquantile` computes the weighted spatial quantile from the
// data matrix `X_data_original` with corresponding weights in the variable
// `Weights`. The data matrix `X_data_original` is of dimension `n`-by-`p` and
// contains `n` functional observations which are recorded on a grid of length
// `p`, given in the row vector `t_vector`. The variable `Weights` is a column
// vector of length `n`, whose entries are the weights of the corresponding
// observations in `X_data_original`. The variables `u_index` and `c` are a
// positive integer and a real number in (-1, 1) respectively, which together
// determine the variable `u` for the spatial quantile computation.
// For example, if `u_index` = 1 and `c` = 0.5, then we compute the weighted
// spatial quantile corresponding to `u` equaling to 0.5 times the weighted
// first principal component of the data in `X_static_original`.


// [[Rcpp::export]]
arma::mat spatialquantileconfidenceset(arma::mat Data_original, arma::vec Weights, int u_index, double c, arma::vec t_vector, double alpha){
  
  return ConfidenceSet;
}
