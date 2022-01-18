# include <RcppArmadillo.h>
# include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]









function rankings = wsdrank(X_to_rank, X_data, X_data_weights, t_vector)
  
  % This function calculates the ranks the observations in X_to_rank
  % according to the decreasing value of the weighted spatial depth with
  % respect to the data matrix X_data and corresponding weights in the column
  % vector X_data_weights. The dimension of the data matrix X_data is n-by-p,
  % where n is the sample size and p is the dimension of an observation
  % vector. The dimension of the matrix X_to_rank is m-by-p, where m is the
  % number of observations to rank. X_data_weights is a n-by-1 column vector,
  % and t_vector is a row vector recording the grid on which the values of
  % the functional observations underlying X_data and X_to_rank are recorded.
  % The number p must be greater than 1.