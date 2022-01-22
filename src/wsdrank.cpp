# include <RcppArmadillo.h>
# include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

// This function calculates the ranks the observations in `X_to_rank` according
// to the decreasing value of the weighted spatial depth with respect to the
// data matrix `X_data` and corresponding weights in the column vector
// `X_data_weights`. The dimension of the data matrix `X_data` is `n`-by-`p`,
// where `n` is the sample size and `p` is the dimension of an observation
// vector. The dimension of the matrix `X_to_rank` is `m`-by-`p`, where `m` is
// the number of observations to rank. `X_data_weights` is a `n`-by-1 column
// vector, and `t_vector` is a row vector recording the grid on which the values
// of the functional observations underlying `X_data` and `X_to_rank` are
// recorded. The number `p` must be greater than 1.


// [[Rcpp::export]]
arma::vec wsdrank(arma::mat X_to_rank, arma::mat X_data, X_data_weights, t_vector){
  int number_of_points = X_to_rank.n_rows;
  int n = X_data.n_rows;
  int p = X_data.n_cols;
  
  int i, j, k;
  double temp;
  
  arma::vec wsd(number_of_points);
  arma::vec y(p);
  arma::mat difference(n, p);
  arma::vec norm_difference(n);
  double num_nonzero_norm;
  int index;
  for (i = 0; i < number_of_points; ++i){
    for (j = 0; j < p; ++j){
      y[j] = X_to_rank(i, j);
    }
    
    num_nonzero_norm = 0;
    for (k = 0; k < n; ++k){
      temp = 0;
      
      for (j = 0; j < p; ++j){
        difference(k, j) = y[j] - X_data(k, j);
      }
      
      for (j = 0; j < (p - 1); ++j){
        // Integration computed using trapezoidal rule
        temp = temp + ((pow(difference(k, j + 1), 2) + pow(difference(k, j), 2)) / 2) * (t_vector[j + 1] - t_vector[j]);
      }
      norm_difference[k] = sqrt(temp);
      
      if (norm_difference[k] != 0){
        num_nonzero_norm = num_nonzero_norm + 1;
      }
    }
    
    if (num_nonzero_norm == 0){
      arma::mat weighted_average(n, p);
      for (k = 0; k < n; ++k){
        for (j = 0; j < p; ++j){
          weighted_average(k, j) = 0;
        }
      }
      
      // NEED TO DEBUG HERE!!!!
      
      
      
    }else{
      arma::mat difference_proper(num_nonzero_norm, p);
      arma::vec norm_difference_proper(num_nonzero_norm);
      arma::vec weights_proper(num_nonzero_norm);
      index = 0;
      for (k = 0; k < n; ++k){
        if (norm_difference[k] != 0){
          for (j = 0; j < p; ++j){
            difference_proper(index, j) = difference(k, j);
          }
          
          norm_difference_proper[index] = norm_difference[k];
          
          weights_proper[index] = X_data_weights[k];
          
          index = index + 1;
        }
      }
      
      arma::mat scaled_difference_proper(num_nonzero_norm, p);
      for (k = 0; k < num_nonzero_norm; ++k){
        for (j = 0; j < p; ++j){
          scaled_difference_proper(k, j) = difference_proper(k, j) / norm_difference_proper[k];
        }
      }
      
      
      
      
    }
    
    
    
  }
  
}



function rankings = wsdrank(X_to_rank, X_data, X_data_weights, t_vector)
  
  
difference = (ones(size(X_data,1),1) * y) - X_data;
norm_difference = sqrt(trapz(t_vector, difference.^2, 2));
check_nonzero_norm = (norm_difference ~= 0);
if sum(check_nonzero_norm) == 0
weighted_average = zeros(size(difference));
else
  difference_proper = difference(check_nonzero_norm,:);
norm_difference_proper = norm_difference(check_nonzero_norm,:);
weights_proper = X_data_weights(check_nonzero_norm);
scaled_difference_proper = difference_proper ./ ...
  ( norm_difference_proper * ones(1,size(difference_proper,2)) );
scaled_difference_proper_weighted = ( weights_proper * ones(1,size(difference_proper,2)) )...
                                                                                          .* scaled_difference_proper;
weighted_average = sum(scaled_difference_proper_weighted,1) / sum(X_data_weights);
end
  
  weighted_spatial_depth_y = 1 - sqrt(trapz(t_vector, weighted_average.^2));
wsd(i) = weighted_spatial_depth_y;
end
  
  [~,rankings] = sort(wsd, 'descend');
end