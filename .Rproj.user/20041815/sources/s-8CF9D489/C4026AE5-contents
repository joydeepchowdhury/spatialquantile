# include <RcppArmadillo.h>
# include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec spquantile(arma::mat Data, arma::vec Weights, int u_index, double c, arma::vec t_vector){
  int n = Data.nrow();
  int p = Data.ncol();
  arma::vec Quantile(p);
  
  int i, j, k;
  
  // Checking whether there is only one distinct observation in Data_original
  arma::vec z(p);
  for (i = 0; i < p; ++i){
    z[i] = Data(1, i);
  }
  
  arma::mat Difference(n, p);
  for (i = 0; i < n; ++i){
    for (j = 0; j < p; ++j){
      Difference(i, j) = z[j] - Data(i, j);
    }
  }
  
  arma::vec norm_difference(n);
  double temp1;
  for (i = 0; i < n; ++i){
    temp1 = 0;
    for (j = 0; j < (p - 1); ++j){
      // Integration computed using trapezoidal rule
      temp1 = temp1 + ((pow(Difference(i, j + 1), 2) + pow(Difference(i, j), 2)) / 2) * (t_vector[j + 1] - t_vector[j]);
    }
    norm_difference[i] = sqrt(temp1);
  }
  
  double sum_norm_difference = 0;
  for (i = 0; i < n; ++i){
    sum_norm_difference = sum_norm_difference + norm_difference[i];
  }
  if (sum_norm_difference == 0){
    for (j = 0; j < p; ++j){
      Quantile[j] = z[j];
    }
    return Quantile;
  }
  
  z = Data_original(1,:);
  Difference = ones(size(Data_original,1),1) * z - Data_original;
  norm_Difference = sqrt(trapz(t_vector, Difference.^2, 2));
  if (sum(norm_Difference) == 0){
    Quantile = Data_original;
    return
  }
    
  // Dimension reduction for the input data
      
  n = size(Data_original,1);
  if sum(Weights) ~= 1
  Weights = Weights / sum(Weights);
  end
}
