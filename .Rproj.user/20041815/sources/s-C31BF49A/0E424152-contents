# include <RcppArmadillo.h>
# include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

// This function computes the kernel weights for the observation present in
// the data matrix `X_data` according to the kernel function `Kernel` and
// bandwidth `h`. `X_static` is a n-by-p data matrix consisting of `n` functional
// observations recorded at `p` points specified in the vector `t_vector`. The
// variable `x` is a row vector of length `p`, whose elements are the recorded
// values of the underlying function at the points in `t_vector`. The argument
// named `Kernel` is a function handle specifying the kernel function, and `h`
// is a positive number which is the bandwidth of the kernel function.

// [[Rcpp::export]]
arma::vec kernelweights(arma::vec x, arma::mat X_static, arma::vec t_vector, double h, char Kernel){
  int n = X_static.n_rows;
  int p = X_static.n_cols;
  
  arma::vec Weights(n);
  
  arma::vec Distance_vector(n);
  
  arma::vec norm_difference(n);
  double temp1;
  for (i = 0; i < n; ++i){
    temp1 = 0;
    for (j = 0; j < (p - 1); ++j){
      // Integration computed using trapezoidal rule
      temp1 = temp1 + ((pow(Difference(i, j + 1), 2) + pow(Difference(i, j), 2)) / 2) * (t_vector[j + 1] - t_vector[j]);
    }
    Distance_vector[i] = sqrt(temp1);
  }
  
  
  
  
  Distance_vector = trapz(t_vector, (ones(size(X_static,1),1) * x - X_static).^2, 2);
  
  Weights = Kernel(Distance_vector / h);
  Weights = Weights / sum(Weights);
  
  
  
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
