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
  
  int i, j;
  
  arma::vec Weights(n);
  
  arma::mat Difference(n, p);
  for (i = 0; i < n; ++i){
    for (j = 0; j < p; ++j){
      Difference(i, j) = x[j] - X_static(i, j);
    }
  }
  
  arma::vec Distance_vector(n);
  double temp1;
  for (i = 0; i < n; ++i){
    temp1 = 0;
    for (j = 0; j < (p - 1); ++j){
      // Integration using trapezoidal rule
      temp1 = temp1 + ((pow(Difference(i, j + 1), 2) + pow(Difference(i, j), 2)) / 2) * (t_vector[j + 1] - t_vector[j]);
    }
    Distance_vector[i] = sqrt(temp1);
  }
  
  double sum_Weights = 0;
  for (i = 0; i < n; ++i){
    Weights[i] = kernel((Distance_vector[i] / h), Kernel);
    
    sum_Weights = sum_Weights + Weights[i];
  }
  for (i = 0; i < n; ++i){
    Weights[i] = Weights[i] / sum_Weights;
  }
  
  return Weights;
}

// [[Rcpp::export]]
double kernel(double u, char Kernel){
  if (Kernel == 'g'){
    return ((1 / sqrt(2 * arma::datum::pi)) * exp(- pow(u, 2) / 2));
  })
}
