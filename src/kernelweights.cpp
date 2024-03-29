# include <RcppArmadillo.h>
# include <algorithm>
# include <string.h>

// [[Rcpp::depends(RcppArmadillo)]]

// The function `kernelweights` computes the kernel weights for the observations
// present in the data matrix `X_data` according to the kernel function `Kernel`
// and bandwidth `h`. `X_static` is a `n`-by-`p` data matrix consisting of `n`
// functional observations recorded at `p` points specified in the vector
// `t_vector`. The variable `x` is a row vector of length `p`, whose elements
// are the recorded values of the underlying function at the points in
// `t_vector`. The argument `Kernel` specifies the name of the kernel function,
// and `h` is a positive number which is the bandwidth of the kernel
// function. The implemented kernel functions are:
//    * `gaussian`: (1 / sqrt(2 pi)) exp(- u^2 / 2),
//    * `triangular`: (1 - abs(u)) (abs(u) <= 1),
//    * `epanechnikov`: (3/4) (1 - u^2) (abs(u) <= 1),
//    * `quartic`: (15/16) ((1 - u^2)^2) (abs(u) <= 1),
//    * `triweight`: (35/32) ((1 - u^2)^3) (abs(u) <= 1),
//    * `tricube`: (70/81) (1 - abs(u)^3)^3 (abs(u) <= 1),
//    * `uniform` (default): 0.5 (abs(u) <= 1).


// [[Rcpp::export]]
double kernelvalue(double u, Rcpp::String Kernel){
  if (Kernel == "gaussian"){
    return ((1 / sqrt(2 * arma::datum::pi)) * exp(- pow(u, 2) / 2));
  }else if (Kernel == "triangular"){
    return ((1 - abs(u)) * (double)(abs(u) <= 1));
  }else if (Kernel == "epanechnikov"){
    return ((3/4) * (1 - pow(u, 2)) * (double)(abs(u) <= 1));
  }else if (Kernel == "quartic"){
    return ((15/16) * pow((1 - pow(u, 2)), 2) * (double)(abs(u) <= 1));
  }else if (Kernel == "triweight"){
    return ((35/32) * pow((1 - pow(u, 2)), 3) * (double)(abs(u) <= 1));
  }else if (Kernel == "tricube"){
    return ((70/81) * pow((1 - pow(abs(u), 3)), 3) * (double)(abs(u) <= 1));
  }else if (Kernel == "uniform"){
    return (0.5 * (double)(abs(u) <= 1));
  }else{
    Rcpp::stop("Enter one of the implemented kernels: gaussian, triangular, epanechnikov, quartic, triweight, tricube, uniform.");
  }
}

// [[Rcpp::export]]
arma::vec kernelweights(arma::vec x, arma::mat X_static, arma::vec t_vector, double h, Rcpp::String Kernel){
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
    Weights[i] = kernelvalue((Distance_vector[i] / h), Kernel);

    sum_Weights = sum_Weights + Weights[i];
  }
  for (i = 0; i < n; ++i){
    Weights[i] = Weights[i] / sum_Weights;
  }

  return Weights;
}
