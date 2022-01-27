#ifndef __KERNELWEIGHTS__
#define __KERNELWEIGHTS__

arma::vec kernelweights(arma::vec x, arma::mat X_static, arma::vec t_vector, double h, Rcpp::String Kernel);

#endif