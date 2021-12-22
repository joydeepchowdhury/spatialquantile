# include <RcppArmadillo.h>
# include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
uint64_t choose_rcpp(uint64_t n, uint64_t k){
  if(k == 0) return 1;
  return (n * choose_rcpp(n - 1, k - 1)) / k;
}
