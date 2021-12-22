# include <RcppArmadillo.h>
# include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec spquantile(arma::mat Data, arma::vec Weights, int u_index, double c, arma::vec t_vector){
  // Checking whether there is only one distinct observation in Data_original
  
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
