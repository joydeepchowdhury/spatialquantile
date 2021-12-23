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
  
  // z = Data_original(1,:);
  // Difference = ones(size(Data_original,1),1) * z - Data_original;
  // norm_Difference = sqrt(trapz(t_vector, Difference.^2, 2));
  // if (sum(norm_Difference) == 0){
  //   Quantile = Data_original;
  //   return
  // }
    
  // Dimension reduction for the input data
      
  double sum_Weights = 0;
  for (i = 0; i < n; ++i){
    sum_Weights = sum_Weights + Weights[i];
  }
  if (sum_Weights != 1){
    for (i = 0; i < n; ++i){
      Weights[i] = Weights[i] / sum_Weights;
    }
  }
  
  double t_1 = sqrt(n);
  double t_2 = 2 * pow(n, 1/3);
  double t_3;
  if (t_1 < t_2){
    t_3 = t_1;
  }else{
    t_3 = t_2;
  }
  double d_n = floor(t_3);
  
  Weighted_Mean = mean((Weights * ones(1,size(Data_original,2))) .* Data_original, 1);
  Centred_Data = Data_original - ones(n,1) * Weighted_Mean;
  Weighted_Cov_Matrix = Centred_Data' * diag(Weights) * Centred_Data;
  Weighted_Cov_Matrix = (Weighted_Cov_Matrix + Weighted_Cov_Matrix') / 2;
  [Eigenvectors, Eigenvalues] = eig(Weighted_Cov_Matrix);
  vector_Eigenvalues = diag(Eigenvalues);
  [~, index_Eigenvalues_sorted] = sort(vector_Eigenvalues,'descend');
  Eigenvectors_sorted = Eigenvectors(:,index_Eigenvalues_sorted);
  Coefficient_Matrix = Centred_Data * Eigenvectors_sorted;
  Eigenvectors_sorted_truncated = Eigenvectors_sorted(:,1:d_n);
  Coefficient_Matrix_truncated = Coefficient_Matrix(:,1:d_n);
  Data_reduced = Coefficient_Matrix_truncated;
  
  Data = Data_reduced;
}
