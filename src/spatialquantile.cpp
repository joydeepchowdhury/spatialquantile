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
  
  arma::vec Weighted_Mean(p);
  for (j = 0; j < p; ++j){
    Weighted_Mean[j] = 0;
    for (i = 0; i < n; ++i){
      Weighted_Mean[j] = Weighted_Mean[j] + (Weights[i] * Data(i, j));
    }
  }
  arma::mat Centered_Data(n, p);
  for (i = 0; i < n; ++i){
    for (j = 0; j < p; ++j){
      Centered_Data(i, j) = Data(i, j) - Weighted_Mean[j];
    }
  }
  arma::mat Weighted_Cov_Matrix(p, p);
  for (j = 0; j < p; ++j){
    for (k = 0; k < p; ++k){
      Weighted_Cov_Matrix(j, k) = 0;
      for (i = 0; i < n; ++i){
        Weighted_Cov_Matrix(j, k) = Weighted_Cov_Matrix(j, k) + (Weights[i] * Centered_Data(i, j) * Centered_Data(i, k));
      }
    }
  }
  arma::vec Eigenvalues;
  arma::mat Eigenvectors;
  eig_sym(Eigenvalues, Eigenvectors, Weighted_Cov_Matrix);
  arma::uvec index_Eigenvalues_sorted = arma::sort_index(Eigenvalues, "descend");
  arma::mat Eigenvectors_sorted(p, p);
  for (k = 0; k < p; ++k){
    for (j = 0; j < p; ++j){
      Eigenvectors_sorted(j, k) = Eigenvectors(j, index_Eigenvalues_sorted[k]);
    }
  }
  
  arma::mat Eigenvectors_sorted_truncated(p, d_n);
  for (k = 0; k < d_n; ++k){
    for (j = 0; j < p; ++k){
      Eigenvectors_sorted_truncated(j, k) = Eigenvectors_sorted(j, k);
    }
  }
  
  arma::mat Coefficient_Matrix_truncated(n, d_n);
  for (i = 0; i < n; ++i){
    for (j = 0; j < d_n; ++j){
      Coefficient_Matrix_truncated(i, j) = 0;
      for (k = 0; k < p; ++k){
        Coefficient_Matrix_truncated(i, j) = Coefficient_Matrix_truncated(i, j) +
          (Centered_Data(i, k) * Eigenvectors_sorted_truncated(k, j));
      }
    }
  }
  
  arma::mat Data_reduced = Coefficient_Matrix_truncated;
  
  
  
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
