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
  
  Data = Data_reduced;
  
  arma::vec u(d_n);
  for (j = 0; j < d_n; ++j){
    u[j] = 0;
  }
  if (u_index > 0 && u_index <= d_n){
    u[u_index - 1] = c;
  }
  if (u_index > d_n){
    u[d_n - 1] = c;
  }          // Look here if error handling can be implemented
  
  // Checking whether the weighted quantile is present in the data itself
  
  int Check = 0;
  arma::vec x(p);
  arma::mat U(n, p);
  arma::vec weights_for_i(n), norms_U(n), weighted_norms_U(n);
  for (i = 0; i < n; ++i){
    for (j = 0; j < p; ++j){
      x[j] = Data(i, j);
    }
    
    for (k = 0; k < n; ++k){
      norms_U[k] = 0;
      for (j = 0; j < p; ++j){
        U(k, j) = Data(k, j) - x[j];
        
        norms_U[k] = norms_U[k] + pow(U(k, j), 2);
      }
      norms_U[k] = sqrt(norms_U[k]);
      
      weights_for_i[k] = Weights[k];
      
      weighted_norms_U[k] = weights_for_i[k] * norms_U[k];
    }
    
    
    
    
  }
  for i=1:n
    X = Data;
  x = X(i,:);
  U = X - ones(size(X,1),1) * x;
  weights_for_i = Weights;
  
  weighted_norms_U = weights_for_i .* sqrt(sum(U.^2,2));
  all_indices_for_i = 1:n;
  J_i = all_indices_for_i(weighted_norms_U == 0);
  J_i_complement = setdiff(all_indices_for_i, J_i);
  J_i = setdiff(J_i, i);
  
  U_new = U(J_i_complement,:);
  U_new = U_new ./ ( sqrt(sum(U_new.^2,2)) * ones(1,size(U_new,2)) );
  weights_proper = weights_for_i(J_i_complement);
  V = sum((weights_proper * ones(1,size(U_new,2))) .* U_new, 1) +...
    sum(weights_proper) * u;
  
  if sqrt(sum(V.^2)) <= (1 + sqrt(sum(u.^2))) * (sum(weights_for_i(J_i)))
    Quantile_coefficients = x;
  Check = 1;
  break
    end
    end
}
