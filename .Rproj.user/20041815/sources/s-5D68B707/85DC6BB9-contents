# include <RcppArmadillo.h>
# include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]

// The function `spquantile` computes the weighted spatial quantile from the
// data matrix `X_data_original` with corresponding weights in the variable
// `Weights`. The data matrix `X_data_original` is of dimension `n`-by-`p` and
// contains `n` functional observations which are recorded on a grid of length
// `p`, given in the row vector `t_vector`. The variable `Weights` is a column
// vector of length `n`, whose entries are the weights of the corresponding
// observations in `X_data_original`. The variables `u_index` and `c` are a
// positive integer and a real number in (-1, 1) respectively, which together
// determine the variable `u` for the spatial quantile computation.
// For example, if `u_index` = 1 and `c` = 0.5, then we compute the weighted
// spatial quantile corresponding to `u` equaling to 0.5 times the weighted
// first principal component of the data in `X_static_original`.


// [[Rcpp::export]]
arma::mat spatialquantileconfidenceset(arma::mat Data, arma::vec Weights,
                                       int u_index, double c, arma::vec t_vector, double alpha){
  
  int n = Data.n_rows;
  int p = Data.n_cols;
  arma::mat ConfidenceSet(2, p);
  
  int i, j, k, l;
  
  // Checking whether there is only one distinct observation in Data
  arma::vec z(p);
  for (i = 0; i < p; ++i){
    z[i] = Data(0, i);
  }
  
  arma::mat Difference(n, p);
  for (i = 0; i < n; ++i){
    for (j = 0; j < p; ++j){
      Difference(i, j) = z[j] - Data(i, j);
    }
  }
  
  arma::vec norm_sq_difference(n);
  double temp1;
  for (i = 0; i < n; ++i){
    temp1 = 0;
    for (j = 0; j < (p - 1); ++j){
      // Integration computed using trapezoidal rule
      temp1 = temp1 + ((pow(Difference(i, j + 1), 2) + pow(Difference(i, j), 2)) / 2) * (t_vector[j + 1] - t_vector[j]);
    }
    norm_sq_difference[i] = temp1;
  }
  
  double sum_norm_sq_difference = 0;
  for (i = 0; i < n; ++i){
    sum_norm_sq_difference = sum_norm_sq_difference + norm_sq_difference[i];
  }
  if (sum_norm_sq_difference == 0){
    for (j = 0; j < p; ++j){
      ConfidenceSet(0, j) = z[j];
      ConfidenceSet(1, j) = z[j];
    }
    return ConfidenceSet;
  }
  
  // Main computation begins
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
  double t_3 = min(t_1, t_2);
  double d_n = floor(t_3);
  
  arma::vec Weighted_Mean;
  for (j = 0; j < p; ++j){
    Weighted_Mean[j] = 0;
    for (i = 0; i < n; ++i){
      Weighted_Mean[j] = Weighted_Mean[j] + (Weights[i] * Data(i, j));
    }
  }
  
  arma::mat Centred_Data(n, p);
  for (i = 0; i < n; ++i){
    for (j = 0; j < p; ++j){
      Centred_Data(i, j) = Data(i, j) - Weighted_Mean[j];
    }
  }
  
  arma::mat Weighted_Cov_Matrix(p, p);
  for (j = 0; j < p; ++j){
    for (k = 0; k < p; ++k){
      Weighted_Cov_Matrix(j, k) = 0;
      for (i = 0; i < n; ++i){
        Weighted_Cov_Matrix(j, k) = Weighted_Cov_Matrix(j, k) + (Weights[i] * Centred_Data(i, j) * Centred_Data(i, k));
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
  
  arma::vec quantile_shifted = spatialquantile(Data, Weights, u_index, c, t_vector);
  for (j = 0; j < p; ++j){
    quantile_shifted[j] = quantile_shifted[j] - Weighted_Mean[j];
  }
  
  arma::vec quantile(d_n);
  for (j = 0; j < d_n; ++j){
    quantile[j] = 0;
    for (k = 0; k < p; ++k){
      quantile[j] = quantile[j] + (quantile_shifted[k] * Eigenvectors_sorted_truncated(k, j));
    }
  }
  
  arma::mat Hessian(d_n, d_n);
  for (j = 0; j < d_n; ++j){
    for (k = 0; k < d_n; ++k){
      Hessian(j, k) = 0;
      for (i = 0; i < n; ++i){
        temp1 = 0;
        for (l = 0; l < d_n; ++l){
          temp1 = temp1 + pow(quantile[l] - Data_reduced(i, l), 2);
        }
        temp1 = sqrt(temp);
        
        if (j == k){
          Hessian(j, k) = Hessian(j, k) + (((1 / temp1) - ((1 / pow(temp1, 3)) *
            (quantile[j] - Data_reduced(i, j)) * (quantile[k] - Data_reduced(i, k)))) * Weights[i]);
        }else{
          Hessian(j, k) = Hessian(j, k) - (((1 / pow(temp1, 3)) *
            (quantile[j] - Data_reduced(i, j)) * (quantile[k] - Data_reduced(i, k))) * Weights[i]);
        }
      }
    }
  }
  
  arma::mat t1matrix(d_n, d_n);
  arma::vec t2vector(d_n);
  for (j = 0; j < d_n; ++j){
    for (k = 0; k < d_n; ++k){
      t1matrix(j, k) = 0;
      for (i = 0; i < n; ++i){
        temp1 = 0;
        for (l = 0; l < d_n; ++l){
          temp1 = temp1 + pow(quantile[l] - Data_reduced(i, l), 2);
        }
        
        t1matrix(j, k) = t1matrix(j, k) + ((1 / temp1) *
          (quantile[j] - Data_reduced(i, j)) * (quantile[k] - Data_reduced(i, k)) * Weights[i]);
      }
    }
  }
  for (j = 0; j < d_n; ++j){
    t2vector[j] = 0;
    for (i = 0; i < n; ++i){
      temp1 = 0;
      for (l = 0; l < d_n; ++l){
        temp1 = temp1 + pow(quantile[l] - Data_reduced(i, l), 2);
      }
      
      t2vector[j] = t2vector[j] + (((quantile[j] - Data_reduced(i, j)) / sqrt(temp1)) * Weights[i]);
    }
  }
  
  arma::mat CovMatrix(d_n, d_n);
  for (j = 0; j < d_n; ++j){
    for (k = 0; k < d_n; ++k){
      CovMatrix(j, k) = t1matrix(j, k) - (t2vector[j] * t2vector[k]);
    }
  }
  
  int num_Weights_positive = 0;
  for (i = 0; i < n; ++i){
    if (Weights[i] > 0){
      num_Weights_positive = num_Weights_positive + 1;
    }
  }
  double sum_Weights_square = 0;
  for (i = 0; i < n; ++i){
    sum_Weights_square = sum_Weights_square + pow(Weights[i], 2);
  }
  double sum_Weights = 1;
  E_2 = sum_Weights_square / num_Weights_positive;
  E_1 = sum_Weights / num_Weights_positive;
  
  arma::mat tempmatrix = arma::trans(arma::solve(arma::trans(Hessian), arma::trans(CovMatrix)));
  arma::mat CovQuantileMatrix = (E_2 / pow(E_1, 2)) * arma::solve(Hessian, tempmatrix);
  
  arma::vec eigenvalues_CovQuantileMatrix = arma::eig_sym(CovQuantileMatrix)
  
  
  
  return ConfidenceSet;
}


function ConfidenceSet = spatialquantileconfidenceset(Data_original, Weights, u_index, c, t_vector, alpha)
  
    E_2 = sum(Weights.^2) / sum((Weights > 0));
    E_1 = sum(Weights) / sum((Weights > 0));
    CovQuantileMatrix = (E_2 / (E_1^2)) * ( Hessian \ (CovMatrix / Hessian) );
    eigenCovQuantileMatrix = eig(CovQuantileMatrix);
    
    probvector = 1 - (1 - alpha).^( 1 ./ (2.^(1:d_n)) );
    UpperCutoffs = sqrt(eigenCovQuantileMatrix)' .* norminv((1 - (probvector / 2)), 0, 1);
    LowerCutoffs = sqrt(eigenCovQuantileMatrix)' .* norminv((probvector / 2), 0, 1);
    
    ConfidenceSet = (1 / sqrt(n)) * [UpperCutoffs; LowerCutoffs] * Eigenvectors_sorted_truncated';
    
    end
      