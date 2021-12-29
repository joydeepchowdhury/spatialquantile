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
  
  arma::vec Quantile_coefficients(d_n);
  
  // Checking whether the weighted quantile is present in the data itself
  
  int Check = 0, count_weighted_norms_i = 0;
  arma::vec x(d_n);
  arma::mat U(n, d_n);
  arma::vec weights_for_i(n), norms_U(n), weighted_norms_U(n), J_i_tally(n);
  for (i = 0; i < n; ++i){
    for (j = 0; j < d_n; ++j){
      x[j] = Data(i, j);
    }
    
    for (k = 0; k < n; ++k){
      norms_U[k] = 0;
      for (j = 0; j < d_n; ++j){
        U(k, j) = Data(k, j) - x[j];
        
        norms_U[k] = norms_U[k] + pow(U(k, j), 2);
      }
      norms_U[k] = sqrt(norms_U[k]);
      
      weights_for_i[k] = Weights[k];
      
      weighted_norms_U[k] = weights_for_i[k] * norms_U[k];
      if (weighted_norms_U[k] == 0){
        count_weighted_norms_i = count_weighted_norms_i + 1;
        J_i_tally[k] = 0;
      }else{
        J_i_tally[k] = 1;
      }
    }
    
    arma::vec J_i(count_weighted_norms_i - 1), J_i_complement(n - count_weighted_norms_i);
    int J_i_index = 0, J_i_complement_index = 0;
    for (k = 0; k < n; ++k){
      if (J_i_tally[k] == 0){
        if (k != i){
          J_i[J_i_index] = k;
          J_i_index = J_i_index + 1;
        }
      }else{
        J_i_complement[J_i_complement_index] = k;
        J_i_complement_index = J_i_complement_index + 1;
      }
    }
    
    arma::mat U_new(n - count_weighted_norms_i, d_n);
    arma::vec weights_proper(n - count_weighted_norms_i);
    double sum_weights_proper = 0;
    for (k = 0; k < n - count_weighted_norms_i; ++k){
      for (j = 0; j < d_n; ++j){
        U_new(k, j) = U(J_i_complement[k], j) / norms_U[J_i_complement[k]];
      }
      weights_proper[k] = weights_for_i[J_i_complement[k]];
      
      sum_weights_proper = sum_weights_proper + weights_proper[k];
    }
    arma::vec V(d_n);
    for (j = 0; j < d_n; ++j){
      V[j] = 0;
      for (k = 0; k < n - count_weighted_norms_i; ++k){
        V[j] = V[j] + (weights_proper[k] * U_new(k, j));
      }
      V[j] = V[j] + (sum_weights_proper * u[j]);
    }
    
    double sum_weights_for_i_J_i = 0;
    for (k = 0; k < count_weighted_norms_i - 1; ++k){
      sum_weights_for_i_J_i = sum_weights_for_i_J_i + weights_for_i[J_i[k]];
    }
    
    double norm_V = 0; norm_u = 0;
    for (j = 0; j < d_n; ++j){
      norm_V = norm_V + pow(V[j], 2);
      norm_u = norm_u + pow(u[j], 2);
    }
    norm_V = sqrt(norm_V);
    norm_u = sqrt(norm_u);
    
    if (norm_V <= (1 + norm_u) * sum_weights_for_i_J_i){
      for (j = 0; j < d_n; ++j){
        Quantile_coefficients[j] = x[j];
      }
      Check = 1;
      break
    }
  }
  
  // Checking whether the data lie on a straight line, and computing the quantile in that case
  
  for (j = 0; j < d_n; ++j){
    x[j] = Data(1, j);
  }
  
  arma::vec direction_vector(d_n);
  double norm_direction_vector;
  for (i = 1; i < n; ++i){
    norm_direction_vector = 0;
    for (j = 0; j < d_n; ++j){
      direction_vector[j] = Data(i, j) - x[j];
      norm_direction_vector = norm_direction_vector + pow(direction_vector[j], 2);
    }
    norm_direction_vector = sqrt(norm_direction_vector);
    
    if (norm_direction_vector > 0){
      break;
    }
  }
  
  arma::vec s(n), s_vector(d_n);
  arma::vec z(d_n);
  bool Check_uniqueness;
  int sum_Check_if_linear = 0;
  for (i = 0; i < n; ++i){
    for (j = 0; j < d_n; ++j){
      s_vector[j] = (Data(i, j) - x[j]) / direction_vector[j];
    }
    
    Check_uniqueness = true;
    for (j = 1; j < d_n; ++j){
      if (s_vector[j] != s_vector[0]){
        Check_uniqueness = false;
        break;
      }
    }
    if (Check_uniqueness){
      sum_Check_if_linear = sum_Check_if_linear + 1;
      s[i] = s_vector[0];
    }else{
      s[i] = 0;
    }
  }
  
  if (sum_Check_if_linear == n){
    double projected_u = 0;
    for (j = 0; j < d_n; ++j){
      projected_u = projected_u + (u[j] * direction_vector[j]);
    }
    projected_u = projected_u / norm_direction_vector;
    
    double alpha = (projected_u + 1) / 2;
    
    arma::uvec s_sorted_index = arma::sort_index(s, "ascend");
    arma::vec s_sorted(n);
    for (i = 0; i < n; ++i){
      s_sorted[i] = s[s_sorted_index[i]];
    }
    
    arma::vec weights_sorted_index(n);
    for (i = 0; i < n; ++i){
      weights_sorted_index[i] = Weights[s_sorted_index[i]];
    }
    arma::vec cumulative_weights_sorted_index(n);
    int index_weighted_quantile = -1;
    for (i = 0; i < n; ++i){
      cumulative_weights_sorted_index[i] = weights_sorted_index[i];
      if (i > 0){
        cumulative_weights_sorted_index[i] = cumulative_weights_sorted_index[i] + cumulative_weights_sorted_index[i - 1];
      }
      
      if (cumulative_weights_sorted_index[i] >= alpha && index_weighted_quantile == -1){
        index_weighted_quantile = i;
      }
    }
    
    double s_weighted_quantile = s_sorted[index_weighted_quantile];
    
    for (j = 0; j < d_n; ++j){
      Quantile_coefficients[j] = x[j] + (s_weighted_quantile * direction_vector[j]);
    }
    Check = 1;
  }
  
  // Iteration procedure when the weighted quantile is not present in the data, or the data is not linear
  
  if (Check == 0){
    arma::mat X(n, d_n);
    for (i = 0; i < n; ++i){
      for (j = 0; j < d_n; ++j){
        X(i, j) = Data(i, j);
      }
    }
    
    arma::vec Q_1(d_n);
    arma::vec vector_concerned(n), vector_concerned_sorted(n);
    arma::vec weights_sorted_index(n), cumulative_weights_sorted_index(n);
    int index_weighted_quantile = -1;
    for (j = 0; j < d_n; ++j){
      for (i = 0; i < n; ++i){
        vector_concerned[i] = X(i, j);
      }
      
      arma::uvec vector_concerned_sorted_index = arma::sort_index(vector_concerned, "ascend");
      for (i = 0; i < n; ++i){
        vector_concerned_sorted[i] = vector_concerned[vector_concerned_sorted_index[i]];
        weights_sorted_index[i] = Weights[vector_concerned_sorted_index[i]];
        
        cumulative_weights_sorted_index[i] = weights_sorted_index[i];
        if (i > 0){
          cumulative_weights_sorted_index[i] = cumulative_weights_sorted_index[i] + cumulative_weights_sorted_index[i - 1];
        }
        
        if (cumulative_weights_sorted_index[i] >= (u[j] + 1) / 2 && index_weighted_quantile == -1){
          index_weighted_quantile = i;
        }
      }
      
      Q_1[j] = vector_concerned_sorted[index_weighted_quantile];
    }
    
    arma::vec Q_best_till_now(d_n);
    for (j = 0; j < d_n; ++j){
      Q_best_till_now[j] = Q_1[j];
    }
    
    double g_best_till_now = g_function_weighted(Data, Q_best_till_now, Weights, u);
    
    arma::mat Phi(d_n, d_n);
    arma::vec t1(d_n);
    double norm_sq_t1;
    for (j = 0; j < d_n; ++j){
      for (k = 0; k < d_n; ++k){
        Phi(j, k) = 0;
      }
    }
    for (i = 0; i < n; ++i){
      norm_sq_t1 = 0;
      for (j = 0; j < d_n; ++j){
        t1[j] = X(i, j) - Q_1[j];
        
        norm_sq_t1 = norm_sq_t1 + pow(t1[j], 2);
      }
      
      if (norm_sq_t1 > 0){
        for (j = 0; j < d_n; ++j){
          for (k = 0; k < d_n; ++k){
            Phi(j, k) = Phi(j, k) + (Weights[i] * (((double)(j == k) - ((t1[j] * t1[k]) / norm_sq_t1)) / sqrt(norm_sq_t1)));
          }
        }
      }
    }
    
    double sum_Weights = 0;
    for (i = 0; i < n; ++i){
      sum_Weights = sum_Weights + Weights[i];
    }
    
    double Threshold = 0.001;
    int iteration_number = 1;
    int maximum_iteration_number = 10000;
    arma::vec Delta(d_n);
    while (iteration_number <= maximum_iteration_number){
      for (i = 0; i < n; ++i){
        for (j = 0; j < d_n; ++j){
          X(i, j) = Data(i, j);
        }
      }
      
      for (j = 0; j < d_n; ++j){
        Delta[j] = 0;
        for (k = 0; k < d_n; ++k){
          Phi(j, k) = 0;
        }
      }
      
      for (i = 0; i < n; ++i){
        norm_sq_t1 = 0;
        for (j = 0; j < d_n; ++j){
          t1[j] = X(i, j) - Q_1[j];
          
          norm_sq_t1 = norm_sq_t1 + pow(t1[j], 2);
        }
        
        if (norm_sq_t1 > 0){
          for (j = 0; j < d_n; ++j){
            Delta[j] = Delta[j] + (Weights[i] * (t1[j] / sqrt(norm_sq_t1)));
            for (k = 0; k < d_n; ++k){
              Phi(j, k) = Phi(j, k) + (Weights[i] * (((double)(j == k) - ((t1[j] * t1[k]) / norm_sq_t1)) / sqrt(norm_sq_t1)));
            }
          }
        }
      }
      
      for (j = 0; j < d_n; ++j){
        Delta[j] = Delta[j] + (sum_Weights * u[j]);
      }
      
      arma::vec Q_2 = arma::solve(Phi, Delta);
      for (j = 0; j < d_n; ++j){
        Q_2[j] = Q_1[j] + Q_2[j];
      }
      
      double norm_sq_Q_1 = 0, norm_sq_Q_2 = 0, difference_relative = 0;
      for (j = 0; j < d_n; ++j){
        norm_sq_Q_1 = norm_sq_Q_1 + pow(Q_1[j], 2);
        norm_sq_Q_2 = norm_sq_Q_2 + pow(Q_2[j], 2);
        
        difference_relative = difference_relative + pow((Q_2[j] - Q_1[j]), 2);
      }
      difference_relative = sqrt(difference_relative) / max(sqrt(norm_sq_Q_1), sqrt(norm_sq_Q_2));
      
      if (difference_relative < Threshold){
        for (j = 0; j < d_n; ++j){
          Quantile_coefficients[j] = Q_2[j];
        }
        Check = 1;
        break;
      }else{
        double g_at_Q_1 = g_function_weighted(Data, Q_1, Weights, u);
        double g_at_Q_2 = g_function_weighted(Data, Q_2, Weights, u);
        
        if (g_at_Q_2 <= g_at_Q_1){
          if (g_at_Q_2 <= g_best_till_now){
            arma::vec Q_best_till_now(d_n);
            for (j = 0; j < d_n; ++j){
              Q_best_till_now[j] = Q_2[j];
            }
            double g_best_till_now = g_function_weighted(Data, Q_best_till_now, Weights, u);
          }
        }else{
          for (j = 0; j < d_n; ++j){
            Q_2[j] = ((g_at_Q_2 * Q_1[j]) + (g_at_Q_1 * Q_2[j])) / (g_at_Q_1 + g_at_Q_2);
          }
          g_at_Q_2 = g_function_weighted(Data, Q_2, Weights, u);
          if (g_at_Q_2 <= g_best_till_now){
            for (j = 0; j < d_n; ++j){
              Q_best_till_now[j] = Q_2[j];
            }
            g_best_till_now = g_function_weighted(Data, Q_best_till_now, Weights, u);
          }
        }
        
        int iteration_counter = 1;
        arma::mat Phi_temp(d_n, d_n);
        arma::vec t2(d_n);
        double norm_sq_t2;
        while (true){
          for (j = 0; j < d_n; ++j){
            for (k = 0; k < d_n; ++k){
              Phi_temp(j, k) = 0;
            }
          }
          
          for (i = 0; i < n; ++i){
            norm_sq_t2 = 0;
            for (j = 0; j < d_n; ++j){
              t2[j] = X(i, j) - Q_2[j];
              
              norm_sq_t2 = norm_sq_t2 + pow(t2[j], 2);
            }
            
            if (norm_sq_t2 > 0){
              for (j = 0; j < d_n; ++j){
                for (k = 0; k < d_n; ++k){
                  Phi_temp(j, k) = Phi_temp(j, k) +
                    (Weights[i] * (((double)(j == k) - ((t2[j] * t2[k]) / norm_sq_t2)) / sqrt(norm_sq_t2)));
                }
              }
            }
          }
          
          if ((arma::cond(Phi_temp) <= 10) || (iteration_counter > 5)){
            break;
          }else{
            for (j = 0; j < d_n; ++j){
              Q_2[j] = (Q_1[j] + Q_2[j]) / 2;
            }
            
            iteration_counter = iteration_counter + 1;
          }
        }
        
        for (j = 0; j < d_n; ++j){
          Q_1[j] = Q_2[j];
        }
      }
      
      iteration_number = iteration_number + 1;
    }
  }
  
  if (Check == 0){
    for (j = 0; j < d_n; ++j){
      Quantile_coefficients[j] = Q_best_till_now[j];
    }
      
      %% Calculating the weighted quantile
      
      Quantile = (Quantile_coefficients * Eigenvectors_sorted_truncated') + Weighted_Mean;
  }
  
  
  
  
  
  
}
