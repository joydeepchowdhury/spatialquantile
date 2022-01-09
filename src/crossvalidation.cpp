# include <RcppArmadillo.h>
# include <algorithm>
# include <string.h>

// [[Rcpp::depends(RcppArmadillo)]]

// The function `crossvalidation` calculates the optimum bandwidth or
// neighborhood size using cross validation. If `method_for_h` = 1, it returns
// the optimum bandwidth, else it returns the optimum neighborhood size.
// `X_static` is the `n`-by-`p` data matrix of `n` covariate observations,
// `Y_static` is the `n`-by-`q` data matrix of the corresponding response
// observations. Three types of estimators can be used for cross validations,
// viz., the weighted coordinate-wise mean, the weighted coordinate-wise median
// and the weighted spatial median. To use the weighted coordinate-wise mean,
// put `type` = 'coord_mean', to use the weighted coordinate-wise median, put
// `type` = 'coord_median', and to use the weighted spatial median, put `type`
// = 'spatial_median'. The numbers `n`, `p` and `q` must be greater than 1.
// `t_vector_X` and `t_vector_Y` are the respective grids on which the
// observations in `X_static` and `Y_static` are recorded. The argument `Kernel`
// specifies the name of the kernel function, The following kernels are
// implemented:
//    * `gaussian`: (1 / sqrt(2 pi)) exp(- u^2 / 2),
//    * `triangular`: (1 - abs(u)) (abs(u) <= 1),
//    * `epanechnikov`: (3/4) (1 - u^2) (abs(u) <= 1),
//    * `quartic`: (15/16) ((1 - u^2)^2) (abs(u) <= 1),
//    * `triweight`: (35/32) ((1 - u^2)^3) (abs(u) <= 1),
//    * `tricube`: (70/81) (1 - abs(u)^3)^3 (abs(u) <= 1),
//    * `uniform` (default): 0.5 (abs(u) <= 1).


// [[Rcpp::export]]
double Lp_norm(arma::vec t_vector, arma::vec X, double p){
  int q = X.n_elem;
  
  int i, j;
  double lp_norm, temp;
  
  temp = 0;
  if (p < arma::datum::inf){
    for (j = 0; j < (q - 1); ++j){
      temp = temp + ((pow(abs(X[j + 1]), p) + pow(abs(X[j]), p)) / 2) * (t_vector[j + 1] - t_vector[j]);
    }
    lp_norm = pow(temp, (1/p));
  }else{
    for (j = 0; j < (q - 1); ++j){
      temp = max(temp, abs(X[j]));
    }
    lp_norm = temp;
  }
  
  return lp_norm;
}

// [[Rcpp::export]]
double crossvalidation(arma::vec t_vector_X, arma::mat X_static,
                       arma::vec t_vector_Y, arma::mat Y_static,
                       int method_for_h, char *type, char *Kernel){
  double optimum_h_or_neighborhood_size;
  
  int sample_size = X_static.n_rows;
  int q_cov = X_static.n_cols, q_res = Y_static.n_cols;
  int p_covariate = 2, p_response = 2;
  
  int i, j, k, l;
  double inf = arma::datum::inf;
  
  arma::mat X_distance(sample_size, sample_size);
  double min_X_distance = 0, max_X_distance = 0;
  arma::vec t1(q_cov);
  for (i = 0; i < sample_size; ++i){
    for (j = 0; j < sample_size; ++j){
      if (i < j){
        for (k = 0; k < q_cov; ++k){
          t1[k] = X_static(i, k) - X_static(j, k);
        }
        X_distance(i, j) = Lp_norm(t_vector_X, t1, p_covariate);
      }else if (i == j){
        X_distance(i, j) = 0;
      }else{
        X_distance(i, j) = X_distance(j, i);
      }
      min_X_distance = min(min_X_distance, X_distance(i, j));
      max_X_distance = max(max_X_distance, X_distance(i, j));
    }
  }
  
  if (method_for_h == 1){
    int h_vector_length = 100;
    int h_vector_length_proper = h_vector_length;
    double h;
    arma::vec h_vector(h_vector_length);
    arma::vec h_vector_check(h_vector_length);
    for (i = 0; i < h_vector_length; ++i){
      h_vector[i] = min_X_distance + (i * ((max_X_distance - min_X_distance) / (h_vector_length - 1)));
      
      h_vector_check[i] = 1;
    }
    
    arma::vec Indices_zero_row_sum_proper(sample_size);
    for (i = 0; i < h_vector_length; ++i){
      h = h_vector[i];
      
      for (j = 0; j < sample_size; ++j){
        Indices_zero_row_sum_proper[j] = 0;
        for (k = 0; k < sample_size; ++k){
          if (X_distance(j, k) <= h){
            Indices_zero_row_sum_proper[j] = Indices_zero_row_sum_proper[j] + 1;
          }
        }
        Indices_zero_row_sum_proper[j] = Indices_zero_row_sum_proper[j] - 1;
        
        if (Indices_zero_row_sum_proper[j] == 0){
          h_vector_check[i] = 0;
          
          h_vector_length_proper = h_vector_length_proper - 1;
          
          break;
        }
      }
    }
    
    arma::mat h_vector_proper(h_vector_length_proper);
    i = 0;
    j = 0;
    while (j < h_vector_length){
      if (h_vector_check[j] > 0){
        h_vector_proper[i] = h_vector[j];
        i = i + 1;
        j = j + 1;
      }else{
        j = j + 1;
      }
    }
    
    int X_distance_h_count, local_values_current_index;
    arma::vec target_Y(q_res), target_X(q_cov);
    arma::vec X_distance_check(sample_size);
    arma::vec Error_Type_temp_average(h_vector_length_proper);
    arma::mat Type_temp(sample_size, q_res);
    arma::vec Error_Type_temp(sample_size);
    for (i = 0; i < h_vector_length_proper; ++i){
      h = h_vector_proper[i];
      
      for (j = 0; j < sample_size; ++j){
        for (l = 0; l < q_res; ++l){
          target_Y[l] = Y_static(j, l);
        }
        for (l = 0; l < q_cov; ++l){
          target_X[l] = X_static(j, l);
        }
        
        X_distance_h_count = 0;
        for (k = 0; k < sample_size; ++k){
          if (X_distance(j, k) <= h && k != j){
            X_distance_check[k] = 1;
            
            X_distance_h_count = X_distance_h_count + 1;
          }else{
            X_distance_check[k] = 0;
          }
        }
        
        arma::mat local_Y_values(X_distance_h_count, q_res);
        arma::mat local_X_values(X_distance_h_count, q_cov);
        local_values_current_index = 0;
        for (k = 0; k < sample_size; ++k){
          if (X_distance_check[k] == 1){
            for (l = 0; l < q_res; ++l){
              local_Y_values(local_values_current_index, l) = Y_static(k, l);
            }
            for (l = 0; l < q_cov; ++l){
              local_X_values(local_values_current_index, l) = X_static(k, l);
            }
            
            local_values_current_index = local_values_current_index + 1;
          }
        }
        
        Weights = kernelweights(target_X, local_X_values, t_vector_X, h, Kernel);
        
        if (strncmp(type, "coord_median", 20) == 0){
          arma::vec weighted_median(q_res);
          arma::vec vector_concerned(X_distance_h_count), vector_concerned_sorted(X_distance_h_count);
          arma::vec weights_by_sorted_index(X_distance_h_count), cumulative_weights_by_sorted_index(X_distance_h_count);
          for (k = 0; k < q_res; ++k){
            for (l = 0; l < X_distance_h_count; ++l){
              vector_concerned[l] = local_Y_values(l, k);
            }
            arma::uvec vector_concerned_sorted_index = arma::sort_index(vector_concerned, "ascend");
            for (l = 0; l < X_distance_h_count; ++l){
              vector_concerned_sorted[l] = vector_concerned[vector_concerned_sorted_index[l]];
              weights_by_sorted_index[l] = Weights[vector_concerned_sorted_index[l]];
              
              cumulative_weights_by_sorted_index[l] = weights_by_sorted_index[l];
              if (l > 0){
                cumulative_weights_by_sorted_index[l] = cumulative_weights_by_sorted_index[l] +
                  cumulative_weights_by_sorted_index[l - 1];
              }
              if (cumulative_weights_by_sorted_index[l] >= 0.5){
                index_weighted_quantile = l;
                weighted_median[k] = vector_concerned_sorted[index_weighted_quantile];
                break;
              }
            }
          }
          
          for (k = 0; k < q_res; ++k){
            Type_temp(j, k) = weighted_median[k];
          }
        }else if (strncmp(type, "coord_mean", 20) == 0){
          arma::vec weighted_mean(q_res);
          double sum_Weights;
          for (k = 0; k < q_res; ++k){
            weighted_mean[k] = 0;
            sum_Weights = 0;
            for (l = 0; l < X_distance_h_count; ++l){
              weighted_mean[k] = weighted_mean[k] + (Weights[l] * local_Y_values(l, k));
              sum_Weights = sum_Weights + Weights[l];
            }
            weighted_mean[k] = weighted_mean[k] / sum_Weights;
          }
          
          for (k = 0; k < q_res; ++k){
            Type_temp(j, k) = weighted_mean[k];
          }
        }else{    // strncmp(type, "spatial_median", 20) == 0
          arma::vec spatialmedian = spatialquantile(local_Y_values, Weights, 0, 0, t_vector_Y);
          
          for (k = 0; k < q_res; ++k){
            Type_temp(j, k) = spatialmedian[k];
          }
        }
        
        arma::vec y(q_res);
        for (k = 0; k < q_res; ++k){
          y[k] = target_Y[k] - Type_temp(j, k);
        }
        
        Error_Type_temp[j] = Lp_norm(t_vector_Y, y, p_response);
      }
      
      Error_Type_temp_average[i] = 0;
      for (j = 0; j < sample_size; ++j){
        Error_Type_temp_average[i] = Error_Type_temp_average[i] + Error_Type_temp[j];
      }
      Error_Type_temp_average[i] = Error_Type_temp_average[i] / sample_size;
    }
    
    int optimum_h_index = 0;
    for (i = 0; i < h_vector_length_proper; ++i){
      if (Error_Type_temp_average[i] < Error_Type_temp_average[optimum_h_index]){
        optimum_h_index = i;
      }
    }
    
    double optimum_h = h_vector_proper[optimum_h_index];
    
    optimum_h_or_neighborhood_size = optimum_h;
  }else{
    int Number_lower = 5;
    double Portion_higher = 0.5;
    double h;
    int Number_higher = (int)ceil(Portion_higher * sample_size);
    int length_Neighborhood_size_vector = Number_higher - Number_lower + 1;
    arma::vec Neighborhood_size_vector(length_Neighborhood_size_vector);
    for (i = 0; i < length_Neighborhood_size_vector; ++i){
      Neighborhood_size_vector[i] = Number_lower + i;
    }
    
    arma::vec Error_Type_temp_average(length_Neighborhood_size_vector);
    int neighborhood_size;
    arma::mat Type_temp(sample_size, q_res);
    arma::vec Error_Type_temp(sample_size);
    arma::vec target_Y(q_res), target_X(q_cov);
    arma::vec distances(sample_size - 1);
    int X_distance_h_count, local_values_current_index;
    arma::vec X_distance_check(sample_size);
    for (i = 0; i < length_Neighborhood_size_vector; ++i){
      neighborhood_size = Neighborhood_size_vector[i];
      
      for (j = 0; j < sample_size; ++j){
        for (l = 0; l < q_res; ++l){
          target_Y[l] = Y_static(j, l);
        }
        for (l = 0; l < q_cov; ++l){
          target_X[l] = X_static(j, l);
        }
        
        for (k = 0; k < (j - 1); ++k){
          distances[k] = X_distance(j, k);
        }
        for (k = j; k < sample_size; ++k){
          distances[k - 1] = X_distance(j, k);
        }
        
        arma::vec sorted_Distance_X = arma::sort(distances, "ascend");
        
        h = sorted_Distance_X[neighborhood_size - 1];
        
        X_distance_h_count = 0;
        for (k = 0; k < sample_size; ++k){
          if (X_distance(j, k) <= h && k != j){
            X_distance_check[k] = 1;
            
            X_distance_h_count = X_distance_h_count + 1;
          }else{
            X_distance_check[k] = 0;
          }
        }
        
        arma::mat local_Y_values(X_distance_h_count, q_res);
        arma::mat local_X_values(X_distance_h_count, q_cov);
        local_values_current_index = 0;
        for (k = 0; k < sample_size; ++k){
          if (X_distance_check[k] == 1){
            for (l = 0; l < q_res; ++l){
              local_Y_values(local_values_current_index, l) = Y_static(k, l);
            }
            for (l = 0; l < q_cov; ++l){
              local_X_values(local_values_current_index, l) = X_static(k, l);
            }
            
            local_values_current_index = local_values_current_index + 1;
          }
        }
        
        Weights = kernelweights(target_X, local_X_values, t_vector_X, h, Kernel);
        
        if (strncmp(type, "coord_median", 20) == 0){
          arma::vec weighted_median(q_res);
          arma::vec vector_concerned(X_distance_h_count), vector_concerned_sorted(X_distance_h_count);
          arma::vec weights_by_sorted_index(X_distance_h_count), cumulative_weights_by_sorted_index(X_distance_h_count);
          for (k = 0; k < q_res; ++k){
            for (l = 0; l < X_distance_h_count; ++l){
              vector_concerned[l] = local_Y_values(l, k);
            }
            arma::uvec vector_concerned_sorted_index = arma::sort_index(vector_concerned, "ascend");
            for (l = 0; l < X_distance_h_count; ++l){
              vector_concerned_sorted[l] = vector_concerned[vector_concerned_sorted_index[l]];
              weights_by_sorted_index[l] = Weights[vector_concerned_sorted_index[l]];
              
              cumulative_weights_by_sorted_index[l] = weights_by_sorted_index[l];
              if (l > 0){
                cumulative_weights_by_sorted_index[l] = cumulative_weights_by_sorted_index[l] +
                  cumulative_weights_by_sorted_index[l - 1];
              }
              if (cumulative_weights_by_sorted_index[l] >= 0.5){
                index_weighted_quantile = l;
                weighted_median[k] = vector_concerned_sorted[index_weighted_quantile];
                break;
              }
            }
          }
          
          for (k = 0; k < q_res; ++k){
            Type_temp(j, k) = weighted_median[k];
          }
        }else if (strncmp(type, "coord_mean", 20) == 0){
          arma::vec weighted_mean(q_res);
          double sum_Weights;
          for (k = 0; k < q_res; ++k){
            weighted_mean[k] = 0;
            sum_Weights = 0;
            for (l = 0; l < X_distance_h_count; ++l){
              weighted_mean[k] = weighted_mean[k] + (Weights[l] * local_Y_values(l, k));
              sum_Weights = sum_Weights + Weights[l];
            }
            weighted_mean[k] = weighted_mean[k] / sum_Weights;
          }
          
          for (k = 0; k < q_res; ++k){
            Type_temp(j, k) = weighted_mean[k];
          }
        }else{    // strncmp(type, "spatial_median", 20) == 0
          arma::vec spatialmedian = spatialquantile(local_Y_values, Weights, 0, 0, t_vector_Y);
          
          for (k = 0; k < q_res; ++k){
            Type_temp(j, k) = spatialmedian[k];
          }
        }
        
        arma::vec y(q_res);
        for (k = 0; k < q_res; ++k){
          y[k] = target_Y[k] - Type_temp(j, k);
        }
        
        Error_Type_temp[j] = Lp_norm(t_vector_Y, y, p_response);
      }
        
      }
      
    }
    
  }
}
      
      
      for j=1:1:sample_size
        
      
      if strcmp(type, 'coord_median') == 1
      weighted_median = zeros(1,size(local_Y_values,2));
      for k=1:size(local_Y_values,2)
        vector_concerned = local_Y_values(:,k);
      [vector_concerned_sorted, vector_concerned_sorted_index] = ...
        sort(vector_concerned,'ascend');
      weights_by_sorted_index = Weights(vector_concerned_sorted_index);
      cumulative_weights_by_sorted_index = cumsum(weights_by_sorted_index);
      index_weighted_quantile = find(cumulative_weights_by_sorted_index >= 0.5, 1);
      weighted_median(k) = vector_concerned_sorted(index_weighted_quantile);
      end
        Type_temp(j,:) = weighted_median;
      elseif strcmp(type, 'coord_mean') == 1
      local_Y_values_weighted = (Weights * ones(1,size(local_Y_values,2)))...
                                                                          .* local_Y_values;
      Type_temp(j,:) = mean(local_Y_values_weighted,1);
      elseif strcmp(type, 'spatial_median') == 1
      Type_temp(j,:) = spatialquantile(local_Y_values, Weights, 0, 0, t_vector_Y);
      else
        error('error: Enter correct type.')
        end
        
        Error_Type_temp(j) = Lp_norm(t_vector_Y, (target_Y - Type_temp(j,:)), p_response);
      end
        Error_Type_temp_average(i) = mean(Error_Type_temp);
      end
        [~,optimum_neighborhood_size_index] = min(Error_Type_temp_average);
      optimum_neighborhood_size = Neighborhood_size_vector(optimum_neighborhood_size_index);
      
      optimum_h_or_neighborhood_size = optimum_neighborhood_size;
      end
        
        end
        
        