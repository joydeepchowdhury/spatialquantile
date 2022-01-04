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
    arma::vec h_vector(h_vector_length);
    for (i = 0; i < h_vector_length; ++i){
      h_vector[i] = min_X_distance + (i * ((max_X_distance - min_X_distance) / (h_vector_length - 1)));
    }
    
    arma::vec h_vector_check(h_vector_length);
    for (i = 0; i < h_vector_length; ++i){
      h_vector_check[i] = 1;
    }
  }
}
      
      if method_for_h == 1
      h_vector = linspace(X_distance_limits(1), X_distance_limits(2),100);
      h_vector_length = length(h_vector);
      
      h_vector_check = ones(size(h_vector));
      for i=1:h_vector_length
        h = h_vector(i);
      
      Indices_zero = (X_distance <= h);
      Indices_zero_row_sum = sum(Indices_zero,2);
      Indices_zero_row_sum_proper = Indices_zero_row_sum - 1;
      if sum(Indices_zero_row_sum_proper == 0) > 0
      h_vector_check(i) = 0;
      end
        end
        h_vector_proper = h_vector(h_vector_check > 0);
      h_vector_length_proper = length(h_vector_proper);
      
      Error_Type_temp_average = zeros(1,h_vector_length_proper);
      for i=1:h_vector_length_proper
        h = h_vector_proper(i);
      
      Type_temp = zeros(sample_size,size(Y_static,2));
      Error_Type_temp = zeros(1,sample_size);
      for j=1:1:sample_size
        Y = Y_static;
      X = X_static;
      distances = X_distance(j,:);
      distances(j) = [];
      target_Y = Y(j,:);
      target_X = X(j,:);
      Y(j,:) = [];
      X(j,:) = [];
      
      local_Y_values = Y((distances <= h),:);
      local_X_values = X((distances <= h),:);
      t_vector = 1:1:size(X,2);
      Weights = kernelweights(target_X, local_X_values, t_vector, h, Kernel);
      
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
        [~,optimum_h_index] = min(Error_Type_temp_average);
      optimum_h = h_vector_proper(optimum_h_index);
      
      optimum_h_or_neighborhood_size = optimum_h;
      else
        Number_lower = 5;
      Portion_higher = 0.5;
      Number_higher = ceil(Portion_higher * sample_size);
      Neighborhood_size_vector = (Number_lower:1:Number_higher);
      
      Error_Type_temp_average = zeros(1,length(Neighborhood_size_vector));
      for i=1:length(Neighborhood_size_vector)
        neighborhood_size = Neighborhood_size_vector(i);
      
      Type_temp = zeros(sample_size,size(Y_static,2));
      Error_Type_temp = zeros(1,sample_size);
      for j=1:1:sample_size
        Y = Y_static;
      X = X_static;
      distances = X_distance(j,:);
      distances(j) = [];
      target_Y = Y(j,:);
      target_X = X(j,:);
      Y(j,:) = [];
      X(j,:) = [];
      
      sorted_Distance_X = sort(distances, 'ascend');
      h = sorted_Distance_X(neighborhood_size);
      
      local_Y_values = Y((distances <= h),:);
      local_X_values = X((distances <= h),:);
      t_vector = 1:1:size(X,2);
      Weights = kernelweights(target_X, local_X_values, t_vector, h, Kernel);
      
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
        
        