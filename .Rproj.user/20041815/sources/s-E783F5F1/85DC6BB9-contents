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
  
  int i, j, k;
  
  // Checking whether there is only one distinct observation in Data
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
  
  return ConfidenceSet;
}


function ConfidenceSet = spatialquantileconfidenceset(Data_original, Weights, u_index, c, t_vector, alpha)
  
  z = Data_original(1,:);
Difference = ones(size(Data_original,1),1) * z - Data_original;
norm_Difference = sqrt(trapz(t_vector, Difference.^2, 2));
if sum(norm_Difference) == 0
ConfidenceSet = ones(2,1) * z;
return
  end

  n = size(Data_original,1);
if sum(Weights) ~= 1
Weights = Weights / sum(Weights);
end
  
  t_1 = sqrt(n);
t_2 = 2 * n^(1/3);
t_3 = min(t_1,t_2);
d_n = floor(t_3);

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

Quantile = (spatialquantile(Data_original, Weights, u_index, c, t_vector) - Weighted_Mean)...
  * Eigenvectors_sorted_truncated;
  
  Hessian = zeros(d_n);
  for i=1:n
    Hessian = Hessian + ( (1 / sqrt( sum((Quantile - Data(i,:)).^2) )) * eye(d_n)...
      - (1 / ( sqrt( sum((Quantile - Data(i,:)).^2) ) )^3) *...
      ( (Quantile - Data(i,:))' * (Quantile - Data(i,:)) ) ) * Weights(i);
  end
    Hessian = Hessian / sum(Weights);
  
  t1matrix = zeros(d_n);
  t2vector = zeros(1, d_n);
  for i=1:n
    t1matrix = t1matrix + ( (1 / sum((Quantile - Data(i,:)).^2)) *...
      ( (Quantile - Data(i,:))' * (Quantile - Data(i,:)) ) ) * Weights(i);
  
  t2vector = t2vector + ( (Quantile - Data(i,:)) / sqrt( sum((Quantile - Data(i,:)).^2) ) )...
    * Weights(i);
    end
      t1matrix = t1matrix / sum(Weights);
    t2vector = t2vector / sum(Weights);
    CovMatrix = t1matrix - (t2vector' * t2vector);
    
    E_2 = sum(Weights.^2) / sum((Weights > 0));
    E_1 = sum(Weights) / sum((Weights > 0));
    CovQuantileMatrix = (E_2 / (E_1^2)) * ( Hessian \ (CovMatrix / Hessian) );
    eigenCovQuantileMatrix = eig(CovQuantileMatrix);
    
    probvector = 1 - (1 - alpha).^( 1 ./ (2.^(1:d_n)) );
    UpperCutoffs = sqrt(eigenCovQuantileMatrix)' .* norminv((1 - (probvector / 2)), 0, 1);
    LowerCutoffs = sqrt(eigenCovQuantileMatrix)' .* norminv((probvector / 2), 0, 1);
    
    ConfidenceSet = (1 / sqrt(n)) * [UpperCutoffs; LowerCutoffs] * Eigenvectors_sorted_truncated';
    
    end
      