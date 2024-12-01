#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
arma::cube setSlice(arma::cube M, const NumericMatrix& M1, int i) {
  // Retrieve dimensions of M
  int dim1 = M.n_rows;
  int dim2 = M.n_cols;
  int dim3 = M.n_slices;

  // Validate slice index
  if(i < 1 || i > dim1){
    stop("Index i is out of bounds. It must be between 1 and %d.", dim1);
  }

  // Validate dimensions of M1
  if(M1.nrow() != dim2 || M1.ncol() != dim3){
    stop("Dimensions of M1 must match (dim2 x dim3). Provided M1 has dimensions (%d x %d), but expected (%d x %d).",
         M1.nrow(), M1.ncol(), dim2, dim3);
  }

  // Convert NumericMatrix M1 to arma::mat
  arma::mat M1_arma = as<arma::mat>(M1);


  // Adjust for 0-based indexing in C++
  int i_idx = i - 1;

  // Assign M1_arma to the i-th slice of M
  // Perform element-wise assignment to the i-th slice
  for(int j = 0; j < dim2; j++){
    for(int k = 0; k < dim3; k++){
      M(i_idx, j, k) = M1_arma(j, k);
    }
  }

  return M;
}

NumericMatrix invert_matrix(NumericMatrix mat) {
  // Convert R matrix to Armadillo matrix
  arma::mat arma_mat = as<arma::mat>(mat);

  // Invert the matrix using Armadillo
  arma::mat arma_inv = arma::inv(arma_mat);

  // Convert back to R matrix
  NumericMatrix inv_mat = wrap(arma_inv);

  return inv_mat;
}

NumericMatrix subset_matrix(const NumericMatrix& mat,
                            const IntegerVector& row_indices,
                            const IntegerVector& col_indices) {
  int nrows = row_indices.size();
  int ncols = col_indices.size();

  NumericMatrix result(nrows, ncols);

  // Copy values
  for (int i = 0; i < nrows; ++i) {
    for (int j = 0; j < ncols; ++j) {
      result(i, j) = mat(row_indices[i], col_indices[j]);
    }
  }

  // Set row names if they exist
  CharacterVector original_rownames = rownames(mat);
  if (original_rownames.size() > 0) {
    CharacterVector new_rownames(nrows);
    for (int i = 0; i < nrows; ++i) {
      new_rownames[i] = original_rownames[row_indices[i]];
    }
    rownames(result) = new_rownames;
  }

  // Set column names if they exist
  CharacterVector original_colnames = colnames(mat);
  if (original_colnames.size() > 0) {
    CharacterVector new_colnames(ncols);
    for (int j = 0; j < ncols; ++j) {
      new_colnames[j] = original_colnames[col_indices[j]];
    }
    colnames(result) = new_colnames;
  }

  return result;
}

NumericMatrix na_matrix(int n, int k){
  NumericMatrix m(n,k) ;
  std::fill( m.begin(), m.end(), NumericVector::get_na() ) ;
  return m ;
}

NumericMatrix create_matrix(NumericVector X) {
  int n = X.size();
  NumericMatrix mat(n, n);

  for(int l = 0; l < n; ++l) {
    for(int ll = 0; ll < n; ++ll) {
      mat(l, ll) = X[l] * X[ll];
    }
  }

  return mat;
}

int which_character(StringVector X, String a) {
  for(int i = 0; i < X.size(); ++i) {
    if(X[i] == a) {
      return i; // R is 1-indexed
    }
  }
  return NA_INTEGER; // Return NA if no match is found
}

NumericMatrix matrixMultiply(NumericMatrix M1, NumericMatrix M2) {
  int n1_rows = M1.nrow();
  int n1_cols = M1.ncol();
  int n2_rows = M2.nrow();
  int n2_cols = M2.ncol();

  // Check if the number of columns in M1 equals the number of rows in M2
  if (n1_cols != n2_rows) {
    stop("Number of columns in M1 must equal number of rows in M2 for multiplication.");
  }

  // Initialize the output matrix with dimensions n1_rows x n2_cols
  NumericMatrix M_output(n1_rows, n2_cols);

  // Perform matrix multiplication
  for(int i = 0; i < n1_rows; ++i) {
    for(int j = 0; j < n2_cols; ++j) {
      double sum = 0.0;
      for(int k = 0; k < n1_cols; ++k) {
        sum += M1(i, k) * M2(k, j);
      }
      M_output(i, j) = sum;
    }
  }

  return M_output;
}

NumericMatrix matrixSubtract(NumericMatrix M1, NumericMatrix M2) {
  int n1_rows = M1.nrow();
  int n1_cols = M1.ncol();
  int n2_rows = M2.nrow();
  int n2_cols = M2.ncol();

  // Check if dimensions match
  if (n1_rows != n2_rows || n1_cols != n2_cols) {
    stop("Matrices M1 and M2 must have the same dimensions for subtraction.");
  }

  // Initialize the output matrix with the same dimensions
  NumericMatrix M_sub(n1_rows, n1_cols);

  // Perform matrix subtraction
  for(int i = 0; i < n1_rows; ++i) {
    for(int j = 0; j < n1_cols; ++j) {
      M_sub(i, j) = M1(i, j) - M2(i, j);
    }
  }

  return M_sub;
}

NumericMatrix matrixAdd(NumericMatrix M1, NumericMatrix M2) {
  int n1_rows = M1.nrow();
  int n1_cols = M1.ncol();
  int n2_rows = M2.nrow();
  int n2_cols = M2.ncol();

  // Check if dimensions match
  if (n1_rows != n2_rows || n1_cols != n2_cols) {
    stop("Matrices M1 and M2 must have the same dimensions for adding.");
  }

  // Initialize the output matrix with the same dimensions
  NumericMatrix M_sub(n1_rows, n1_cols);

  // Perform matrix subtraction
  for(int i = 0; i < n1_rows; ++i) {
    for(int j = 0; j < n1_cols; ++j) {
      M_sub(i, j) = M1(i, j) + M2(i, j);
    }
  }

  return M_sub;
}

NumericMatrix columnSumsMatrixCpp(NumericMatrix M) {
  int n_rows = M.nrow();
  int n_cols = M.ncol();

  // Initialize the output matrix with n_cols rows and 1 column
  NumericMatrix M_sub(n_cols, 1);

  for(int j = 0; j < n_cols; ++j) {          // Iterate over columns
    double sum = 0.0;
    for(int i = 0; i < n_rows; ++i) {      // Iterate over rows
      sum += M(i, j);
    }
    M_sub(j, 0) = sum;
  }

  // Optionally, assign a column name
  M_sub.attr("dimnames") = List::create(
    R_NilValue,                // No row names
    CharacterVector::create("Sum") // Column name
  );

  return M_sub;
}

NumericMatrix transposeMatrix(NumericMatrix M) {
  int n_rows = M.nrow();
  int n_cols = M.ncol();

  // Initialize the transposed matrix with dimensions n_cols x n_rows
  NumericMatrix M_transposed(n_cols, n_rows);

  // Perform the transposition
  for(int i = 0; i < n_rows; ++i){
    for(int j = 0; j < n_cols; ++j){
      M_transposed(j, i) = M(i, j);
    }
  }

  return M_transposed;
}

// [[Rcpp::export]]

List palm_multi_rcpp(
    List null_obj,
    List SUB_id,
    CharacterVector study_ID,
    CharacterVector feature_ID,
    List Cov_int_info,
    List Sample_info
) {
  int K = feature_ID.size();

  List summary_stat_study;

  for (int idx_d = 0; idx_d < study_ID.size(); ++idx_d) {
    String d = study_ID[idx_d];
    CharacterVector cov_int_nm = Cov_int_info[d];
    LogicalVector kp_indices = Sample_info[d];
    NumericMatrix est_mat = na_matrix(K, cov_int_nm.size());
    arma::cube cov_mat(K, cov_int_nm.size(), cov_int_nm.size(), arma::fill::zeros);

    colnames(est_mat) = cov_int_nm;
    rownames(est_mat) = feature_ID;

    List null_obj_d = null_obj[d];
    NumericMatrix Y_R = as<NumericMatrix>(null_obj_d["Y_R"]);
    NumericMatrix Y_I = as<NumericMatrix>(null_obj_d["Y_I"]);
    NumericMatrix X = as<NumericMatrix>(null_obj_d["Z"]);
    CharacterVector SUBid = as<CharacterVector>(SUB_id[d]);

    IntegerVector kp_indices_int = seq_along(kp_indices) - 1; // indices for subsetting
    IntegerVector kp_indices_int_sub = kp_indices_int[kp_indices];


    NumericMatrix Y_R_sub = subset_matrix(Y_R, kp_indices_int_sub, Range(0,Y_R.ncol()-1));
    NumericMatrix Y_I_sub = subset_matrix(Y_I, kp_indices_int_sub, Range(0,Y_I.ncol()-1));
    NumericMatrix X_sub = subset_matrix(X,kp_indices_int_sub, Range(0,X.ncol()-1));
    CharacterVector SUBid_sub = SUBid[kp_indices_int_sub];

    CharacterVector uniq_SUBid = unique(SUBid_sub);

    for (int k = 0; k < Y_I_sub.ncol(); ++k) {
      NumericMatrix s_i_lst(uniq_SUBid.size(), X_sub.ncol());
      NumericMatrix I_mat(X.ncol(), X_sub.ncol());

      for (int l = 0; l < uniq_SUBid.size(); ++l) {
        NumericVector s_i_SUB(X.ncol());
        for (int i = 0; i < SUBid_sub.size(); ++i) {
          if (uniq_SUBid[l] == SUBid_sub[i]) {
            s_i_SUB += Y_R_sub(i, k) * X_sub(i, _);
            NumericMatrix XX = create_matrix(X_sub(i, _));
            I_mat += Y_I_sub(i, k) * XX;
          }
        }
        s_i_lst(l, _) = s_i_SUB;
      }

      NumericMatrix s_i_lst_t = transposeMatrix(s_i_lst);
      IntegerVector beta_name = seq(cov_int_nm.size(), X_sub.ncol()) - 1;
      IntegerVector gamma_name = seq(0, X_sub.ncol() - beta_name.size() - 1);

      NumericMatrix I_gamma = subset_matrix(I_mat, gamma_name, gamma_name);
      NumericMatrix I_left = subset_matrix(I_mat, beta_name, gamma_name);
      NumericMatrix I_right = subset_matrix(I_mat, gamma_name, beta_name);
      NumericMatrix Inv_I_gamma = invert_matrix(I_gamma);

      NumericMatrix I_left_inv = matrixMultiply(I_left, Inv_I_gamma);
      NumericMatrix I_sub = matrixMultiply(I_left_inv, I_right);
      NumericMatrix I_beta = subset_matrix(I_mat, beta_name, beta_name);
      I_beta = matrixSubtract(I_beta, I_sub);
      NumericMatrix Inv_I_beta = invert_matrix(I_beta);

      IntegerVector indice_s =  seq(0, s_i_lst.nrow() - 1);
      NumericMatrix est_prod = columnSumsMatrixCpp(subset_matrix(s_i_lst, indice_s, beta_name));
      NumericMatrix est_d = matrixMultiply(Inv_I_beta, est_prod);
      est_mat(k, _) = est_d;

      NumericMatrix core_U(beta_name.size(), beta_name.size());
      for(int i = 0; i < s_i_lst.nrow(); ++i){
        IntegerVector idx_i = seq(i, i);
        NumericMatrix tmp_U = matrixSubtract(subset_matrix(s_i_lst_t, beta_name, idx_i),
                                             matrixMultiply(I_left_inv,
                                                            subset_matrix(s_i_lst_t, gamma_name, idx_i)));
        core_U = matrixAdd(core_U, matrixMultiply(tmp_U, transposeMatrix(tmp_U)));
      }

      NumericMatrix cov_d = matrixMultiply(matrixMultiply(Inv_I_beta, core_U), Inv_I_beta);

      cov_mat = setSlice(cov_mat, cov_d, k + 1);
    }

    NumericVector M = wrap(cov_mat);
    M.attr("dim") = Dimension(K, cov_int_nm.size(), cov_int_nm.size());
    List dim_names = List::create(feature_ID, cov_int_nm, cov_int_nm);
    M.attr("dimnames") = dim_names;

    List summary_stat_study_one = List::create(
      Named("est") = est_mat,
      Named("var") = M,
      Named("n") = est_mat.nrow()
    );

    summary_stat_study[d] = summary_stat_study_one;
  }

  return summary_stat_study;
}


