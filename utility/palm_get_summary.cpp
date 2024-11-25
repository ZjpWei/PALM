#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

NumericMatrix invert_matrix(NumericMatrix mat) {
  // Convert R matrix to Armadillo matrix
  arma::mat arma_mat = as<arma::mat>(mat);

  // Invert the matrix using Armadillo
  arma::mat arma_inv = arma::inv(arma_mat);

  // Convert back to R matrix
  NumericMatrix inv_mat = wrap(arma_inv);

  return inv_mat;
}

NumericMatrix subset_matrix(const NumericMatrix& mat, const IntegerVector& row_indices, const IntegerVector& col_indices) {
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

// [[Rcpp::export]]

List palm_rcpp(
    List null_obj,
    List covariate_interest, // assuming covariate_interest is also passed as input
    List SUB_id,
    CharacterVector study_ID,
    CharacterVector feature_ID,
    CharacterVector cov_int_nm,
    List Sample_info
) {
  int K = feature_ID.size();

  List summary_stat_study;

  // for (int idx_d = 0; idx_d < study_ID.size(); ++idx_d) {
    int idx_d = 0;
    String d = study_ID[idx_d];
    LogicalMatrix kp_sample_id = Sample_info[d];
    NumericMatrix est_mat = na_matrix(K, cov_int_nm.size());
    NumericMatrix cov_mat = na_matrix(K, cov_int_nm.size());

    colnames(est_mat) = cov_int_nm;
    rownames(est_mat) = feature_ID;

    colnames(cov_mat) = cov_int_nm;
    rownames(cov_mat) = feature_ID;
    NumericMatrix cov_int_d =  as<NumericMatrix>(covariate_interest[d]);

    for (int cov_name_idx = 0; cov_name_idx < cov_int_nm.size(); ++cov_name_idx) {
      String cov_names = cov_int_nm(cov_name_idx);

      List null_obj_d = null_obj[d];
      std::cout << "No need to store this string";
      NumericMatrix Y_R = as<NumericMatrix>(null_obj_d["Y_R"]);
      std::cout << "No need to store this string";
      NumericMatrix Y_I = as<NumericMatrix>(null_obj_d["Y_I"]);
      NumericMatrix X = as<NumericMatrix>(null_obj_d["Z"]);
      CharacterVector SUBid = as<CharacterVector>(SUB_id[d]);

      LogicalVector kp_indices = kp_sample_id(_, cov_name_idx);
      IntegerVector kp_indices_int = seq_along(kp_indices) - 1; // indices for subsetting
      IntegerVector kp_indices_int_sub = kp_indices_int[kp_indices];

      NumericMatrix Y_R_sub = subset_matrix(Y_R, kp_indices_int_sub, Range(0,Y_R.ncol()-1));
      NumericMatrix Y_I_sub = subset_matrix(Y_I, kp_indices_int_sub, Range(0,Y_I.ncol()-1));
      NumericMatrix X_sub = subset_matrix(X,kp_indices_int_sub, Range(0,X.ncol()-1));
      CharacterVector SUBid_sub = SUBid[kp_indices_int_sub];
      CharacterVector uniq_SUBid = unique(SUBid_sub);

      NumericVector cov_col = cov_int_d[cov_name_idx];

      for (int i = 0; i < kp_indices_int_sub.size(); ++i) {
        X_sub(i, X_sub.ncol() - 1) = cov_col[i];
      }

      NumericVector ests, covs;

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

        int beta_name = X_sub.ncol() - 1;
        IntegerVector gamma_name = seq(0, X_sub.ncol() - 2);

        NumericMatrix I_gamma = subset_matrix(I_mat, gamma_name, gamma_name);
        NumericMatrix I_beta(1, 1);
        I_beta(0, 0) = I_mat(beta_name, beta_name);

        NumericMatrix Inv_I_gamma = invert_matrix(I_gamma);

        for (int l1 = 0; l1 < gamma_name.size(); ++l1) {
          for (int l2 = 0; l2 < gamma_name.size(); ++l2) {
            I_beta(0, 0) -= I_mat(beta_name, l1) * Inv_I_gamma(l1, l2) * I_mat(beta_name, l2);
          }
        }

        double est = sum(s_i_lst(_, beta_name)) / I_beta(0, 0);
        ests.push_back(est);

        double core_U = 0;
        for (int i = 0; i < s_i_lst.nrow(); ++i) {
          double tmp_U = s_i_lst(i, beta_name);
          for (int l1 = 0; l1 < gamma_name.size(); ++l1) {
            for (int l2 = 0; l2 < gamma_name.size(); ++l2) {
              tmp_U -= I_mat(beta_name, l1) * Inv_I_gamma(l1, l2) * s_i_lst(i, l2);
            }
          }
          core_U += tmp_U * tmp_U;
        }

        double cov = core_U / (I_beta(0, 0) * I_beta(0, 0));
        covs.push_back(cov);
      }

      CharacterVector feature_idd = colnames(Y_I_sub);
      for (int j = 0; j < Y_I_sub.ncol(); ++j) {
        String a = feature_idd[j];
        int tmp_id = which_character(feature_ID, a);
        est_mat(tmp_id, cov_name_idx) = ests[j];
        cov_mat(tmp_id, cov_name_idx) = covs[j];
      }
    }

    List summary_stat_study_one = List::create(
      Named("est") = est_mat,
      Named("var") = cov_mat,
      Named("n") = est_mat.nrow()
    );

    summary_stat_study[d] = cov_int_d;
  // }

  return summary_stat_study;
}
