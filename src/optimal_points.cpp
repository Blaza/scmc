#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

NumericVector eigen_values(arma::mat M) {
  NumericVector result = wrap(arma::eig_sym(M));
  result.attr("dim") = R_NilValue;
  return result;
}

arma::mat scmc_tridiag(NumericVector alpha, NumericVector beta) {
  NumericVector sbeta = sqrt(beta);
  arma::mat out(alpha.size(), alpha.size(), arma::fill::zeros);

  for (int i = 0; i < alpha.size(); i++) {
    for (int j = 0; j <= i; j++) {
      if (i == j)
        out(i, j) = alpha[j];

      if (i - j == 1)
        out(i, j) = out(j, i) = sbeta[j];
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericVector calc_optimal_points(NumericVector moments, int N) {
  // Generate moments matrix
  arma::mat M(N + 1, N + 1, arma::fill::zeros);
  for (int i = 0; i < N + 1; i++) {
    for (int j = 0; j <= i; j++) {
      M(i,j) = M(j, i) = moments[i + j];
    }
  }

  NumericVector ev = eigen_values(M);
  if (min(ev) < 0) {
    NumericVector negatives = ev[ev < 0];

    M.diag() -= min(negatives);
  }

  // Get R matrix from cholesky decomposition
  arma::mat R = arma::chol(M);

  // Generate alpha vector
  NumericVector alpha(N);
  for (int j = 0; j < N; j++) {
    if (j == 0) {
      alpha[j] = R(j, j+1)/R(j, j);
    }
    else {
      alpha[j] = R(j, j+1)/R(j, j) - R(j-1, j)/R(j-1, j-1);
    }
  }


  // Generate beta vector
  NumericVector beta(N - 1);
  for (int j = 0; j < N - 1; j++) {
      beta[j] = R(j+1, j+1)/R(j, j) * R(j+1, j+1)/R(j, j); // squaring
  }

  // Create the tridiagonal matrix and return the eigenvalues
  return eigen_values(scmc_tridiag(alpha, beta));
}


// [[Rcpp::export]]
NumericVector sample_optimal_points(NumericVector x, int N) {
  // Generate moments vector
  NumericVector moments(2*N + 1);
  moments[0] = 1; // zeroth moment is 1
  NumericVector last_moment(x.size(), 1);
  for (int i = 1; i < 2*N + 1; i++) {
    last_moment = last_moment * x;
    moments[i] = mean(last_moment);
  }

  return calc_optimal_points(moments, N);
}