#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector C_lagrange_weights(List col_pts, IntegerMatrix indices) {
  List lambdas(clone(col_pts));
  NumericVector lambda_prods(indices.nrow());

  // Calculate lambdas as in Grzelak's paper (page 11)
  for (int z = 0; z < col_pts.size(); z++) {
    NumericVector z_vec = col_pts[z];
    NumericVector lambda_vec = lambdas[z];

    for (int l = 0; l < z_vec.size(); l++) {
      double prod = 1;
      for (int i = 0; i < z_vec.size(); i++) {
        prod *= (i != l) ? z_vec[l] - z_vec[i] : 1;
      }
      lambda_vec[l] = 1 / prod;
    }
    lambdas[z] = lambda_vec;
  }

  NumericVector lambda;
  for (int k = 0;  k < indices.nrow(); k++) {
    double prod = 1;

    for (int i = 0; i < indices.ncol(); i++) {
      lambda = lambdas[i];
      prod *= lambda[indices(k, i)];
    }

    lambda_prods[k] = prod;
  }

  return lambda_prods;
}
