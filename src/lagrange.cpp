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

// [[Rcpp::export]]
NumericVector lagrange_eval(NumericMatrix var_vecs, NumericMatrix nodes,
                            NumericVector y, NumericVector lambda_prods) {

  NumericVector evaluated(var_vecs.nrow());
  NumericVector vars(var_vecs.ncol());
  int n = var_vecs.nrow();
  int ysize = y.size();
  int var_count = var_vecs.ncol();

  for (int i = 0; i != n; i++) {
    vars = var_vecs(i, _);
    double bary_nominator = 0;
    double bary_denominator = 0;

    for (int j = 0; j != ysize; j++) {

      double common_element = lambda_prods[j];

      for (int k = 0; k != var_count; k++) {
        if (vars[k] != nodes(j, k)) {
          common_element /= vars[k] - nodes(j, k);
        }
        else {
          // Monkey patch! We should find an elegant way to switch to lower
          // dimensional interpolation at nodes. Right now, this works ok, as
          // the probability of evaluating at exactly a node point is near zero.
          common_element /= vars[k] - nodes(j, k) + 1e-10;
        }
      }

      bary_nominator += common_element * y[j];
      bary_denominator += common_element;
    }

    evaluated[i] = bary_nominator / bary_denominator;
  }

  return evaluated;
}
