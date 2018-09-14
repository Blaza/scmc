#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector monom_eval(NumericVector x, NumericVector coefs) {
  NumericVector eval(x.size());
  int n = eval.size();
  int last_coef_index = coefs.size() - 1;

  for (int k = 0; k != n; k++) {
    for (int i = last_coef_index; i != 0; i--) {
      eval[k] = (eval[k] + coefs[i]) * x[k];
    }
    eval[k] += coefs[0]; // add the intercept
  }

  return eval;
}