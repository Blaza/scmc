#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix divided_diffs(NumericVector xs, NumericVector ys) {
  NumericMatrix diffs(xs.size(), xs.size());

  diffs(_, 0) = ys;

  for (int j = 1; j < xs.size(); j++) {
    for (int i = 0; i < xs.size() - j; i++) {
      diffs(i, j) = (diffs(i + 1, j - 1) - diffs(i, j - 1)) /
                    (xs[i + j] - xs[i]);
    }
  }

  return diffs;
}
