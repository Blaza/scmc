#' @export
monomial <- function(FUN, col_pts) {
  if (is.numeric(col_pts))
    col_pts <- list(col_pts)

  node_matrix <- as.matrix(do.call(expand.grid, col_pts))

  if (is.numeric(FUN)) {
    yvals <- FUN
  } else {
    yvals <- apply(node_matrix, 1, function(row)
      do.call(FUN, unname(as.list(row))))
  }

  xs <- col_pts[[1]]
  vm_mat <- sapply(seq_along(xs) - 1, function(k) xs^k)

  coefs <- solve(vm_mat, yvals)


  function(x) {
    monom_eval(x, coefs)
  }
}