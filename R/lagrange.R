#' Generate a Lagrange polynomial
#'
#' Generates a function which evaluates the Lagrange polynomial of the given function and given node
#' points
#'
#' @param FUN The function to interpolate.
#' @param col_pts The node points at which to calculate FUN, where the poly must equal FUN. Should
#'   be a list of node points for multivariate polynomial.
#' @return The function which evaluates the Lagrange polynomial.
#' @export
lagrange <- function(FUN, col_pts) {
  if (is.numeric(col_pts))
    col_pts <- list(col_pts)

  node_matrix <- as.matrix(do.call(expand.grid, col_pts))

  if (is.numeric(FUN)) {
    yvals <- FUN
  } else {
    yvals <- apply(node_matrix, 1, function(row)
      do.call(FUN, unname(as.list(row))))
  }

  lm_data <- cbind(col_grid, y = yvals)

  function (...) {
    args <- list(...)
    var_vecs <- do.call(cbind, args)

    lagrange_eval(var_vecs, node_matrix, yvals, lambda_prods)
  }
}