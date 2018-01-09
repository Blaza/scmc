#' Generate a Lagrange polynomial
#'
#' Generates a function which evaluates the Lagrange polynomial of the given function and given node
#' points
#'
#' @param FUN The function to interpolate.
#' @param col_pts The node points at which to calculate FUN, where the poly must equal FUN. Should
#'   be a list of node points for multivariate polynomial.
#' @param compile_c Logical indicating whether to compile the function using Rcpp, resulting in
#'   faster evaluation later on.
#' @return The function which evaluates the Lagrange polynomial.
#' @export
lagrange <- function(FUN, col_pts, compile_c = FALSE) {
  if (is.numeric(col_pts))
    col_pts <- list(col_pts)

  if (is.numeric(FUN)) {
    yvals <- FUN
  } else {
    yvals <- apply(do.call(expand.grid, col_pts), 1, function(row)
      do.call(FUN, unname(as.list(row))))
  }

  indices <- as.matrix(do.call(expand.grid, lapply(col_pts, seq_along)) - 1)
  lambda_prods <- C_lagrange_weights(col_pts, indices)

  # initialise main polynomial
  poly <- list(c(coef = 0))
  # pre-generate mpoly objects which represent variables, so they can be operated on as polynomials
  # this is to improve speed, as calling mp(...) takes a long time so it should be in the loop
  var_polys <- list()
  for (d in seq_along(col_pts)) {
    var_polys[[d]] <- mp(paste0("x", d))
  }

  for (i in seq_along(yvals)) {
    poly_elem <- yvals[i] * lambda_prods[i]
    for (d in seq_along(col_pts)) {
      for (l in seq_along(col_pts[[d]])) {
        if (l - 1 != indices[i, d])
          poly_elem <- poly_elem * (var_polys[[d]] + (-col_pts[[d]][l]))
      }
    }
    poly <- poly + poly_elem
  }

  mpoly2function(poly, compile_c = compile_c)
}