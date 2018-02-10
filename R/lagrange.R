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
#' @param digits The number of decimal places to round the coefficients to.
#' @return The function which evaluates the Lagrange polynomial.
#' @export
lagrange <- function(FUN, col_pts, compile_c = FALSE, digits = 16) {
  if (is.numeric(col_pts))
    col_pts <- list(col_pts)

  col_grid <- do.call(expand.grid, col_pts)

  if (is.numeric(FUN)) {
    yvals <- FUN
  } else {
    yvals <- apply(col_grid, 1, function(row)
      do.call(FUN, unname(as.list(row))))
  }

  lm_data <- cbind(col_grid, y = yvals)

  monom_degrees <- do.call(expand.grid, lapply(col_pts, seq_along)) - 1
  formula_members <- apply(monom_degrees, 1, function(degs) {
    degs <- degs[degs > 0] # remove zero degree variables
    varnames <- names(degs)
    if(length(degs) == 0)
      return("1")

    monom_string <- paste0(ifelse(degs == 0, "", paste0(varnames, "^", degs)), collapse = " * ")
    monom_string <- gsub("\\^1", "", monom_string) # remove ^1 if present

    if (length(degs) == 1 && degs == 1)
      monom_string
    else
      sprintf("I(%s)", monom_string)
  })

  lm_formula <- sprintf("y ~ %s", paste(formula_members, collapse = " + "))

  poly_lm <- lm(lm_formula, data = lm_data)

  poly <- as.matrix(cbind(monom_degrees, coef = poly_lm$coefficients))
  #poly[abs(poly[,"coef"]) < eps, "coef"] <- 0
  poly[, "coef"] <- round(poly[, "coef"], digits)

  mpoly2function(poly, compile_c = compile_c)
}