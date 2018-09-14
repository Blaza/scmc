#' @export
newton <- function(FUN, col_pts, compile_c = FALSE) {
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

  div_diffs <- divided_diffs(xs, yvals)
  coefs <- div_diffs[1, ]

  #newton2function(xs, coefs)
  function(x) {
    newton_eval(x, xs, coefs)
  }
}
