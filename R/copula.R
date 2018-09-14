#' Generate a sample from a copula
#'
#' Generates a sample from the specified copula object
#'
#' @param copula An R object of class "Copula"
#' @param N An integer indicating the number of collocation points to use for each dimension
#' @return A function which takes one parameter ('n') and returns a matrix with dimensions (n,
#'   dim(copula)) containing the sample from the copula, each row being one observation.
#' @export
copula_sampler <- function(copula, N, gss = 1.37) {
  if (!requireNamespace("copula", quietly = TRUE)) {
    stop("Package copula is needed to work with copulas. Please install it.",
         call. = FALSE)
  }

  col_pts <- optimal_points(N, norm_moment_fun)

  x_vals <- do.call(expand.grid, rep(list(col_pts), dim(copula)))
  u_vals <- pnorm(as.matrix(x_vals), sd = gss)
  y_vals <- log(-log(copula::cCopula(u_vals, copula, inverse = TRUE)))

  polys <- lapply(2:dim(copula), function(i) {
    lagrange(y_vals[1 : (N^i), i], rep(list(col_pts), i))
  })

  # sampler function
  function(n) {
    copula_sample <- matrix(numeric(n*dim(copula)), nrow = n)
    copula_sample[,1] <- runif(n)
    cond_samples <- list(qnorm(copula_sample[,1], sd = gss))

    for (i in 2:dim(copula)) {
      margin_sample <- exp(-exp(do.call(polys[[i-1]], c(cond_samples, list(rnorm(n, sd = gss))))))
      copula_sample[,i] <- margin_sample
      cond_samples <- c(cond_samples, list(qnorm(margin_sample, sd = gss)))
    }

    copula_sample
  }
}