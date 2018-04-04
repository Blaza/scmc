#' Calculate optimal collocation points from (sample) moments of a distribution
#'
#' Calculates optimal collocation points from (sample) moments of the approximating distribution
#'
#' Note: Using sample moments for a large sample from the approximating random variable gives
#' similar results as using theoretical. For example, using sample = rnorm(1e6) gives similar points
#' as using moment_fun = scmc::normal_moment_fun.
#'
#' @param N The number of collocation points to calculate
#' @param moment_fun A function which takes one (unsigned) integer k and return the k-th moment of
#'   the distribution to use for collocation. Only 2*N moments are needed.
#' @param sample If moment_fun not specified, the sample from which to calculate the collocation
#'   points using the 2*N sample moments.
#' @return The N optimal collocation points for interpolation used in SCMC.
#' @export
optimal_points <- function(N, moment_fun = NULL, sample = NULL) {
  if (!is.null(sample))
    return(sample_optimal_points(sample, N))

  if (is.null(moment_fun))
    stop("Either a moment function or a sample must be provided.")

  moments <- c(1, sapply(1:(2*N), moment_fun))
  calc_optimal_points(moments, N)
}


#' Theoretical moments calculation
#'
#' @param k The moment which to calculate
#' @param sd The standard deviation of the normal distibution (mean is assumed 0)
#' @param rate The rate parameter of the exponential distribution
#' @param supp The vector (a, b) indicating the support of the uniform distribution
#' @return The k-th moment of the distribution.
#'
#' @export
norm_moment_fun <- function(k, sd = 1) {
  if (k %% 2 == 1)
    return(0)

  n <- k/2
  sd^k * 2^n * gamma(n + 0.5) /sqrt(pi)
}

#' @rdname norm_moment_fun
#' @export
exp_moment_fun <- function(k, rate = 1) {
  factorial(k) / (rate^k)
}

#' @rdname norm_moment_fun
#' @export
unif_moment_fun <- function(k, supp = c(0, 1)) {
  a <- supp[1]
  b <- supp[2]

  n <- k+1
  (b^n - a^n)/(n*(b-a))
}


#' Chebyshev Nodes
#'
#' Calculate Chebyshev points for interpolation
#'
#' @param n The number of points to generate
#' @param interval The interval over which to calculate Chebyshev nodes
#' @return A numeric vector containing n Chebyshev nodes in the interval 'interval'
#' @export
chebyshev_nodes <- function(n, interval = c(-1, 1)) {
  k <- seq_len(n)
  a <- interval[1]
  b <- interval[2]

  nodes <- (a+b)/2 + cos(pi * (2*k - 1) / (2*n)) * (b-a)/2
  nodes
}
