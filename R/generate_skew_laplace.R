#' @name generate_skew_laplace
#' @title Generate Multivariate Skew Laplace Random Vectors
#'
#' @description Generates a matrix of independent draws from a multivariate skew Laplace distribution with location vector \eqn{\alpha}, scale matrix \eqn{\Sigma}, and skewness parameter \eqn{\gamma}. Each observation is constructed as \eqn{\epsilon = \alpha + \gamma V + \sqrt{V} Z}, where \eqn{V \sim \text{Gamma}((p+1)/2, 1/2)} and \eqn{Z \sim N_p(0, \Sigma)} independent of \eqn{V}.
#'
#' @param n Number of observations (rows) to generate.
#' @param p Dimension of each observation (columns).
#' @param alpha Numeric vector of length \code{p} giving the location (shift) parameter. Can be a scalar if \code{p = 1}.
#' @param gamma Numeric vector of length \code{p} giving the skewness parameter for each coordinate. Can be a scalar if \code{p = 1}.
#' @param Sigma A \eqn{p \times p} positive definite covariance matrix for the Gaussian component \eqn{Z}.
#'
#' @return An \code{n x p} matrix where each row is a draw from the specified multivariate skew Laplace distribution.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rgamma
#' @export
#'
#' @examples
#' set.seed(123)
#' p <- 3
#' alpha <- c(0.1, -0.2, 0.0)
#' gamma <- c(0.5, 0.3, -0.4)
#' Sigma <- diag(p)
#' X <- generate_skew_laplace(100, p, alpha, gamma, Sigma)
#' pairs(X)
#'
generate_skew_laplace <- function(n, p, alpha, gamma, Sigma) {
  epsilon <- matrix(0, nrow = n, ncol = p)
  shape_param <- (p + 1)/2
  for (i in 1:n) {
    V <- rgamma(1, shape = shape_param, rate = 0.5)
    Z <- MASS::mvrnorm(1, mu = rep(0, p), Sigma = Sigma)
    epsilon[i, ] <- alpha + gamma * V + sqrt(V) * Z
  }
  return(epsilon)
}