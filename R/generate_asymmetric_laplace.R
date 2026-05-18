#' @name generate_asymmetric_laplace
#' @title Generate Multivariate Asymmetric Laplace Random Vectors
#'
#' @description Generates a matrix of independent draws from a multivariate asymmetric Laplace distribution. Each row is \eqn{\epsilon = \alpha W + Z\sqrt{W}}, where \eqn{W \sim \text{Exp}(1)} independent of \eqn{Z \sim N_p(0, \Sigma)}. The parameter \eqn{\alpha} controls the asymmetry (skewness) in each coordinate.
#'
#' @param n Number of observations (rows) to generate.
#' @param p Dimension of each observation (columns).
#' @param alpha A numeric vector of length \code{p} giving the asymmetry (shift) parameters for each coordinate. Can be a scalar if \code{p = 1}.
#' @param Sigma A \eqn{p \times p} positive definite covariance matrix for the Gaussian component \eqn{Z}.
#'
#' @return An \code{n x p} matrix where each row is a draw from the multivariate asymmetric Laplace distribution with shift vector \code{alpha} and scale matrix \code{Sigma}.
#'
#' @importFrom MASS mvrnorm
#' @export
#'
#' @examples
#' set.seed(123)
#' p <- 3
#' alpha <- c(0.5, -0.3, 0.2)
#' Sigma <- diag(p)
#' X <- generate_asymmetric_laplace(100, p, alpha, Sigma)
#' pairs(X)
#' @importFrom stats rexp
generate_asymmetric_laplace <- function(n, p, alpha, Sigma) {
  W <- rexp(n, rate = 1)
  Z <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  epsilon <- matrix(rep(alpha, each = n), nrow = n) * W + Z * sqrt(W)
  return(epsilon)
}  