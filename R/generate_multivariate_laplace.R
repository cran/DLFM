#' @name generate_multivariate_laplace
#' @title Generate Multivariate Symmetric Laplace Random Vectors
#'
#' @description Generates a matrix of independent draws from a multivariate symmetric Laplace distribution. Each row is an independent observation of a \eqn{p}-dimensional random vector \eqn{\epsilon = Z\sqrt{W}}, where \eqn{Z \sim N_p(0, \Sigma_L)} and \eqn{W \sim \text{Exp}(1)} independent of \eqn{Z}.
#'
#' @param n Number of observations (rows) to generate.
#' @param p Dimension of each observation (columns).
#' @param Sigma_L A \eqn{p \times p} positive definite covariance matrix for the Gaussian component \eqn{Z}.
#'
#' @return An \code{n x p} matrix where each row is a draw from the multivariate symmetric Laplace distribution with scale matrix \eqn{\Sigma_L}.
#'
#' @importFrom MASS mvrnorm
#' @export
#'
#' @examples
#' set.seed(123)
#' p <- 3
#' Sigma_L <- diag(p)
#' X <- generate_multivariate_laplace(100, p, Sigma_L)
#' pairs(X)
#' @importFrom stats rexp
generate_multivariate_laplace <- function(n, p, Sigma_L) {
  W <- rexp(n, rate = 1)
  Z <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma_L)
  epsilon <- Z * sqrt(W)
  return(epsilon)
}