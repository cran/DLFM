#' @name LFM
#' @title Generate Laplace factor models
#' @description The function is to generate Laplace factor model data.
#' The function supports various distribution types for generating the data,
#' including:
#' - `truncated_laplace`: Truncated Laplace distribution
#' - `log_laplace`: Univariate Symmetric Log-Laplace distribution
#' - `Asymmetric Log_Laplace`: Log-Laplace distribution
#' - `Skew-Laplace`: Skew-Laplace distribution
#' @usage LFM(n, p, m, distribution_type)
#' @param n An integer specifying the sample size.
#' @param p An integer specifying the sample dimensionality or the number of variables.
#' @param m An integer specifying the number of factors in the model.
#' @param distribution_type A character string indicating the type of distribution to use
#' for generating the data.
#' @return A list containing the following elements:
#' \item{data}{A numeric matrix of the generated data.}
#' \item{A}{A numeric matrix representing the factor loadings.}
#' \item{D}{A numeric matrix representing the uniquenesses, which is a diagonal matrix.}
#' @examples
#' n <- 1000
#' p <- 10
#' m <- 5
#' sigma1 <- 1
#' sigma2 <- matrix(c(1,0.7,0.7,1), 2, 2)
#' distribution_type <- "Asymmetric Log_Laplace"
#' results <- LFM(n, p, m, distribution_type)
#' print(results)
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom LaplacesDemon rllaplace
#' @importFrom LaplacesDemon rallaplace
#' @importFrom LaplacesDemon rslaplace
#' @importFrom LaplacesDemon rlaplace
#' @importFrom LaplacesDemon lower.triangle
#' @importFrom LaplacesDemon is.positive.definite
#' @importFrom LaplacesDemon is.symmetric.matrix
#' @importFrom LaplacesDemon upper.triangle
#' @importFrom LaplacesDemon is.square.matrix
#' @importFrom matrixcalc frobenius.norm
#' @importFrom stats runif
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats pcauchy
#' @importFrom stats qcauchy
LFM <- function(n, p, m, distribution_type) {
  mu <- t(matrix(rep(runif(p, 0, 1000), n), p, n))
  mu0 <- as.matrix(runif(m, 0))
  sigma0 <- diag(runif(m, 1))
  sigma1 <- 1
  sigma2 <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
  F <- matrix(mvrnorm(n, mu0, sigma0), nrow = n)
  A <- matrix(runif(p * m, -1, 1), nrow = p)
  lanor <- numeric(0)
  if (distribution_type == "truncated_laplace") {
    mu1 <- c(0, 1)
    lower <- c(-2, -3)
    upper <- c(3, 3)
    truncated_normal <- function(n_samples, mean, sd, lower_bound, upper_bound) {
      p_lower <- pnorm(lower_bound, mean, sd)
      p_upper <- pnorm(upper_bound, mean, sd)
      u <- runif(n_samples, p_lower, p_upper)
      return(qnorm(u, mean, sd))
    }
    norm_dim1 <- truncated_normal(n * p / 2, mu1[1], 1, lower[1], upper[1])
    norm_dim2 <- truncated_normal(n * p / 2, mu1[2], 1, lower[2], upper[2])
    laplace_from_normal <- function(norm_values) {
      cauchy_values <- qcauchy(pnorm(norm_values, 0, 1), location = 0, scale = 1)
      return(cauchy_values * 0.7)
    }
    
    laplace_dim1 <- laplace_from_normal(norm_dim1)
    laplace_dim2 <- laplace_from_normal(norm_dim2)
    
    lanor <- c(laplace_dim1, laplace_dim2)
    
  } else if (distribution_type == "log_laplace") {
    lanor <- rllaplace(n * p, 0, sigma1)
  } else if (distribution_type == "Asymmetric Log_Laplace") {
    lanor <- rallaplace(n * p, 0, sigma1, 1)
  } else if (distribution_type == "Skew-Laplace") {
    lanor <- rslaplace(n * p, 0, 1, 1)
  } else {
    stop("Unknown distribution type")
  }
  
  epsilon <- matrix(lanor, nrow = n)
  D <- diag(t(epsilon) %*% epsilon)
  data <- mu + F %*% t(A) + epsilon
  
  return(list(data = data, A = A, D = D))
}