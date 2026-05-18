#' @name estimate_AD_from_OSDR
#' @title Estimate Factor Loadings and Uniquenesses from OSDR Output
#'
#' @description Extracts the estimated factor loading matrix \eqn{A} and the diagonal noise covariance matrix \eqn{D} from the output of \code{\link{online_sir_lfm}}. The subspace estimate \code{B_hat} is truncated or padded to have exactly \code{m} columns, used as \eqn{A}, and \eqn{D} is obtained by subtracting the factor covariance contribution from the sample covariance matrix, with a lower bound on the diagonal entries to ensure positivity.
#'
#' @param data The original data matrix (n x p) used in the online SIR procedure.
#' @param res A list returned by \code{\link{online_sir_lfm}}, containing at least the component \code{B_hat} (p x K_est matrix).
#' @param m The target number of factors (structural dimension). If \code{B_hat} has fewer columns, the missing directions are filled with zeros; if more, only the first \code{m} columns are retained.
#'
#' @return A list with components:
#' \item{A_est}{Estimated factor loading matrix of dimension p x m.}
#' \item{D_est}{Estimated diagonal noise covariance matrix of dimension p x p, with diagonal entries forced to be at least 0.01.}
#'
#' @export
#'
#' @importFrom stats cov
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' n <- 500; p <- 20; m <- 3
#' B_true <- qr.Q(qr(matrix(rnorm(p * m), p, m)))
#' f <- matrix(rnorm(n * m), n, m)
#' eps <- matrix(rexp(n * p, rate = 1) - 1, n, p) # Laplace noise
#' X <- f %*% t(B_true) + eps
#'
#' # Run online SIR
#' res <- online_sir_lfm(X, K_true = m, method = "gradient", verbose = FALSE)
#'
#' # Estimate A and D from the result
#' AD <- estimate_AD_from_OSDR(X, res, m)
#' print(head(AD$A_est))
#' print(diag(AD$D_est)[1:5])
#' }
#'
estimate_AD_from_OSDR <- function(data, res, m) {
  p <- ncol(data)
  B_hat <- res$B_hat  
  if (ncol(B_hat) < m) {
    B_hat <- cbind(B_hat, matrix(0, nrow = p, ncol = m - ncol(B_hat)))
  } else if (ncol(B_hat) > m) {
    B_hat <- B_hat[, 1:m, drop = FALSE]
  }
  A_est <- B_hat
  S_sample <- cov(data)
  S_factors <- A_est %*% t(A_est)
  D_est <- diag(diag(S_sample - S_factors))
  diag(D_est) <- pmax(diag(D_est), 0.01)
  return(list(A_est = A_est, D_est = D_est))
}