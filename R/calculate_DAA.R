#' @name calculate_DAA
#' @title Distance Between Column Spaces of Two Loading Matrices
#'
#' @description Computes a normalized distance between the subspaces spanned by the columns of two factor loading matrices \code{Ahat} (estimated) and \code{A1} (true/reference). The metric is based on the projection matrix of \code{A1}, regularized to avoid singularity, and is scaled by the number of factors \code{m} and squared number of variables \code{p^2} to lie in \eqn{[0,1]}.
#'
#' @param Ahat Estimated factor loading matrix of dimension \code{p x m}, or at least having \code{m} columns.
#' @param A1 True or reference factor loading matrix of dimension \code{p x m} (or \code{p x m1}, but typically \code{m} columns).
#' @param m Number of factors (structural dimension) used in the scaling factor.
#' @param p Number of variables (observed dimensions).
#'
#' @return A scalar in \eqn{[0,1]} representing the discrepancy between the column spaces of \code{Ahat} and \code{A1}. A value of 0 indicates identical subspaces; values closer to 1 indicate greater discrepancy.
#'
#' @details The distance is defined as
#' \deqn{D(A_{\mathrm{est}}, A_{\mathrm{true}}) = \sqrt{1 - \frac{1}{m p^2} \, \mathrm{Re}\left( \mathrm{tr}\left( A_{\mathrm{est}}^\top (A_{\mathrm{true}} A_{\mathrm{true}}^\top + \epsilon I)^{-1} A_{\mathrm{true}} A_{\mathrm{true}}^\top A_{\mathrm{est}} \right) \right)},}
#' where \eqn{\epsilon = 10^{-6}} is a small ridge parameter to ensure invertibility of \eqn{A_{\mathrm{true}} A_{\mathrm{true}}^\top}.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' p <- 10; m <- 3
#' # True loading matrix (orthonormal columns)
#' A_true <- qr.Q(qr(matrix(rnorm(p * m), p, m)))
#' # Slightly perturbed estimate
#' A_est <- A_true + matrix(rnorm(p * m, sd = 0.05), p, m)
#' A_est <- qr.Q(qr(A_est))  # re-orthogonalize
#' 
#' calculate_DAA(A_est, A_true, m, p)
#'
calculate_DAA <- function(Ahat, A1, m, p) {
  epsilon <- 1e-6  
  A1tA1 <- A1 %*% t(A1) + epsilon * diag(nrow(A1))
  A1tA1_inv <- solve(A1tA1)
  result <- (1 - 1 / (m * (p)^2) * Re(t(Ahat) %*% A1tA1_inv %*% A1 %*% t(A1) %*% Ahat))^(1/2)
  return(result)
}