#' @name frobenius.norm
#' @title Compute the Frobenius Norm of a Matrix
#'
#' @description Calculates the Frobenius norm of a matrix, defined as the square root of the sum of squares of all its entries. This is the standard Euclidean norm when the matrix is treated as a vector.
#'
#' @param M A numeric matrix.
#'
#' @return A non-negative scalar representing the Frobenius norm of \code{M}.
#'
#' @export
#'
#' @examples
#' M <- matrix(1:6, nrow = 2)
#' frobenius.norm(M)
#'
frobenius.norm <- function(M) {
  sqrt(sum(M^2))
}