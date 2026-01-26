#' @name online_sir_lfm
#' @title Online Sufficient Dimension Reduction for Laplace Factor Model (LFM)
#'
#' @description Implements an online SIR algorithm tailored for LFM data, using a proxy response constructed from the current subspace estimate and robust updates to handle heavy-tailed noise. The algorithm supports two optimization methods: gradient-based updates and perturbation-based updates.
#'
#' @param X A matrix or data stream of size n x p (rows = observations, cols = features). Can be processed row-by-row in streaming setting.
#' @param K_true Optional true dimension (for monitoring). If NULL, will estimate online via BIC-like criterion.
#' @param K_max Maximum candidate dimension for online selection (default = min(10, ncol(X))).
#' @param c_robust Robustness scale for tanh transformation (default = 1.345, approx. 0.95 efficiency for Gaussian).
#' @param eta Learning rate schedule: either a function of t, or "auto" for 1/t.
#' @param method Optimization method: "gradient" for gradient-based updates with learning rate, or "perturbation" for direct eigenvector computation of the moment matrix (default = "gradient").
#' @param verbose Logical; if TRUE, prints progress and estimated K at each step.
#'
#' @return A list with:
#' \item{B_hat}{Final estimated basis matrix (p x K_est)}
#' \item{K_est}{Estimated structural dimension}
#' \item{B_path}{List of B estimates over time (optional, for debugging)}
#' \item{loss}{Reconstruction loss trace (optional)}
#' \item{method_used}{The optimization method actually used}
#' @export
#' @examples
#' set.seed(123)
#' n <- 500; p <- 20; m <- 3
#' B_true <- qr.Q(qr(matrix(rnorm(p * m), p, m)))
#' f <- matrix(rnorm(n * m), n, m)
#' eps <- matrix(rexp(n * p, rate = 1) - 1, n, p) # Asymmetric Laplace-like noise
#' X <- f %*% t(B_true) + eps
#'
#' # Using gradient method (default)
#' out_grad <- online_sir_lfm(X, K_true = m, verbose = TRUE)
#'
#' # Using perturbation method
#' out_pert <- online_sir_lfm(X, K_true = m, method = "perturbation", verbose = TRUE)
#'
#' @importFrom stats rnorm
#'
online_sir_lfm <- function(X, K_true = NULL, K_max = NULL, c_robust = 1.345,
                           eta = "auto", method = "gradient",  
                           verbose = FALSE) {
  
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  
  if (is.null(K_max)) K_max <- min(10, p)
  if (K_max > p) K_max <- p
  
  B_hat <- qr.Q(qr(matrix(rnorm(p * K_max), p, K_max)))
  x_bar <- rep(0, p)
  M_hat <- matrix(0, p, p)
  lambda_sq_sum <- rep(0, K_max)
  
  B_path <- list()
  loss_trace <- numeric(n)
  
  if (is.character(eta) && eta == "auto") {
    get_eta <- function(t) 1 / t
  } else if (is.function(eta)) {
    get_eta <- eta
  } else {
    stop("eta must be 'auto' or a function of iteration index t")
  }
  
  valid_methods <- c("gradient", "perturbation")
  if (!method %in% valid_methods) {
    warning(paste("method", method, "not recognized. Using 'gradient'"))
    method <- "gradient"
  }
  
  for (t in 1:n) {
    x_new <- X[t, , drop = TRUE]
    
    x_bar <- ((t - 1) / t) * x_bar + (1 / t) * x_new
    
    proj <- B_hat %*% t(B_hat) %*% x_new
    y_t <- sqrt(sum(proj^2))
    
    resid <- x_new - x_bar
    r_robust <- tanh(resid / c_robust)
    
    outer_rr <- tcrossprod(r_robust)
    M_hat <- ((t - 1) / t) * M_hat + (1 / t) * outer_rr * y_t
    
    if (method == "gradient") {
      eta_t <- get_eta(t)
      grad <- M_hat %*% B_hat
      B_temp <- B_hat + eta_t * grad
      qr_temp <- qr(B_temp)
      B_hat <- qr.Q(qr_temp)
      
    } else if (method == "perturbation") {
      B_temp <- M_hat %*% B_hat
      qr_temp <- qr(B_temp)
      B_hat <- qr.Q(qr_temp)
    }
    
    Mt_proj <- t(B_hat) %*% M_hat %*% B_hat
    eig_Mt <- eigen(Mt_proj, symmetric = TRUE, only.values = TRUE)$values
    eig_Mt <- pmax(eig_Mt, 0)
    total_var <- sum(eig_Mt)
    
    if (total_var > 1e-12) {
      D_vals <- numeric(K_max)
      for (k in 1:K_max) {
        ratio <- sum(eig_Mt[1:k]^2) / total_var
        penalty <- (t^0.5) * k * (k + 1) / (2 * t)
        D_vals[k] <- ratio - penalty
      }
      K_est_t <- which.max(D_vals)
    } else {
      K_est_t <- 1
    }
    
    if (!is.null(K_true)) {
      B_trunc <- B_hat[, 1:K_est_t, drop = FALSE]
      recon <- B_trunc %*% (t(B_trunc) %*% x_new)
      loss_trace[t] <- sum((x_new - recon)^2)
    }
    
    if (verbose && (t %% 50 == 0 || t == n)) {
      cat(sprintf("Step %d: Method=%s, Est K=%d\n", t, method, K_est_t))
    }
  }
  
  final_K <- K_est_t
  B_final <- B_hat[, 1:final_K, drop = FALSE]
  
  return(list(
    B_hat = B_final,
    K_est = final_K,
    B_path = B_path,
    loss = if (!is.null(K_true)) loss_trace else NULL,
    method_used = method
  ))
}