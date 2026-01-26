#' @name osdr_lfm
#' @title Online Sufficient Dimension Reduction for Laplace Factor Models (OSDR-LFM)
#'
#' @description Implements an online SIR-based sufficient dimension reduction method tailored for Laplace Factor Models (LFM) with symmetric, asymmetric, or skewed error structures.
#' Supports distributed deployment via local updates and global aggregation.
#'
#' @param X numeric matrix (n x p), observations in rows.
#' @param Y optional numeric vector (n) of proxy responses (e.g., factor scores). 
#'        If NULL, uses norm of projection as proxy (unsupervised LFM mode).
#' @param laplace_type character; one of "symmetric", "asymmetric", or "skewed".
#' @param K_max integer; maximum candidate dimension (default = min(10, p)).
#' @param H integer; number of slices for SIR (default = max(5, floor(sqrt(n)))).
#' @param method_svd character; "perturbation" or "gradient" (default = "gradient").
#' @param is_distributed logical; if TRUE, simulate distributed node behavior.
#' @param node_id integer; node identifier (only used if is_distributed = TRUE).
#' @param sync_interval integer; how often to "aggregate" in distributed mode (ignored if not distributed).
#' @param verbose logical; print progress.
#'
#' @return list with B_hat (p x K_est), K_est, lambda_trace, and (if distributed) local_B.
#'
#' @examples
#' set.seed(42)
#' n <- 600; p <- 30; m <- 4
#' A <- qr.Q(qr(matrix(rnorm(p * m), p, m)))
#' F <- matrix(rnorm(n * m), n, m)
#' eps <- matrix(rexp(n * p) - rexp(n * p), n, p)
#' X <- F %*% t(A) + eps
#' 
#' out <- osdr_lfm(X, laplace_type = "asymmetric", K_max = 6, verbose = TRUE)
#' cat("Estimated K:", out$K_est, "\n")
#' 
#' @importFrom stats quantile rnorm
#' @export
osdr_lfm <- function(X, Y = NULL, laplace_type = c("symmetric", "asymmetric", "skewed"),
                     K_max = NULL, H = NULL, method_svd = c("gradient", "perturbation"),
                     is_distributed = FALSE, node_id = 1, sync_interval = 50,
                     verbose = FALSE) {
  
  laplace_type <- match.arg(laplace_type)
  method_svd <- match.arg(method_svd)
  
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  
  if (is.null(K_max)) K_max <- min(10, p)
  if (K_max > p) K_max <- p
  if (is.null(H)) H <- max(5L, floor(sqrt(n)))
  
  # --- Step 0: Initialize ---
  B_t <- qr.Q(qr(matrix(rnorm(p * K_max), p, K_max)))  # p x K_max
  x_bar <- rep(0, p)
  Sigma_t <- diag(p) * 1e-6  # online covariance (regularized)
  M_t <- matrix(0, p, p)     # kernel matrix
  Gamma_t <- matrix(0, p, p) # running average of M (for perturbation method)
  
  # Parameters for asymmetric/skewed Laplace
  alpha_t <- rep(0, p)  # location shift (asymmetric)
  gamma_t <- rep(0, p)  # skewness (skewed)
  
  # Storage
  lambda_trace <- matrix(0, n, K_max)
  local_B_list <- list()
  
  # Slice boundaries (fixed in advance for stability)
  if (is.null(Y)) {
    # Unsupervised: use rolling projection norm as proxy
    Y_proxy <- numeric(n)
  } else {
    Y_proxy <- Y
  }
  
  for (t in 1:n) {
    x_new <- X[t, , drop = TRUE]
    
    # ---- Update mean and covariance (with Laplace adaptation) ----
    x_bar_old <- x_bar
    x_bar <- ((t - 1) / t) * x_bar + (1 / t) * x_new
    
    # Robust residual based on Laplace type
    resid <- x_new - x_bar
    if (laplace_type == "asymmetric") {
      # Estimate location shift alpha_t online
      alpha_t <- ((t - 1) / t) * alpha_t + (1 / t) * resid
      resid_adj <- resid - alpha_t
    } else if (laplace_type == "skewed") {
      # Simple skewness proxy: use sign-weighted magnitude
      skew_weight <- ifelse(resid > 0, 1.5, 0.5)
      gamma_t <- ((t - 1) / t) * gamma_t + (1 / t) * skew_weight * resid
      resid_adj <- resid - 2 * gamma_t
    } else {
      resid_adj <- resid
    }
    
    # Online covariance (rank-1 update)
    diff <- x_new - x_bar_old
    Sigma_t <- ((t - 1) / t) * Sigma_t + (1 / t) * tcrossprod(diff)
    
    # ---- Proxy response construction (if unsupervised) ----
    if (is.null(Y)) {
      proj_norm <- sqrt(sum((B_t %*% crossprod(B_t, x_new))^2))
      Y_proxy[t] <- proj_norm
    }
    y_t <- Y_proxy[t]
    
    # ---- Slice assignment (quantile-based for robustness) ----
    if (t < H) {
      slice_id <- 1
    } else {
      quantiles <- quantile(Y_proxy[1:t], probs = seq(0, 1, length.out = H + 1), na.rm = TRUE)
      slice_id <- findInterval(y_t, quantiles, all.inside = TRUE)
      slice_id <- max(1, min(H, slice_id))
    }
    
    # ---- Compute slice mean incrementally (simplified: single active slice) ----
    m_h <- solve(Sigma_t + diag(p) * 1e-6, resid_adj)  # approx \Sigma^{-1}(x - \bar{x})
    
    # Update kernel matrix: M_t = sum_h m_h m_h^T (approx with current h)
    M_t <- ((t - 1) / t) * M_t + (1 / t) * tcrossprod(m_h)
    
    # ---- Online SVD update ----
    if (method_svd == "gradient") {
      eta_t <- 1 / t
      grad <- M_t %*% B_t
      B_temp <- B_t + eta_t * grad
      B_t <- qr.Q(qr(B_temp))
    } else if (method_svd == "perturbation") {
      Gamma_t <- ((t - 1) / t) * Gamma_t + (1 / t) * M_t
      # Perturbation requires eigen-decomp at each step (costly for large p)
      eig_M <- eigen(M_t, symmetric = TRUE)
      lambda_t <- eig_M$values[1:K_max]
      beta_t <- eig_M$vectors[, 1:K_max, drop = FALSE]
      B_t <- beta_t
      Gamma_t <- M_t  # simplify: treat Gamma_t = M_t
    }
    
    # ---- Online eigenvalue trace (for BIC) ----
    Mt_proj <- crossprod(B_t, M_t %*% B_t)
    eig_vals <- eigen(Mt_proj, symmetric = TRUE, only.values = TRUE)$values
    eig_vals <- pmax(eig_vals, 0)
    lambda_trace[t, ] <- eig_vals
    
    # ---- Distributed simulation: store local B periodically ----
    if (is_distributed && (t %% sync_interval == 0 || t == n)) {
      local_B_list[[length(local_B_list) + 1]] <- B_t
    }
    
    if (verbose && (t %% 100 == 0 || t == n)) {
      cat(sprintf("Step %d | Slice: %d | Laplace: %s\n", t, slice_id, laplace_type))
    }
  }
  
  # ---- Final dimension selection (BIC-like) ----
  total_var <- rowSums(lambda_trace^2)
  D_vals <- matrix(0, n, K_max)
  for (k in 1:K_max) {
    ratio <- rowSums(lambda_trace[, 1:k, drop = FALSE]^2) / pmax(total_var, 1e-12)
    penalty <- (sqrt(1:n) * k * (k + 1)) / (2 * (1:n))
    D_vals[, k] <- ratio - penalty
  }
  K_est <- which.max(D_vals[n, ])
  B_final <- B_t[, 1:K_est, drop = FALSE]
  
  # Ensure orthonormality
  B_final <- qr.Q(qr(B_final))
  
  result <- list(
    B_hat = B_final,
    K_est = K_est,
    lambda_trace = lambda_trace,
    laplace_type = laplace_type
  )
  
  if (is_distributed) {
    result$local_B_history <- local_B_list
    result$node_id <- node_id
  }
  
  return(result)
}