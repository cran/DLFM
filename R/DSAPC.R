#' The distributed stochastic approximation principal component for handling online data sets with highly correlated data across multiple nodes.
#'
#' @param data is a highly correlated online data set
#' @param m is the number of principal component 
#' @param eta is the proportion of online data to total data
#' @param n1 is the length of each data subset
#' @param K is the number of nodes
#'
#' @return Asa, Dsa (lists containing results from each node)
#' @export
#'
#' @examples
#' library(LaplacesDemon)
#' library(MASS)
#' n=1000
#' p=10
#' m=5
#' mu=t(matrix(rep(runif(p,0,1000),n),p,n))
#' mu0=as.matrix(runif(m,0))
#' sigma0=diag(runif(m,1))
#' F=matrix(mvrnorm(n,mu0,sigma0),nrow=n)
#' A=matrix(runif(p*m,-1,1),nrow=p)
#' lanor <- rlaplace(n*p,0,1)
#' epsilon=matrix(lanor,nrow=n)
#' D=diag(t(epsilon)%*%epsilon)
#' data=mu+F%*%t(A)+epsilon
#' DSAPC(data=data, m=3, eta=0.8, n1=128, K=2)
DSAPC <- function(data, m, eta, n1, K) {
  n <- nrow(data)
  p <- ncol(data)
  
  Asa_list <- list()
  Dsa_list <- list()
  
  for (i in 1:K) {
    L <- matrix(rep(0, K * n1), ncol = n1)
    R <- matrix(0, n1, n)
    L[i, ] <- sample(1:n, n1, replace = FALSE)
    r <- matrix(c(1:n1, L[i, ]), ncol = n1, byrow = TRUE)
    R[t(r)] <- 1
    
    X_subset <- R %*% as.matrix(data)
    
    sapc_result <- SAPC(data = X_subset, m = m, eta = eta)
    
    Asa_list[[i]] <- sapc_result$Asa
    Dsa_list[[i]] <- sapc_result$Dsa
  }
  
  return(list(Asa = Asa_list, Dsa = Dsa_list))
}