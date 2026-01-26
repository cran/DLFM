#' @name PPC
#' @title Projection principal component 
#' @param data is a total  data set
#' @param m is the number of principal component
#' @return Apro, Dpro, Sigmahatpro
#' @export
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
#' PPC(data=data,m=5)
#' @importFrom stats cor cov
PPC=function(data,m){
  X=scale(data)
  n=nrow(X)
  p=ncol(X)
  P=as.matrix(diag(c(0,1),n,n))
  Xpro=scale(P%*%X)
  Sigmahatpro<-cov(Xpro)
  eig<-eigen(Sigmahatpro)
  lambdahat =eig$values[1:m]
  ind<-order(lambdahat,decreasing=T)
  lambdahat<-lambdahat[ind]
  Q <- eig$vectors
  Q<-Q[,ind]
  Qhat<-Q[,1:m]
  Apro <- matrix(0, nrow = p, ncol = m)
  for (j in 1:m) {Apro[, j] <- sqrt(lambdahat[j]) * Qhat[, j]}; Apro
  hpro <- diag(Apro %*% t(Apro))
  Dpro <- diag(Sigmahatpro - hpro)
  return(list(Apro=Apro,Dpro=Dpro,Sigmahatpro=Sigmahatpro))}

