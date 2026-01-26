#' The stochastic approximation principal component can handle online data sets with highly correlated.
#'
#' @param data is a highly correlated online data set
#' @param m is the number of principal component 
#' @param eta is the proportion of online data to total data
#'
#' @return Asa,Dsa
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
#' SAPC(data=data,m=3,eta=0.8) 
SAPC<-function(data,m,eta){
  X<-as.matrix(scale(data))
  S<-cov(X)   
  n<-nrow(X)          
  n0<-round(eta*n)  
  p<-ncol(X)           
  Xbar<-colMeans(X[1:n0,])      
  eig1<-eigen(cov(X[1:n0,]-Xbar)) 
  lambda<-eig1$values[1:m]       
  V<-eig1$vectors[,1:m]      
  V1<-V
  
  for (i in (n0+1):n) {
    Xbar<-((i-1)/i)*Xbar+(1/i)*X[i,]
    Xcenter<-t(t(X[i,]-Xbar))
    gamma<-1/(i^(0.73))
    V<-V+gamma*Xcenter%*%t(Xcenter)%*%V
    V<-qr.Q(qr(V))
    lambda<- lambda+gamma*(t(V)%*%Xcenter%*%t(Xcenter)%*%V)
  }
  V2<-V
  Asa<-matrix(0,nrow=p,ncol=m)
  for (j in 1:m){
    Asa[,j]<-sqrt(lambda[j])*V2[,j]
  } 
  h2<-diag(Asa%*%t(Asa))
  Dsa<-diag(S-h2) 
  return(list(Asa=Asa,Dsa=Dsa))
}
