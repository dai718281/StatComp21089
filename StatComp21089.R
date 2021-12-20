#' @title Benchmark R functions
#' @name benchmarks
#' @description The gradient descent functions and GnS function to deal with the impact of truncation data
#' @importFrom stats integrate
#' @useDynLib StatComp21089
#' @examples
#' \dontrun{
#' x<-matrix(iris[1:50,1],50,1)
#' y<-matrix(iris[1:50,2],50,1)
#' GradientDescent(x,y,1e-14,1e4,stepmethod=T,step=0.001,alpha=0.25,beta=0.8)
#' GradientDescent(x,y,1e-14,1e4,stepmethod=F,step=0.001,alpha=0.25,beta=0.8)
#' US<- matrix(1,nrow=100,ncol=3)
#' U<-V<-Y<-rep(-1,100)
#' US[,1]<- Y;US[,2]<- U;US[,3]<- V
#' TS<- GnS(US)
#' }
NULL

#' @title Gradient Descent
#' @param x a matrix
#' @param y Response variable
#' @param error Termination condition, the range of two adjacent search results
#' @param maxiter The maximum number of iterations
#' @param stepmethod If True,function uses backward descent method to find the step length,else uses fixed step
#' @param step Fixed step
#' @param alpha Parameter of the retrospective descent method
#' @param beta Parameter of the retrospective descent method
#' @return coefficient,iter and error of gradient descent functions
#' @examples
#' \dontrun{
#' x<-matrix(iris[1:50,1],50,1)
#' y<-matrix(iris[1:50,2],50,1)
#' GradientDescent(x,y,1e-14,1e4,stepmethod=T,step=0.001,alpha=0.25,beta=0.8)
#' GradientDescent(x,y,1e-14,1e4,stepmethod=F,step=0.001,alpha=0.25,beta=0.8)
#' }
#' @export
GradientDescent<-function(x,y,error,maxiter,stepmethod=T,step=0.001,alpha=0.25,beta=0.8)
{
  m<-nrow(x)
  x<-cbind(matrix(1,m,1),x)
  n<-ncol(x)
  theta<-matrix(rep(0,n),n,1)  #theta initial value is set to 0
  iter<-1
  newerror<-1
  
  while((newerror>error)|(iter<maxiter)){
    iter<-iter+1
    h<-x%*%theta  
    des<-t(t(h-y)%*%x)  #gradient
    #Backward descent method to find the step length
    if(stepmethod==T){
      sstep=1
      new_theta<-theta-sstep*des
      new_h<-x%*%new_theta
      costfunction<-t(h-y)%*%(h-y)  #Least squares loss function
      new_costfunction<-t(new_h-y)%*%(new_h-y)
      #Backtracking descent method to find the step size sstep
      while(new_costfunction>costfunction-alpha*sstep*sum(des*des)){
        sstep<-sstep*beta
        new_theta<-theta-sstep*des
        new_h<-x%*%new_theta
        new_costfunction<-t(new_h-y)%*%(new_h-y)  
      }
      newerror<-t(theta-new_theta)%*%(theta-new_theta)       
      theta<-new_theta     
    }
    
    #Set fixed step directly
    if(stepmethod==F){        
      new_theta<-theta-step*des
      new_h<-x%*%new_theta
      # new_costfunction<-t(new_h-y)%*%(new_h-y)
      newerror<-t(theta-new_theta)%*%(theta-new_theta)
      theta<-new_theta 
    }
    
  }
  costfunction<-t(x%*%theta-y)%*%(x%*%theta-y)
  result<-list(theta,iter,costfunction)
  names(result)<-c('coefficient','iter','error')
  result
}

#' @title Deal with the impact of truncation data on Y
#' @param x a matrix(x=(Y,U,V)) Y is respond variable and (U,V) is double truncated about Y.
#' @return a vector
#' @examples
#' \dontrun{
#' US<- matrix(1,nrow=100,ncol=3)
#' U<-V<-Y<-rep(-1,100)
#' US[,1]<- Y;US[,2]<- U;US[,3]<- V
#' TS<- GnS(US)
#' }
#' @export
GnS=function(x)    # x=(Y,U,V)
{
  n<- nrow(x)
  H<- h1<- h2<- 1:n
  F<- f<- 1:n
  I=function(t)
  {
    if(t[2]<=t[1]&t[1]<=t[3]) {return(1)}
    else {return(0)}
  }
  for(i in 1:n) {h1[i]<- 1/n}
  for(i in 1:n)
  {
    H[i]<- 0
    for(j in 1:n) {H[i]<- H[i]+h1[j]*I(c(x[j,1],x[i,2],x[i,3]))}
  }
  while(abs(max(h1-h2))>=10^(-6))
  {
    h2<- h1
    for(j in 1:n)
    {
      f[j]<- 0
      for(i in 1:n) {f[j]<- f[j]+1/H[i]}
      f[j]<- (1/f[j])*(1/H[j])
    }
    for(i in 1:n)
    {
      F[i]<- 0
      for(j in 1:n){F[i]<- F[i]+f[j]*I(c(x[i,1],x[j,2],x[j,3]))}
    }
    for(j in 1:n)
    {
      h1[j]<- 0
      for(i in 1:n){h1[j]<- h1[j]+1/F[i]}
      h1[j]<- (1/h1[j])*(1/F[j])
    }
    for(i in 1:n)
    {
      H[i]<- 0
      for(j in 1:n){H[i]<- H[i]+h1[j]*I(c(x[j,1],x[i,2],x[i,3]))}
    }
  }
  return(f)
}