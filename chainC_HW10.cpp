#include <Rcpp.h>
using namespace Rcpp;

//' @title Use the Gibbs sampler to generate a chain by C++
//' @description Use the Gibbs sampler to generate a chain by C++
//' @param N the length of chain
//' @param a Initial parameters
//' @param b Initial parameters
//' @param n Initial parameters
//' @return a Gibbs sample chain
//' @examples
//' \dontrun{
//' N <- 1e4;a<-4;b<-4;n<-10 #Initial value
//' chainsC<-chainC_HW10(N,a,b,n)
//' }
//' @export 
// [[Rcpp::export]]
NumericMatrix chainC_HW10(int N, int a,int b,int n) {
  NumericMatrix Z(N, 2);
  Z(0, 0) = n-1,Z(0, 1) = runif(1)[0];
  double x = 0, y = 0;
  for(int i = 1; i < N; i++) {
    y = Z(i-1,1);
    Z(i, 0) = rbinom(1,n,y)[0];
    x = Z(i, 0);
    Z(i,1) = rbeta(1,x+a,n-x+b)[0];
  }
  return(Z);
}