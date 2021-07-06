#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector INFIGARCH11dpsi(double alpha, double beta, double d, int lagMax){
  NumericMatrix res(lagMax, 6);
  
  res(0,0) = d;   // Xi
  res(0,1) = alpha + d; // Psi
  
  res(0,2) = 1; // Psi nach alpha_1
  res(0,3) = 0; // Psi nach beta_1
  res(0,4) = 1; // Xi nach d
  res(0,5) = 1; // Psi nach d
  
  
  for(int i = 1; i < lagMax; i++){
    res(i,0) = ((i - d)/(i + 1))*res(i - 1,0);
    res(i,1) = beta*res(i - 1,1) + res(i,0) - (alpha + beta)*res(i - 1,0);
    res(i,2) = beta*res(i-1,2) - res(i-1,0);
    res(i,3) = res(i-1,1) + beta*res(i-1,3) - res(i-1,0);
    res(i,4) = (i-d)/(i+1)*res(i-1,4) - res(i-1,0)/(i+1);
    res(i,5) = beta*res(i-1,5) + res(i,4) - (alpha + beta)*res(i-1,4);
  }
  return(res);
}