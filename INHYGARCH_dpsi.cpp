#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector INHYGARCH11dpsi(double alpha, double beta, double d, double eta, int lagMax){
  NumericMatrix res(lagMax, 7);
  
  res(0,0) = d;   // Xi
  res(0,1) = alpha + eta*d; // Psi
  
  res(0,2) = 1; // Psi nach alpha
  res(0,3) = 0; // Psi nach beta
  res(0,4) = 1; // Xi nach d
  res(0,5) = eta; // Psi nach d
  res(0,6) = d; // Psi nach eta
  
  for(int i = 1; i < lagMax; i++){
    res(i,0) = ((i-d)/(i+1))*res(i-1,0); // Xi
    res(i,1) = beta*res(i-1,1) + eta*res(i,0) - eta*(alpha + beta)*res(i - 1,0); // Psi
    
    res(i,2) = beta*res(i-1,2) - eta*res(i-1,0);
    res(i,3) = beta*res(i-1,3) + res(i-1,1) - eta*res(i-1,0);
    res(i,4) = (i-d)/(i+1)*res(i-1,4) - res(i-1,0)/(i+1);
    res(i,5) = beta*res(i-1,5) + eta*res(i,4) - eta*(alpha + beta)*res(i-1,4);
    res(i,6) = beta*res(i-1,6) + res(i,0) - (alpha + beta)*res(i-1,0);
  }
  return(res);
}