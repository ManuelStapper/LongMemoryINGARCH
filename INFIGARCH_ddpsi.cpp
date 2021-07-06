#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector INFIGARCH11ddpsi(double alpha, double beta, double d, int lagMax){
  NumericMatrix res(lagMax, 16);
  
  res(0,0) = d;   // xi
  res(0,1) = alpha + d; // psi
  
  res(0,2) = 1; // Psi nach alpha_1
  res(0,3) = 0; // Psi nach beta_1
  res(0,4) = 1; // Xi nach d
  res(0,5) = 1; // Psi nach d
  
  res(0,6) = 0; // Xi nach d^2
  res(0,7) = 0; // psi nach alpha^2
  res(0,8) = 0; // psi nach beta^2
  res(0,9) = 0; // psi nach d^2
  res(0,10) = 0; // psi nach alpha beta
  res(0,11) = 0; // psi nach alpha d
  res(0,12) = 0; // psi nach beta alpha
  res(0,13) = 0; // psi nach beta d
  res(0,14) = 0; // psi nach d alpha
  res(0,15) = 0; // psi nach d beta
  
  for(int i = 1; i < lagMax; i++){
    res(i,0) = ((i - d)/(i + 1))*res(i - 1,0);
    res(i,1) = beta*res(i - 1,1) + res(i,0) - (alpha + beta)*res(i - 1,0);
    
    res(i,2) = beta*res(i-1,2) - res(i-1,0);
    res(i,3) = res(i-1,1) + beta*res(i-1,3) - res(i-1,0);
    res(i,4) = (i - d)/(i + 1)*res(i-1,4) - res(i-1,0)/(i+1);
    res(i,5) = beta*res(i-1,5) + res(i,4) - (alpha + beta)*res(i-1,4);
    
    res(i,6) = (i-d)/(i+1)*res(i-1,6) - 2/(i+1)*res(i-1,4);
    res(i,7) = beta*res(i-1,7);
    res(i,8) = beta*res(i-1,8) + 2*res(i-1,3);
    res(i,9) = beta*res(i-1,9) + res(i,6) - (alpha+beta)*res(i-1,6);
    res(i,10) = beta*res(i-1,10) + res(i-1,2);
    res(i,11) = beta*res(i-1,11) - res(i-1,4);
    res(i,12) = beta*res(i-1,12) + res(i-1,2);
    res(i,13) = beta*res(i-1,13) + res(i-1,5) - res(i-1,4);
    res(i,14) = beta*res(i-1,14) - res(i-1,4);
    res(i,15) = beta*res(i-1,15) + res(i-1,5) - res(i-1,4);
  }
  return(res);
}