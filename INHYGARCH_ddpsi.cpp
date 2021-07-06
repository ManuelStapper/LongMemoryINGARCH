#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector INHYGARCH11ddpsi(double alpha, double beta, double d, double eta, int lagMax){
  NumericMatrix res(lagMax, 24);
  
  res(0,0) = d;   // Xi
  res(0,1) = alpha + eta*d; // Psi
  
  res(0,2) = 1; // Psi nach alpha
  res(0,3) = 0; // Psi nach beta
  res(0,4) = 1; // Xi nach d
  res(0,5) = eta; // Psi nach d
  res(0,6) = d; // Psi nach eta
  
  res(0,7) = 0; // Xi nach d^2
  res(0,8) = 0; // Psi nach alpha^2
  res(0,9) = 0; // Psi nach beta^2
  res(0,10) = 0; // Psi nach d^2
  res(0,11) = 0; // Psi nach eta^2
  res(0,12) = 0; // Psi nach alpha beta
  res(0,13) = 0; // Psi nach alpha d
  res(0,14) = 0; // Psi nach alpha eta
  res(0,15) = 0; // Psi nach beta alpha
  res(0,16) = 0; // Psi nach beta d
  res(0,17) = 0; // Psi nach beta eta
  res(0,18) = 0; // Psi nach d alpha
  res(0,19) = 0; // Psi nach d beta
  res(0,20) = 1; // Psi nach d eta
  res(0,21) = 0; // Psi nach eta alpha
  res(0,22) = 0; // Psi nach eta beta
  res(0,23) = 1; // Psi nach eta d
  
  for(int i = 1; i < lagMax; i++){
    res(i,0) = ((i-d)/(i+1))*res(i-1,0); // Xi
    res(i,1) = beta*res(i-1,1) + eta*res(i,0) - eta*(alpha + beta)*res(i - 1,0); // Psi
    
    res(i,2) = beta*res(i-1,2) - eta*res(i-1,0);
    res(i,3) = beta*res(i-1,3) + res(i-1,1) - eta*res(i-1,0);
    res(i,4) = (i-d)/(i+1)*res(i-1,4) - res(i-1,0)/(i+1);
    res(i,5) = beta*res(i-1,5) + eta*res(i,4) - eta*(alpha + beta)*res(i-1,4);
    res(i,6) = beta*res(i-1,6) + res(i,0) - (alpha + beta)*res(i-1,0);
    
    res(i,7) = (i-d)/(i+1)*res(i-1,7) - 2/(i+1)*res(i-1,4);
    res(i,8) = beta*res(i-1,8);
    res(i,9) = beta*res(i-1,9) + 2*res(i-1,3);
    res(i,10) = beta*res(i-1,10) + eta*res(i,7) - eta*(alpha + beta)*res(i-1,7);
    res(i,11) = beta*res(i-1,11);
    res(i,12) = beta*res(i-1,12) + res(i-1,2);
    res(i,13) = beta*res(i-1,13) - eta*res(i-1,4);
    res(i,14) = beta*res(i-1,14) - res(i-1,0);
    res(i,15) = beta*res(i-1,15) + res(i-1,3);
    res(i,16) = beta*res(i-1,16) + res(i-1,5) - eta*res(i-1,4);
    res(i,17) = beta*res(i-1,17) + res(i-1,6) - res(i-1,0);
    res(i,18) = beta*res(i-1,18) + eta*res(i-1,4);
    res(i,19) = beta*res(i-1,19) + res(i-1,5) - eta*res(i-1,4);
    res(i,20) = beta*res(i-1,20) + res(i,4) - (alpha + beta)*res(i-1,4);
    res(i,21) = beta*res(i-1,21) + res(i-1,0);
    res(i,22) = beta*res(i-1,22) + res(i-1,6) - res(i-1,0);
    res(i,23) = beta*res(i-1,23) + res(i,4) - (alpha + beta)*res(i-1,4);
  }
  return(res);
}