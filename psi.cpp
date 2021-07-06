#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector INFIGARCH11coef(double alpha, double beta, double d, int lagMax){
  NumericVector zeta(lagMax);
  NumericVector psi(lagMax);
  double phi = alpha + beta;
  psi(0) = phi - beta + d;
  zeta(0) = d;
  for(int i = 1; i < lagMax; i++){
    zeta(i) = ((i - d)/(i + 1))*zeta(i - 1);
    psi(i) = beta*psi(i - 1) + zeta(i) - phi*zeta(i - 1);
  }
  return(psi);
}

// [[Rcpp::export]]
NumericVector INHYGARCH11coef(double alpha, double beta, double d, double eta, int lagMax){
  NumericVector xi(lagMax);
  NumericVector psi(lagMax);
  double phi = alpha + beta;
  psi(0) = phi - beta + eta*d;
  xi(0) = d;
  for(int i = 1; i < lagMax; i++){
    xi(i) = (i - d)/(i + 1)*xi(i - 1);
    psi(i) = beta*psi(i - 1) + eta*(xi(i) - phi*xi(i - 1));
  }
  return(psi);
}