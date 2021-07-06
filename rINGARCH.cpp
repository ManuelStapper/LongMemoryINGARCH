#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

// lambda and sigma need to be calculated in the wrapper
// [[Rcpp::export]]
NumericVector rINGARCHcppIDpois(int n, NumericVector betas, NumericVector alphas,
                                double lambda, int n_start){
  int q = betas.size() - 1;
  int p = alphas.size();
  NumericVector tsRaw(n + n_start);
  NumericVector lambdas(n + n_start);
  NumericVector out(n);
  int pq = std::max<int>(p, q);
  
  for(int r = 0; r < pq; r++){
    lambdas(r) = lambda;
    tsRaw(r) = lambda;
  }
  for(int r = pq; r < n + n_start; r++){
    lambdas(r) = betas(0);
    for(int s = 0; s < p; s++) lambdas(r) += alphas(s)*tsRaw(r - s - 1);
    for(int s = 0; s < q; s++) lambdas(r) += betas(s+1)*lambdas(r - s - 1);
    tsRaw(r) = R::rpois(lambdas(r));
    if(r >= n_start){
      out(r - n_start) = tsRaw(r);
    }
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector rINGARCHcppLOGpois(int n, NumericVector betas, NumericVector alphas,
                                 double lambda, int n_start){
  int q = betas.size() - 1 ;
  int p = alphas.size();
  NumericVector tsRaw(n + n_start);
  NumericVector lambdas(n + n_start);
  NumericVector out(n);
  int pq = std::max<int>(p, q);
  
  for(int r = 0; r < pq; r++){
    lambdas(r) = lambda;
    tsRaw(r) = lambda;
  }
  for(int r = pq; r < n + n_start; r++){
    lambdas(r) = betas(0);
    for(int s = 0; s < p; s++) lambdas(r) += alphas(s)*log(tsRaw(r - s - 1) + 1);
    for(int s = 0; s < q; s++) lambdas(r) += betas(s+1)*lambdas(r - s - 1);
    tsRaw(r) = R::rpois(exp(lambdas(r)));
    if(r >= n_start){
      out(r - n_start) = tsRaw(r);
    }
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector rINGARCHcppIDnbinom(int n, NumericVector betas, NumericVector alphas,
                                  double lambda, double size, int n_start){
  int q = betas.size() - 1;
  int p = alphas.size();
  NumericVector tsRaw(n + n_start);
  NumericVector lambdas(n + n_start);
  NumericVector out(n);
  int pq = std::max<int>(p, q);
  
  for(int r = 0; r < pq; r++){
    lambdas(r) = lambda;
    tsRaw(r) = lambda;
  }
  for(int r = pq; r < n + n_start; r++){
    lambdas(r) = betas(0);
    for(int s = 0; s < p; s++) lambdas(r) += alphas(s)*tsRaw(r - s - 1);
    for(int s = 0; s < q; s++) lambdas(r) += betas(s+1)*lambdas(r - s - 1);
    tsRaw(r) = R::rnbinom(size, size/(size + lambdas(r)));
    if(r >= n_start){
      out(r - n_start) = tsRaw(r);
    }
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector rINGARCHcppLOGnbinom(int n, NumericVector betas, NumericVector alphas,
                                 double lambda, double size, int n_start){
  int q = betas.size() - 1 ;
  int p = alphas.size();
  NumericVector tsRaw(n + n_start);
  NumericVector lambdas(n + n_start);
  NumericVector out(n);
  int pq = std::max<int>(p, q);
  
  for(int r = 0; r < pq; r++){
    lambdas(r) = lambda;
    tsRaw(r) = lambda;
  }
  for(int r = pq; r < n + n_start; r++){
    lambdas(r) = betas(0);
    for(int s = 0; s < p; s++) lambdas(r) += alphas(s)*log(tsRaw(r - s - 1) + 1);
    for(int s = 0; s < q; s++) lambdas(r) += betas(s+1)*lambdas(r - s - 1);
    tsRaw(r) = R::rnbinom(size, size/(size + exp(lambdas(r))));
    if(r >= n_start){
      out(r - n_start) = tsRaw(r);
    }
  }
  return(out);
}