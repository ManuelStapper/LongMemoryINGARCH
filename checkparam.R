# Input:
# params: vector of parameters, for INFIGARCH (intercept, alpha, beta, d)
# for INHYGARCH (intercept, alpha, beta, d, tau)
# model:  eiter "infigarch" or "inhygarch"

# Output:
# TRUE  if all coefficients positive
# FLASE if at least one negative coefficient

checkparam <- function(params, model = c("infigarch", "inhygarch")){
  if(length(model) == 2) model <- model[1]
  if(length(params) != switch(model, infigarch = 4, inhygarch = 5)) stop("Parameter vector not suitable")
  
  if(params[1] <= 0) return(F)
  params <- params[-1]
  
  if(model == "infigarch"){
    phi <- params[1] + params[2]
    alpha <- params[1]
    beta <- params[2]
    d <- params[3]
    if(abs(beta) >= 1 | abs(phi) >= 1) return(F)
    if(beta == phi) return(F)
    psi <- INFIGARCH11coef(alpha = alpha, beta = beta, d = d, lagMax = 100)
    if(any(psi < 0)) return(F)
    f <- ((1:100) - 1 - d)/(1:100)
    
    if(beta > 0){
      if(psi[1] >= 0 & phi <= f[2]) return(T)
      else{
        kRel <- which((f[2:99] < phi)&(phi <= f[3:100])) + 2
        if(length(kRel) == 0) return(F)
        if(all(psi[kRel - 1] >= 0)) return(T)
      }
      return(F)
    }
    
    else{
      if(psi[1] >= 0 & psi[2] >= 0 & phi <= f[2]*(beta + f[3])/(beta + f[2])) return(T)
      else{
        temp1 <- (f[2:98]*(beta + f[3:99])/(beta + f[2:98]) < phi)
        kRel <- which(temp1&(phi <= f[3:99]*(beta + f[4:100])/(beta + f[3:99]))) + 3
        if(length(kRel) == 0) return(F)
        if(all(psi[kRel - 1] >= 0)) return(T)
        return(F)
      }
    }
  }
  
  if(model == "inhygarch"){
    phi <- params[1] + params[2]
    alpha <- params[1]
    beta <- params[2]
    d <- params[3]
    eta <- params[4]
    if(beta == 0 | abs(beta) >= 1 | eta <= 0 | eta >= 1) return(F)
    if(beta == phi) return(F)
    psi <- INHYGARCH11coef(alpha = alpha, beta = beta, d = d, eta = eta, lagMax = 100)
    if(any(psi < 0)) return(F)
    f <- ((1:100) - 1 - d)/(1:100)
    
    if(beta > 0){
      if(psi[1] >= 0 & phi <= f[2]) return(T)
      else{
        kRel <- which((f[2:99] < phi)&(phi <= f[3:100])) + 2
        if(length(kRel) == 0) return(F)
        if(all(psi[kRel - 1] >= 0)) return(T)
        return(F)
      }
    }
    else{
      if(psi[1]>= 0 & psi[2] >= 0 & abs(beta) <= f[2]) return(T)
      else{
        kRel <- which((f[2:99] < abs(beta)) & (abs(beta) <= f[3:100])) + 2
        if(length(kRel) == 0) return(F)
        if(!any(psi[1:(max(kRel) - 1)] < 0)) return(T)
        return(F)
      }
    }
  }
}

