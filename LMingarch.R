LMingarch <- function(ts, model = c("infigarch", "inhygarch"), lagMax = 1000,
                      parInit, controlInit = list(lagMax4psi = 100, lagMax4comp = 5),
                      gradRoot = T, testfor = NULL, testwith = c("LR", "LM"),
                      covEst = c("hessian", "sandwich"), returnAll = F){
  calcStart <- Sys.time()
  n <- length(ts)
  if(missing(parInit)) parInit <- "MM"
  lagMax4psi <- controlInit$lagMax4psi
  lagMax4comp <- controlInit$lagMax4comp
  if(length(testwith) > 1) testwith <- testwith[1]
  if(length(covEst) > 1) covEst <- covEst[1]
    
  if(length(model) > 1) model <- model[1]
  if(is.character(parInit) & length(parInit) > 1) parInit <- parInit[1]
  
  tsMat <- t(sapply((lagMax + 1):n, function(i) ts[(i-1):(i-lagMax)]))
  theta2psi <- switch(model, infigarch = function(theta){
    INFIGARCH11coef(alpha = theta[2], beta = theta[3],
                    d = theta[4], lagMax = lagMax)
  }, inhygarch = function(theta){
    INHYGARCH11coef(alpha = theta[2], beta = theta[3], d = theta[4],
                    eta = theta[5], lagMax = lagMax)
  })
  
  tf <- function(theta, model){
    if(!checkparam(params = theta, model = model)) return(1e16)
    coefTF <- theta2psi(theta = theta)
    lambdas <- tsMat%*%coefTF + theta[1]/(1 - theta[3])
    res <- sum((ts[(lagMax+1):n]*log(lambdas) - lambdas))
    if(res == -Inf | is.na(res)) return(1e16)
    return(-res)
  }
  
  if(model == "infigarch"){
    gradInner <- function(theta, model){
      INFIGARCH11grad(tsRel = ts[(lagMax+1):n], tsMat = tsMat, beta0 = theta[1], alpha = theta[2],
                      beta = theta[3], d = theta[4], lagMax = lagMax)
    }
    
    initValF <- function(ts, lagMax4psi, lagMax4comp, lagMax){
      MMinfigarch(ts = ts, lagMax4psi = lagMax4psi, lagMax4comp = lagMax4comp, lagMax = lagMax)
    }
  }
  
  if(model == "inhygarch"){
    gradInner <- function(theta, model){
      INHYGARCH11grad(tsRel = ts[(lagMax+1):n], tsMat = tsMat, beta0 = theta[1], alpha = theta[2],
                      beta = theta[3], d = theta[4], eta = theta[5], lagMax = lagMax)
    }
    
    initValF <- function(ts, lagMax4psi, lagMax4comp, lagMax){
      MMinhygarch(ts = ts, lagMax4psi = lagMax4psi, lagMax4comp = lagMax4comp, lagMax = lagMax)
    }
  }
  
  convergence <- -1
  
  if(is.numeric(parInit)){
    optRaw <- optim(par = parInit, fn = tf, gr = gradInner, model = model)
    convergence <- optRaw$convergence
    if(convergence != 0){
      warning("Bad starting values provided - Changed to 'MM'")
      parInit <- "MM"
    }
  }
  
  else{
    if(parInit == "MM"){
      parInit <- initValF(ts = ts, lagMax4psi = lagMax4psi, lagMax4comp = lagMax4comp, lagMax = lagMax)
      durationInit <- as.numeric(difftime(time1 = Sys.time(), time2 = calcStart, units = "sec"))
    }
    if(!is.numeric(parInit)) warning("Provide initial values for estimation")
    
    optRaw <- optim(par = parInit, fn = tf, gr = gradInner, model = model)
    convergence <- optRaw$convergence
    
    if(convergence == 1){
      optRaw <- optim(par = parInit, fn = tf, gr = gradInner, model = model,
                      control = list(maxit = 5000))
      convergence <- optRaw$convergence
    }
    
    if(convergence != 0){
      warning("No convergence - Try different starting values")
    }
  }
  
  if(gradRoot){
    tf2 <- function(theta, model){
      if(!checkparam(params = theta, model = model)) return(1e16)
      sum(gradInner(theta = theta, model = model)^2)
    }
    optRawGrad <- optim(par = optRaw$par, fn = tf2, gr = NULL, model = model, control = list(maxit = 5000))
  }
  
  pVal <- NaN
  tStat <- NaN
  
  if(!is.null(testfor) & testwith == "LM"){
    if(model == "infigarch"){
      if(covEst == "sandwich"){
        Sig <- INFIGARCH11sigma(tsRel = ts[-(1:lagMax)], tsMat = tsMat, beta0 = optRaw$par[1], alpha = optRaw$par[2],
                                beta = optRaw$par[3], d = optRaw$par[4], lagMax = lagMax)
      }
      if(covEst == "hessian"){
        Sig <- INFIGARCH11hess(tsRel = ts[-(1:lagMax)], tsMat = tsMat, beta0 = optRaw$par[1], alpha = optRaw$par[2],
                                  beta = optRaw$par[3], d = optRaw$par[4], lagMax = lagMax)
      }
    }
    
    if(model == "inhygarch"){
      if(covEst == "sandwich"){
        Sig <- INHYGARCH11sigma(tsRel = ts[-(1:lagMax)], tsMat = tsMat, beta0 = optRaw$par[1], alpha = optRaw$par[2],
                                beta = optRaw$par[3], d = optRaw$par[4], eta = optRaw$par[5], lagMax = lagMax)
      }
      if(covEst == "hessian"){
        Sig <- INHYGARCH11hess(tsRel = ts[-(1:lagMax)], tsMat = tsMat, beta0 = optRaw$par[1], alpha = optRaw$par[2],
                                beta = optRaw$par[3], d = optRaw$par[4], eta = optRaw$par[5], lagMax = lagMax)
      }
    }
    
    tStat <- (c(optRaw$par - testfor)%*%solve(Sig)%*%c(optRaw$par - testfor))
    pVal <- pchisq(tStat, switch(model, infigarch = 4, inhygarch = 5), lower.tail = F)
  }
  
  if(!is.null(testfor) & testwith == "LR"){
    tStat <- -2*(tf(optRaw$par, model = model) - tf(testfor, model = model))
    if(tStat < 0) tStat <- 0
    pVal <- pchisq(q = -2*(tf(optRaw$par, model = model) - tf(testfor, model = model)),
                   df = switch(model, infigarch = 4, inhygarch = 5), lower.tail = F)
  }
  
  resNames <- c("(Intercept)", "alpha", "beta", "d",
                switch(model, infigarch = NULL, inhygarch = "eta"), "Time", "tStat", "p-value")
  duration <- as.numeric(difftime(time1 = Sys.time(), time2 = calcStart, units = "sec"))
  
  if(returnAll && !gradRoot){
    return(list(parInit = c(parInit, durationInit), parLL = c(optRaw$par, duration, tStat, pVal)))
  }
  res <- c(optRaw$par, duration, tStat, pVal)
  names(res) <- resNames
  return(res)
}
