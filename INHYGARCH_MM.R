MMinhygarch <- function(ts, lagMax4psi, lagMax4comp, lagMax){
  acfSam <- acf(ts, lag.max = lagMax4comp, plot = F)$acf[-1]
  tf <- function(theta){
    if(!checkparam(model = "inhygarch", params = c(1, theta))) return(1e16)
    psiTF <- INHYGARCH11coef(alpha = theta[1], beta = theta[2], d = theta[3], eta = theta[4], lagMax = lagMax4psi)
    acfTF <- ingarch.acf(intercept = 1, past_obs = psiTF, past_mean = NULL, lag.max = lagMax4comp, plot = F)[-1]
    return(sum((acfSam - acfTF)^2))
  }
  optRaw <- optim(par = c(0.35, 0.2, 0.3, 0.7), fn = tf, control = list(maxit = 5000))
  res <- c(NaN, optRaw$par)
  if(res[5] >= 0.9) res[5] <- 0.9
  if(res[5] <= 0.1) res[5] <- 0.1
  psires <- INHYGARCH11coef(alpha = res[2], beta = res[3], d = res[4], eta = res[5], lagMax = lagMax)
  res[1] <- mean(ts)*(1 - sum(psires))*(1 - res[3])
  return(res)
}