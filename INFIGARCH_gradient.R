INFIGARCH11grad <- function(tsRel, tsMat, beta0, alpha, beta, d, lagMax){
  hugeMat <- INFIGARCH11dpsi(alpha = alpha, beta = beta, d = d, lagMax = lagMax)
  
  lambdas <- c(tsMat%*%hugeMat[,2] + beta0/(1 - beta))
  dlambda <- cbind(1/(1-beta), tsMat%*%hugeMat[,c(3,4,6)])
  dlambda[,3] <- dlambda[,3] + beta0/((1 - beta)^2)
  res <- (tsRel/lambdas - 1)%*%dlambda
  return(res)
}
