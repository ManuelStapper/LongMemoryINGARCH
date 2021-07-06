INFIGARCH11hess <- function(tsRel, tsMat, beta0, alpha, beta, d, lagMax){
  n <- length(tsRel + lagMax)
  hugeMat <- INFIGARCH11ddpsi(alpha = alpha, beta = beta, d = d, lagMax = lagMax)
  
  lambdas <- c(tsMat%*%hugeMat[,2] + beta0/(1 - beta))
  dlambda <- cbind(1/(1-beta), tsMat%*%hugeMat[,c(3,4,6)])
  dlambda[,3] <- dlambda[,3] + beta0/((1 - beta)^2)
  
  dl <- (tsRel/lambdas - 1)%*%dlambda
  
  # ddomega <- rbind(cbind(matrix(c(0, 1/((1-beta)^2), 1/((1-beta)^2), 2*beta0/((1-beta)^3)), nrow = 2, ncol = 2),0,0),0,0)
  
  # ddlambda <- cbind(0, beta0/((1-beta)^2), tsMat%*%hugeMat[,c(8,11,12,13,9,14,15,16,10)])
  ddlambda <- cbind(0, 1/((1-beta)^2), tsMat%*%hugeMat[,c(8,11,12,13,9,14,15,16,10)])
  ddlambda[,7] <- ddlambda[,7] + 2*beta0/((1-beta)^3)
  
  temp1 <- (tsRel/lambdas - 1)%*%ddlambda
  temp1 <- rbind(c(temp1[1],0, temp1[2],0), cbind(c(0, temp1[2],0), matrix(temp1[3:11], ncol = 3)))
  
  temp2 <- t(apply(dlambda, 1, function(x) outer(x, x)))
  temp2 <- (tsRel/lambdas^2)%*%temp2
  temp2 <- matrix(temp2, ncol = 4)
  
  res <- -solve(temp1 - temp2)
  
  return(res)
}

INFIGARCH11sigma <- function(tsRel, tsMat, beta0, alpha, beta, d, lagMax){
  n <- length(tsRel) + lagMax
  hugeMat <- INFIGARCH11ddpsi(alpha = alpha, beta = beta, d = d, lagMax = lagMax)
  lambdas <- c(tsMat%*%hugeMat[,2] + beta0/(1 - beta))
  dlambda <- cbind(1/(1-beta), tsMat%*%hugeMat[,c(3,4,6)])
  dlambda[,3] <- dlambda[,3] + beta0/((1 - beta)^2)
  
  dl <- (tsRel/lambdas - 1)%*%dlambda
  
  ddlambda <- cbind(0, 1/((1-beta)^2), tsMat%*%hugeMat[,c(8,11,12,13,9,14,15,16,10)])
  ddlambda[,7] <- ddlambda[,7] + 2*beta0/((1-beta)^3)
  
  temp1 <- (tsRel/lambdas - 1)%*%ddlambda
  temp1 <- rbind(c(temp1[1], 0, temp1[2], 0), cbind(c(0, temp1[2], 0), matrix(temp1[3:11], ncol = 3, byrow = T)))
  
  
  temp2 <- t(apply(dlambda, 1, function(x) outer(x, x)))
  temp2 <- (tsRel/lambdas^2)%*%temp2
  temp2 <- matrix(temp2, ncol = 4)
  
  A <- (temp1 - temp2)
  temp3 <- (tsRel/lambdas - 1)*dlambda
  B <- (t(temp3)%*%temp3)
  res <- solve(A)%*%B%*%solve(A)
  
  return(res)
}
