rINGARCH <-
function(n, param = list(intercept = 1, past_obs = NULL, past_mean = NULL),
         model = list(past_obs = NULL, past_mean = NULL),
         link = c("identity", "log"),
         distr = c("poisson", "nbinom"),
         distrcoefs = NULL,
         n_start = 1000){
  if(length(link) > 1) link <- link[1]
  if(length(distr) > 1) distr <- distr[1]
  
  lambda <- ingarch.mean(intercept = param$intercept, past_obs = param$past_obs,
                                   past_mean = param$past_mean)
  
  if(!is.null(model$past_mean)){
    betas <- numeric(max(model$past_mean) + 1)
    betas[model$past_mean + 1] <- param$past_mean
  }
  else betas <- numeric(1)
  betas[1] <- param$intercept
  
  if(!is.null(model$past_obs)){
    alphas <- numeric(max(model$past_obs))
    alphas[model$past_obs] <- param$past_obs
  }
  else alphas <- 0
  
  if(link == "identity"){
    if(distr == "poisson"){
      res <- rINGARCHcppIDpois(n = n, betas = betas, alphas = alphas, lambda = lambda,
                               n_start = n_start)
    }
    
    if(distr == "nbinom"){
      res <- rINGARCHcppIDnbinom(n = n, betas = betas, alphas = alphas, lambda = lambda,
                                 n_start = n_start, size = distrcoefs)
    }
  }
  
  if(link == "log"){
    if(distr == "poisson"){
      res <- rINGARCHcppLOGpois(n = n, betas = betas, alphas = alphas, lambda = lambda,
                                n_start = n_start)
    }
    
    if(distr == "nbinom"){
      res <- rINGARCHcppLOGnbinom(n = n, betas = betas, alphas = alphas, lambda = lambda,
                                  n_start = n_start, size = distrcoefs)
    }
  }
  
  return(res)
}
