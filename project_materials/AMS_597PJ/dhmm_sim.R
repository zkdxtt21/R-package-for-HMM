
### Generate a DHMM Model
dhmm.model <- function(lambda, mc_type = 'eq_full') {
  
  states <- seq(1, lambda,  length.out = lambda)
  
  if (mc_type == 'eq_full') { transp <- matrix(rep(1 / lambda, lambda^2), nrow = lambda) }
  else if (mc_type == 'rd_full') { transp <- matrix(rep(1 / lambda, lambda^2), nrow = lambda) }
  else if (mc_type == 'count') {
    
    transp <- matrix(rep(1 / lambda, lambda^2), nrow = lambda) 
    
  }
  else {
    transp <- matrix(rep(1 / lambda, lambda^2), nrow = lambda) 
  }
  
  set.seed(123)
  mu <- round((lambda * runif(lambda) / 2), 2)
  set.seed(234)
  sd <- round(runif(lambda), 2)
  
  res <- list(states, transp, mu, sd)
  return (res)
}


### samples DHMM
dhmm.sim <- function(model, sample_size = 1) {
  
  dhmm <- model
  # get states, transp, mu, var
  states <- dhmm[[1]]
  transp <- dhmm[[2]]
  mu <- dhmm[[3]]
  sd <- dhmm[[4]]
  
  # return var augment
  samples <- rep(NA, sample_size)
  obser <- rep(NA, sample_size)
  
  # initialize uniformly
  samples[1] <- sample(states, 1, replace = TRUE)
  obser[1] <- rnorm(1, mu[samples[1]], sd[samples[1]])
  
  for (i in 2:sample_size) {
    
    samples[i] <- sample(states, 1, prob = transp[samples[i-1], ], replace = TRUE)
    obser[i] <- rnorm(1, mu[samples[i]], sd[samples[i]])
  }  
  
  res <- list(samples, obser)
  return (res)
}


dhmm.sim(dhmm.model(5, 'eq_full'), 10)




