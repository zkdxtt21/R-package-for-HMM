
# compute forward procedure given model
dhmm.forward <- function(model, obser) {
  
  # get states, transp, mu, sd
  states <- model[[1]]
  transp <- model[[2]]
  mu <- model[[3]]
  sd <- model[[4]]
  
  time_stp <- length(obser)
  state_num <- length(states)
  
  # create a forward table
  forward <- matrix(rep(NA, time_stp * state_num), nrow = state_num)
  
  # create a scaling vector
  scaling <- rep(NA, time_stp)
  
  # initialize 
  forward[,1] <- rep(1/state_num, state_num) * dnorm(rep(obser[1], state_num), mean=mu, sd=sd)
  scaling[1] <- 1 / sum(forward[,1])
  forward[,1] <- forward[,1] * scaling[1]
  
  for (i in 2:time_stp) {
    
    forward[,i] <- forward[,i-1] %*% 
                   transp * dnorm(rep(obser[i], state_num), mean=mu, sd=sd)
    scaling[i] <- 1 / sum(forward[,i])
    forward[,i] <- forward[,i] * scaling[i]
  }
  
  return ( list(forward, scaling) )
}

# compute forward procedure given model
dhmm.forward <- function(model, obser) {
  
  # get states, transp, mu, sd
  states <- model[[1]]
  transp <- model[[2]]
  mu <- model[[3]]
  sd <- model[[4]]
  
  time_stp <- length(obser)
  state_num <- length(states)
  
  # create a forward table
  forward <- matrix(rep(NA, time_stp * state_num), nrow = state_num)
  
  # create a scaling vector
  scaling <- rep(NA, time_stp)
  
  # initialize 
  forward[,1] <- rep(1/state_num, state_num) * dnorm(rep(obser[1], state_num), mean=mu, sd=sd)
  scaling[1] <- 1 / sum(forward[,1])
  forward[,1] <- forward[,1] * scaling[1]
  
  for (i in 2:time_stp) {
    
    forward[,i] <- forward[,i-1] %*% 
                   transp * dnorm(rep(obser[i], state_num), mean=mu, sd=sd)
    scaling[i] <- 1 / sum(forward[,i])
    forward[,i] <- forward[,i] * scaling[i]
  }
  
  return ( list(forward, scaling) )
}
