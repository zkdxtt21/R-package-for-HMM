
# compute backward procedure given model
dhmm.backward <- function(model, obser, scaling) {
  
  # get states, transp, mu, sd
  states <- model[[1]]
  transp <- model[[2]]
  mu <- model[[3]]
  sd <- model[[4]]
  
  time_stp <- length(obser)
  state_num <- length(states)
  
  # create a backward table
  backward <- matrix(rep(NA, time_stp * state_num), nrow = state_num)
  
  # initialize 
  backward[, time_stp] <- rep(1, state_num) * scaling[time_stp]
  
  for (i in time_stp:2) {
    
    tmp <- backward[, i] * dnorm(rep(obser[i], state_num), mean=mu, sd=sd)
    backward[, i-1] <- transp %*% tmp * scaling[i-1]
  }
  
  return (backward)
}
