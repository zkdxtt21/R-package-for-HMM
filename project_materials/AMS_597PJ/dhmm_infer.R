

source('dhmm_sim.R')

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
  
  # initialize 
  forward[,1] <- rep(1/state_num, state_num) * dnorm(rep(obser[1], state_num), mean=mu, sd=sd)
  
  for (i in 2:time_stp) {
    
    forward[,i] <- forward[,i-1] %*% 
                   transp * dnorm(rep(obser[i], state_num), mean=mu, sd=sd)
    
  }
  
  return (forward)
}


# compute backward procedure given model
dhmm.backward <- function(model, obser) {
  
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
  backward[,time_stp] <- rep(1, state_num)

  for (i in 1:(time_stp-1)) {
    
    tmp <- backward[,time_stp-i+1] * dnorm(rep(obser[time_stp-i+1], state_num), mean=mu, sd=sd)
    backward[,time_stp-i] <- transp %*% t(tmp)
    
  }
  
  return (backward)
}





# estimate the DHMM model given the obser & a initial model
dhmm.estimate <- function(obser, model, lambda, seed, max_iter) {
  # using Baum-Welch algorithm
  
  # intialize a HMM
  ini_dhmm <- model(lambda, seed)
  # get states, transp, mu, var
  states <- model[[1]]
  transp <- model[[2]]
  mu <- model[[3]]
  sd <- model[[4]]
  
  time_stp <- length(obser)
  state_num <- length(states)
  
  r_tbl <- matrix(rep(NA, time_stp * state_num), nrow = state_num)
  xi_tbl <- matrix(rep(NA,  state_num ^ 2), nrow = state_num)
  
  while (i < max_iter) {
    
    forward <- dhmm.forward(ini_dhmm, obser)
    backward <- dhmm.backward(ini_dhmm, obser)
    
    ## -------------  E-steps ------------------
    for (i in 1:time_stp) {
      
      tmp <- forward[, i] * backward[, i]
      r_tbl[, i] <- tmp / sum(tmp)
      
    }

    #################### Tingting 04/13 start ####################      
    xi <- list()
    for (i in 1: (T-1) ) {
      tmp <- backward[,i+1] * dnorm(rep(obser[i+1], state_num), mean=mu, sd=sd)
      xi_i <- diag(forward[,i]) %*% transp %*% diag(tmp)
      xi[[i]] <-  xi_i
    }
    #################### Tingting 04/13 end ####################    
  
    
    ## ------------- M-steps -------------------
    
    #################### Tingting 04/13 start ####################    
    xi_sum <-  matrix(rep(0, n*n), nrow = n)
    for(i in 1:(T-1)){
      xi_sum <- xi_sum + xi[[i]]
    }
    transp <- xi_sum / matrix(rep(rowSums(r_tbl), state_num), nrow = state_num)
    #################### Tingting 04/13 end ####################   
    
    
    # update mu, sd
    mu <- r_tbl %*% obser / rowSums(r_tbl)
    
    # need to be check !!!
    mu_mtx <-  matrix(rep(mu, time_stp), nrow = time_stp, byrow = TRUE)
    obser_mtx <- matrix(rep(obser, state_num), nrow = state_num, byrow = TRUE)
    var_mtx <- (obser_mtx - mu_mtx) ^ 2
    tmp3 <- r_tbl * var_mtx
    sd <- sqrt(rowSums(tmp3) / rowSums(r_tbl))
    
  }
  
  res <- list(transp, mu, sd)
  return (res)
}


### ------------------------------------------------------------------------------------


# given X and estimate model, get esitimate of hidden state sample
dhmm.viterbi <- function(model, obser) {
  
  # get states, transp, mu, var
  states <- model[[1]]
  transp <- model[[2]]
  mu <- model[[3]]
  sd <- model[[4]]
  
  time_stp <- length(obser)
  state_num <- length(states)
  
  # create a prob mtx and best_path matrix
  prob_mtx <- matrix(rep(NA, state_num * time_stp), nrow = state_num)
  best_path <- prob_mtx
  
  prob_mtx[, 1] <- rep(1/state_num, state_num) * 
                   dnorm(rep(obser[1], state_num), mean=mu, sd=sd)
  best_path[, 1] <- 0
  
  for (i in 2:time_stp) {
    for (j in 1:state_num) {
      
      tmp <- prob_mtx[, i-1] %*% transp[, j] * 
             dnorm(obser[i], mean=mu[j], sd=sd[j])
      
      prob_mtx[j, i] <- max(tmp)
      best_path[j, i] <- which.max(tmp)
        
    }
  }
  
  best <- rep(NA, time_stp)
  best[time_stp] <-  which.max(prob_mtx[, time_stp])
  
  for (i in time_stp:2) {
    best[i-1] <- which.max(best_path[best[i], i])
  }
  
  return (best)
}



# inference hidden states sequence given a sequence of observations
dhmm.inference <- function(obser) {
  
  # parameter initialize
  seed <- 123
  lambda <- 5
  max_iter <- 1000
  
  # em to estimate model
  dhmm_model <- dhmm.estimate(obser, model, lambda, seed, max_iter)
  # viterbi algorithm to find best path
  best_gauss <- dhmm.viterbi(dhmm_model, obser)
  
  return (best_gauss)
  
}



