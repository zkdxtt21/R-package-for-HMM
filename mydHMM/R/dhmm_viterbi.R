

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
  scaling <- rep(NA, time_stp)
  
  prob_mtx[, 1] <- colSums(transp)/sum(colSums(transp)) * 
                   dnorm(rep(obser[1], state_num), mean=mu, sd=sd)
  scaling[1] = sum(prob_mtx[, 1])
  prob_mtx[, 1] = prob_mtx[, 1]/scaling[1]
  best_path[, 1] <- 0
  
  for (i in 2:time_stp) {
    for (j in 1:state_num) {
      
      tmp <- prob_mtx[, i-1] * transp[, j] * 
             dnorm(obser[i], mean=mu[j], sd=sd[j])
      
      prob_mtx[j, i] <- max(tmp)
      best_path[j, i] <- which.max(tmp)
    }
    scaling[i] = sum(prob_mtx[, i])
    prob_mtx[, i] = prob_mtx[, i]/scaling[i]
  }
  
  best <- rep(NA, time_stp)
  best[time_stp] <- which.max(prob_mtx[, time_stp])
  
  for (i in time_stp:2) {
    best[i-1] <- best_path[best[i], i]
  }

  # compute the log-likelihood  w
#  log_scale_total = 1
#  for(i in 1:state_num) {
#    log_scale_total = log_scale_total + log(scaling[i])
#  }
#  llk <- log_scale_total + log(max(prob_mtx[, time_stp]), base = exp(1))
#  print(llk)
    
  res <- list(prob_mtx, best)
  return (res)
}

