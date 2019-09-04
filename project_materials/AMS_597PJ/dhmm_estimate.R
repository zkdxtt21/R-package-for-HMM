
source('AMS_597PJ/dhmm_sim.R')
source('AMS_597PJ/dhmm_forward.R')
source('AMS_597PJ/dhmm_backward.R')

# estimate the DHMM model given the obser & a initial model
dhmm.estimate <- function(ini_dhmm, obser, max_iter) {
  # using Baum-Welch algorithm
  
#   transp <- ini_dhmm[[2]]
#   mu <- ini_dhmm[[3]]
#   sd <- ini_dhmm[[4]]
  
  time_stp <- length(obser)
  state_num <- length(ini_dhmm[[1]])
  
  for (i in 1:max_iter) {

    ####################
#    print(ini_dhmm)
    #####################
    
    
    
    forward <- dhmm.forward(ini_dhmm, obser)[[1]]
    scaling <- dhmm.forward(ini_dhmm, obser)[[2]]
    backward <- dhmm.backward(ini_dhmm, obser, scaling)
    
    transp <- ini_dhmm[[2]]
    mu <- ini_dhmm[[3]]
    sd <- ini_dhmm[[4]]
    

    ## -------------  E-steps ------------------
#    r_tbl <- matrix(rep(NA, state_num * time_stp), nrow = state_num)
#    for (i in 1:time_stp) {
      
#      tmp <- forward[, i] * backward[, i]
#      r_tbl[, i] <- tmp / sum(tmp)
#    }
    
    gamma_mtx <- matrix(rep(NA, state_num * time_stp), nrow = state_num)
    xi <- list()
    for (i in 1:(time_stp-1)) {
      
      tmp <- backward[,i+1] * dnorm(rep(obser[i+1], state_num), mean=mu, sd=sd)
      xi_i <- diag(forward[,i]) %*% transp %*% diag(tmp)
      xi[[i]] <- xi_i
      gamma_mtx[ ,i] = rowSums(xi_i)
    }
    gamma_mtx[,time_stp] = forward[,time_stp]
    
    
    ## ------------- M-steps -------------------
    
    # estimate transit matrix 
    xi_sum <- matrix(rep(0, state_num ^ 2), nrow = state_num)
    for (i in 1:(time_stp-1)) {
      
      xi_sum <- xi_sum + xi[[i]]

    }
    
    ini_dhmm[[2]] <- xi_sum / rowSums(gamma_mtx)
    
    # estimate mu, sd
    ini_dhmm[[3]] <- gamma_mtx %*% obser / rowSums(gamma_mtx)
    ini_dhmm[[3]] = c(ini_dhmm[[3]])
    mu_mtx <-  matrix(rep(mu, time_stp), nrow = state_num)
    obser_mtx <- matrix(rep(obser, state_num), nrow = state_num, byrow = TRUE)
    var_mtx <- (obser_mtx - mu_mtx) ^ 2
    ini_dhmm[[4]] <- sqrt(rowSums(gamma_mtx * var_mtx) / rowSums(gamma_mtx))
    
  }
  
  return (ini_dhmm)
}


