
source('AMS_597PJ/dhmm_sim.R')
source('AMS_597PJ/dhmm_estimate.R')
source('AMS_597PJ/dhmm_viterbi.R')

# inference hidden states sequence given a sequence of observations
dhmm.decoding <- function(obser, range = 5, max_iter) {
  
  # range can only take one or two dim vector
  
  bic <- NULL
  model_list <- NULL
  best_list <- NULL

  
  
  for (i in 1:length(range)) {
    
    lambda <- range[i]
    
    # using kmeans method to get the initial value
    ini_dhmm <- dhmm.model(lambda, 'eq_full')
    kmeans_iter = 100
    mu_ini <- matrix(rep(NA, kmeans_iter * lambda), nrow = kmeans_iter)
    sd_ini <- matrix(rep(NA, kmeans_iter * lambda), nrow = kmeans_iter)
    for(m in 1:kmeans_iter){
      km <- kmeans(obser, lambda, iter.max = 10)
      mu_ini[m,] <- c(sort(km$centers))
      ord <- order(km$centers)
      sd_ini[m,] <- sqrt(km$withinss[ord]/(km$size[ord]-1))
    }
    ini_dhmm[[3]] <- apply(mu_ini, 2, mean)
    ini_dhmm[[4]] <- apply(sd_ini, 2, mean)
    
    # em to estimate model
    dhmm_model <- dhmm.estimate(ini_dhmm, obser, max_iter)
    # viterbi algorithm to find best path
    best_gauss <- dhmm.viterbi(dhmm_model, obser)
    
    # compute BIC...
#    bic[i] = -2*best_gauss[[3]] + log(time_stp) * lambda
#    print(bic[i])
  
    # save data to list
    model_list[[i]] <- dhmm_model
    best_list <- cbind(best_list, best_gauss[[2]])

  }
  
  # best_model <- model_list[[which(BIC == min(BIC))]]
  # best_state <- best_list[, which(BIC == min(BIC))]
  # best <- list(best_model, best_state)
  
  # ------------ temp, test code ----------------
  bestv <- list(dhmm_model, best_gauss[[2]])
  
  return (bestv)
}



