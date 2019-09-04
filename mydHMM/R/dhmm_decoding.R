

# inference hidden states sequence given a sequence of observations
dhmm.decoding <- function(obser, range = c(2,3,4,5), epsilon = 0.0001, max_iter = 1000) {
  
  stopifnot(is.numeric(obser))
  stopifnot(all(range == floor(range)))
  stopifnot(all(range > 0))
  stopifnot(epsilon > 0)
  stopifnot(max_iter > 0)

  # range is a postive integer vector
  bic <- NULL
  model_list <- NULL
  best_list <- NULL

  time_stp <- length(obser)
  likeli_list <- list()

  for (i in 1:length(range)) {

    lambda <- range[i]

    # using kmeans method to get the initial mu and sd value
    ini_dhmm <- dhmm.model(lambda, 'eq_full')
    kmeans_iter = 100
    mu_ini <- matrix(rep(NA, kmeans_iter * lambda), nrow = kmeans_iter)
    sd_ini <- matrix(rep(NA, kmeans_iter * lambda), nrow = kmeans_iter)
    for(m in 1:kmeans_iter){
      km <- kmeans(obser, lambda, iter.max = 10)
      mu_ini[m,] <- c(sort(km$centers)) # why c() ???
      ord <- order(km$centers)
      sd_ini[m,] <- sqrt(km$withinss[ord] / (km$size[ord]-1))
    }
    ini_dhmm[[3]] <- apply(mu_ini, 2, mean)
    ini_dhmm[[4]] <- apply(sd_ini, 2, mean)

    # EM to estimate model
    res <- dhmm.estimate(ini_dhmm, obser, epsilon, max_iter)
    dhmm_model <- res[[1]]
    # viterbi algorithm to find best path
    best_gauss <- dhmm.viterbi(dhmm_model, obser)

    ## --------------- BIC compute -----------------##
    likeli_list[[i]] <- res[[2]]
    max_iter <- length(likeli_list[[i]])
    para_num <- lambda * (lambda + 1)
    bic[i] <- (-2) * likeli_list[[i]][max_iter] + log(time_stp) * para_num
    print(bic[i])

    # save data to list
    model_list[[i]] <- dhmm_model
    best_list <- cbind(best_list, best_gauss[[2]])

  }
  best_model <- model_list[[which(bic == min(bic))]]
  names(best_model) <- c("states", "transition matrix", "mean", "sd")
  best_seq <- best_list[, which(bic == min(bic))]
  
  res <- list(best_model, best_seq, likeli_list)
  names(res) <- c("best_model", "best_seq", "likeli_list")
  return (res)
}



