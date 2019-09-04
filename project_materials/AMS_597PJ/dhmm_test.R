
source('AMS_597PJ/dhmm_sim.R')
source('AMS_597PJ/dhmm_estimate.R')
source('AMS_597PJ/dhmm_decoding.R')

source('AMS_597PJ/dhmm_viterbi.R')

#### -------------------- THIS IS A TEST SCRIPT -------------------

## --------------------- generate DHMM -----------------------
size <- 200
model <- dhmm.model(3, 'count')
sim_samples <- dhmm.sim(model, size)

true_states <- sim_samples[[1]]
obser <-  sim_samples[[2]]
transp <- model[[2]]
mu <- model[[3]]
sd <- model[[4]]

mu_col <- mu[true_states]
sd_col <- sd[true_states]

mydata <- data.frame(true_states, obser, mu_col, sd_col)
remove(sim_samples)
remove(mu_col)
remove(sd_col)

## --------------------- DHMM decoding -----------------------

range <- 3
max_iter <- 1000

 best <- dhmm.decoding(obser, range, max_iter)
 error_rate <- sum(best[[2]] != true_states) / size

res <- dhmm.viterbi(model, obser)
  
