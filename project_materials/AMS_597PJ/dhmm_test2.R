############# model 1 #################
size <- 200
model <- dhmm.model(2, 'count')
model[[3]] = c(0,3.5)
model[[2]] = matrix(c(0.2,0.8,0.3,0.7),2,2,byrow = T)

############# model 2 ######££££#######
size <- 600
model <- dhmm.model(4, 'count')
#ini_dhmm = dhmm.model(3, 'count')
#ini_dhmm[[3]] = c( 5,0.5, 5)
model[[3]] = c(0,3,5,7)
model[[2]] <- matrix(c(0.3,0.2,0.3,0.2,0.2,0.2,0.5,0.1,0.3,0.3,0.1,0.3,0.2,0.2,0.3,0.3),4,4, byrow = T)






sim_samples <- dhmm.sim(model, size)
obs = sim_samples[[2]]
obser = obs


range <- c(2,3,4)
max_iter <- 100




result <- dhmm.decoding(obser,4, max_iter)
plot(obs, pch=sim_samples[[1]], col=result[[2]])




km = kmeans(obser, 3)
1-sum(km$withinss)/km$totss
plot(obs, pch=sim_samples[[1]], col=km$cluster)
