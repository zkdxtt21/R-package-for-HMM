\name{dhmm.decoding}
\alias{dhmm.decoding}

\title{
Hidden Markov Chain decoding
}
\description{
Given as input a sequence of observations O = o1,o2,...,oT and possible number of states, find the most probable sequence of states
}
\usage{
dhmm.decoding(obser, range = c(2,3,4,5), epsilon = 0.0001, max_iter = 1000)
}

\arguments{
  \item{obser}{
A sequence of observations
}
  \item{range}{
Possible number of states, default = 3
}
  \item{epsilon}{
precision of log-likelihood
}
  \item{max_iter}{
the maximum number of iterations
}
}

\value{

%%  If it is a LIST, use
\item{best_model}{best estimate of all parameters of hiddem markov model}
\item{best_seq}{best estimate of hidden states}
\item{likeli_list}{each element in the list represent the sequence of log-likihood value in a full iteration }

}

\references{
http://www.cogsci.ucsd.edu/~ajyu/Teaching/Tutorials/hmm.pdf
https://www.cs.sjsu.edu/~stamp/RUA/HMM.pdf
}

\examples{
size <- 400
epsilon <- 0.00001
max_iter <- 1000

model <- dhmm.model(2, 'count')
model[[3]] = c(0,3.5)
model[[2]] = matrix(c(0.2,0.8,0.3,0.7), 2, 2, byrow = T)

sim_samples <- dhmm.sim(model, size)
obser <- sim_samples[[2]]
range <- c(2, 3, 4)

result <- dhmm.decoding(obser, range, epsilon, max_iter)
plot(obser, pch=sim_samples[[1]], col=result[[2]])

}

