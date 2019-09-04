\name{dhmm.estimate}
\alias{dhmm.estimate}

\title{
Estemating Hidden Markov Model given number of states
}
\description{
Given as input initial HMM parameters consist of number of states, transition matrix, parameters of normal distribution and a sequence of observations O = o1,o2,...,oT , find the best parameters of HMM given the number of states.
}
\usage{
dhmm.estimate(ini_dhmm, obser, epsilon = 0.0001, max_iter = 1000)
}
\arguments{
  \item{ini_dhmm}{
inistial HMM paramters consist of number of states, transition matrix, parameters of normal distributions
}
  \item{obser}{
A sequence of observations
}
  \item{epsilon}{
precision of log-likelihood
}
  \item{max_iter}{
the maximum number of iterations
}
}

\value{

\item{ini_dhmm}{Best estimate of all parameters of hiddem markov model given number of states}
\item{log_likeli}{a sequence of log-likihood value in a full iteration}

}
\references{
http://www.cogsci.ucsd.edu/~ajyu/Teaching/Tutorials/hmm.pdf
https://www.cs.sjsu.edu/~stamp/RUA/HMM.pdfyh
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
size <- 400
epsilon <- 0.00001
max_iter <- 1000

model <- dhmm.model(2, 'count')
model[[3]] = c(0,3.5)
model[[2]] = matrix(c(0.2,0.8,0.3,0.7), 2, 2, byrow = T)

sim_samples <- dhmm.sim(model, size)
obser <- sim_samples[[2]]

ini_dhmm <- list(c(1,2), matrix(c(0.5, 0.5, 0.5, 0.5), nrow = 2), c(1, 3), c(0.5, 0.5))
result <- dhmm.estimate(ini_dhmm, obser, epsilon, max_iter)

}

