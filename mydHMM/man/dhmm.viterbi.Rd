\name{dhmm.viterbi}
\alias{dhmm.viterbi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the most prossible path of states
}
\description{
The Viterbi algorithm is a dynamic programming algorithm for computing the most likely sequence given the parameters of Hidden Markov Model consist of number of states, transition matrix, parameters of normal distributions
}
\usage{
dhmm.viterbi(model, obser)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A Hidden Markov Model consist of number of states, transition matrix, parameters of normal distributions
}
  \item{obser}{
A sequence of observations
}
}

\value{
A vector of strings,containing the most probable path of states
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
dhmm.viterbi(model, obser)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

