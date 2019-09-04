\name{dhmm.sim}

\title{
Simulate states and observations for a Hidden Markov Model
}

\description{
Simulates a path of states and observations for a given Hidden Markov Model.
}

\usage{
dhmm.sim(model, sample_size)
}

\arguments{
  \item{model}{
		A hidden Markov Model consist of number of states, transition matrix, parameters of normal distributions
	}
  \item{sample_size}{
		The length of the simulated sequence of observations and states
	}
}

\value{
The function dhmm.sim returns a path of states and associated observations:

  \item{states}{The path of states}
  \item{obser}{The sequence of obeservations}

}
\references{
	http://www.cogsci.ucsd.edu/~ajyu/Teaching/Tutorials/hmm.pdf
	https://www.cs.sjsu.edu/~stamp/RUA/HMM.pdf

}
\examples{
# initialise HMM
model <- dhmm.model(4, 'count')
# simulate from the HMM
size <- 400
sim_samples <- dhmm.sim(model, size)
}
