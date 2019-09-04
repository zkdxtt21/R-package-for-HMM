\name{dhmm.model}
\alias{dhmm.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
genuate random parameters of dhmm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dhmm.model(lambda, mc_type = "eq_full")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{mc_type}{
%%     ~~Describe \code{mc_type} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (lambda, mc_type = "eq_full")
{
    states <- seq(1, lambda, length.out = lambda)
    if (mc_type == "eq_full") {
        transp <- matrix(rep(1/lambda, lambda^2), nrow = lambda)
    }
    else if (mc_type == "rd_full") {
        transp <- matrix(rep(1/lambda, lambda^2), nrow = lambda)
    }
    else if (mc_type == "count") {
        transp <- matrix(rep(1/lambda, lambda^2), nrow = lambda)
    }
    else {
        transp <- matrix(rep(1/lambda, lambda^2), nrow = lambda)
    }
    set.seed(123)
    mu <- round((lambda * runif(lambda)/2), 2)
    set.seed(234)
    sd <- round(runif(lambda), 2)
    res <- list(states, transp, mu, sd)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
