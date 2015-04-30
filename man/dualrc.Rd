\name{dualrc}
\alias{dualrc}
\docType{data}

\title{Simulated dual-right--censored data from an illness-death process}

\description{
  Data from a Markov illness-death process with progression and death
  observed under dual right censoring, simulated according to the
  initial scenario described in Boruvka and Cook (2014).
}

\usage{data(dualrc)}

\format{
  A data frame with 644 observations on the following 7 variables.
  \describe{
    \item{\code{id}}{
      subject identifier.
    }
    \item{\code{from}}{
      originating state index with \code{0} denoting the initial state,
      \code{1} the intermediate state, \code{2} the terminal state and
      \code{NA} an unknown state.
    }
    \item{\code{to}}{
      subsequent state index.
    }
    \item{\code{start}}{
      left endpoint of the time interval at which the subject is known
      to be at risk for a transition between state \code{from} and state
      \code{to}.
    }
    \item{\code{stop}}{
      right endpoint of the at-risk interval.
    }
    \item{\code{status}}{
      indicator that a transition between state \code{from} and state
      \code{to} was observed at \code{stop}.
    }
    \item{\code{z}}{
      a binary covariate.
    }
  }
}

\references{
  Boruvka, A. and Cook, R. J. (2014)
  Sieve estimation in a Markov illness-death process under dual censoring.
}

\seealso{
  \code{\link{coxdual}}
}

\examples{data(dualrc)}

\keyword{datasets}
