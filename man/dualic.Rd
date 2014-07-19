\name{dualic}
\alias{dualic}
\docType{data}

\title{Simulated dual-censored data from an illness-death process}

\description{
  Data from a Markov illness-death process with interval-censored
  progression times, simulated according to the initial scenario
  described in Boruvka and Cook (2014b).
}

\usage{data(dualic)}

\format{
  A data frame with 723 observations on the following 7 variables.
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
      left-endpoint of the time interval at which the subject is known
      to be at risk for a transition between state \code{from} and state
      \code{to}.
    }
    \item{\code{stop}}{
      right-endpoint of the at-risk interval.
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
  Boruvka, A. and Cook, R. J. (2014b)
  Sieve estimation in a Markov illness-death process under dual censoring.
}

\examples{data(dualic)}

\seealso{
  \code{\link{coxic}}
}

\keyword{datasets}
