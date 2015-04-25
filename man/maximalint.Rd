\name{maximalint}
\alias{maximalint}

\title{Find maximal intersections from an interval-type survival object}

\description{
  A utility function that returns the "maximal intersections" from a set
  of censoring intervals in a given \code{type = "interval"} or
  \code{type = "interval2"} \code{\link[survival]{Surv}} object.
}

\usage{maximalint(x, eps = 1e-7)}

\arguments{
  \item{x}{
    a two-column matrix or data frame giving the left- and
    right-endpoints of the censoring intervals or, alternatively, a
    \code{type = "interval"} or \code{type = "interval2"}
    \code{\link[survival]{Surv}} object.  Following the
    \code{type = "interval2"} format, any \code{NA}-valued right
    endpoints are considered the same as \code{Inf}.
  }
  \item{eps}{
    a small value used to break ties in maximal intersection endpoints.
    Note that a large time scale requires a larger epsilon in order to
    ensure that \code{t != t + eps} for a given time value \code{t}.
  }
}

\value{
  \item{int}{
    a two-column matrix whose rows give the maximal
    intersections (Wong and Yu, 1999) of the censoring intervals
    represented in \code{x}
  }
  \item{ind}{
    a \code{nrow(x)} by \code{nrow(int)} matrix that indicates
    overlap between each maximal intersection and censoring interval
  }
}

\details{
  Each censoring interval is assumed to exclude its left
  endpoint.  Analogous to Maathuis (2005, Section 2.1) ties between the
  unique (open) left and (closed) right endpoints among the censoring
  intervals are broken by subtracting a small value from the right
  endpoint value before constructing the maximal intersections.
}

\references{
  Maathuis, MH (2005)
  \href{http://dx.doi.org/10.1198/106186005X48470}{
    Reduction algorithm for the NPMLE for the distribution function
    of bivariate interval-censored data}.
  \emph{Journal of Computational and Graphical Statistics}
  \bold{14}, 352--362.

  Wong, GYC. and Yu, Q. (1999)
  \href{http://dx.doi.org/10.1006/jmva.1998.1807}{
    Generalized MLE of a joint distribution function with multivariate
    interval-censored data}.
  \emph{Journal of Multivariate Analysis} \bold{69}, 155--166.
}

\examples{
s <- with(cosmesis[1:10, ], Surv(left, right, type = "interval2"))
s
maximalint(s)
}

\keyword{survival}
