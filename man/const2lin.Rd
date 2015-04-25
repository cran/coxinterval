\name{const2lin}
\alias{const2lin}

\title{Accumulate a piecewise constant function of time}

\description{
  A utility function for \code{\link{coxdual}} that returns the integrates
  of a piecewise constant function of time, evaluated at given values.
}

\usage{const2lin(const, time = "time", stratum = NULL)}

\arguments{
  \item{const}{
    A matrix whose columns give the (possibly multivariate) piecewise
    constant function values, corresponding time points and strata.
  }
  \item{time}{
    an integer or character value indicating the column index or name of
    the time variable in the matrix \code{const}.  Defaults to \code{"time"}.
  }
  \item{stratum}{
    integer or character value for the column index or name of the
    stratifying variable in the matrix \code{const}.  If \code{stratum}
    is \code{NULL} (default), \code{const} is presumed unstratified.
  }
}

\value{A matrix of the same dimension as \code{const}.}

\seealso{\code{lin2const}}

\examples{
fit <- coxdual(Surv(start, stop, status) ~ cluster(id)
               + trans(from, to) + z, data = dualrc, init.coxph = TRUE)
haz <- lin2const(fit$basehaz, stratum = 3)
Haz <- const2lin(haz, stratum = 3)
all(Haz == fit$basehaz)
}

\keyword{survival}
