\name{lin2const}
\alias{lin2const}

\title{Slopes from values of a piecewise linear function of time}

\description{
  A utility function for \code{\link{coxdual}} that gives the slope
  corresponding to values from a piecewise linear function of time.
}

\usage{lin2const(lin, time = "time", stratum = NULL)}

\arguments{
  \item{lin}{
    a matrix whose columns give the (possibly multivariate) piecewise
    linear function values, time points and strata.
  }
  \item{time}{
    an integer or character value indicating the column index or name of
    the time variable in the matrix \code{lin}.  Defaults to \code{"time"}.
  }
  \item{stratum}{
    integer or character value for the column index or name of the
    stratifying variable in the matrix \code{lin}.  If \code{stratum}
    is \code{NULL} (default), \code{lin} is presumed unstratified.
  }
}

\value{
  A matrix of the same dimension as \code{lin}.
}

\details{
  This is a utility function for \code{\link{coxdual}}.
}

\seealso{
  \code{const2lin}
}

\examples{
data(dualrc)
fit <- coxdual(Surv(start, stop, status) ~ cluster(id)
               + trans(from, to) + z, data = dualrc, init.coxph = TRUE)
fit$basehaz
haz <- lin2const(fit$basehaz, stratum = 3)
Haz <- const2lin(haz, stratum = 3)
all(Haz == fit$basehaz)
}

\keyword{survival}
