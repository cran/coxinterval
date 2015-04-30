\name{jump2step}
\alias{jump2step}

\title{Accumulate increments in values from a step function of time}

\description{
  A utility function for \code{\link{coxdual}} that returns values from
  a step function given its increments over time.
}

\usage{jump2step(jump, time = "time", stratum = NULL)}

\arguments{
  \item{jump}{
    a matrix whose columns give increments in the (possibly
    multivariate) step function, corresponding time points and strata.
  }
  \item{time}{
    an integer or character value indicating the column index or name of
    the time variable in the matrix \code{jump}.  Defaults to \code{"time"}.
  }
  \item{stratum}{
    integer or character value for the column index or name of the
    stratifying variable in the matrix \code{jump}.  If \code{stratum}
    is \code{NULL} (default), \code{jump} is presumed unstratified.
  }
}

\value{
  A matrix of the same dimension as \code{jump}.
}

\seealso{
  \code{\link{step2jump}}
}

\examples{
fit <- coxdual(Surv(start, stop, status) ~ cluster(id)
               + trans(from, to) + z, data = dualrc, init.coxph = TRUE)
haz <- step2jump(fit$coxph$basehaz, stratum = 3)
Haz <- jump2step(haz, stratum = 3)
all(Haz == fit$coxph$basehaz)
}

\keyword{survival}
