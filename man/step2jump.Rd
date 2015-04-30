\name{step2jump}
\alias{step2jump}

\title{Increments from a step function of time}

\description{
  A utility function for \code{\link{coxdual}} that returns increments
  from values of a step function over time.
}

\usage{step2jump(step, time = "time", stratum = NULL)}

\arguments{
  \item{step}{
    a matrix whose columns give the (possibly multivariate) step
    function values, corresponding time points and strata.
  }
  \item{time}{
    an integer or character value indicating the column index or name of
    the time variable in the matrix \code{step}.  Defaults to \code{"time"}.
  }
  \item{stratum}{
    integer or character value for the column index or name of the
    stratifying variable in the matrix \code{step}.  If \code{stratum}
    is \code{NULL} (default), \code{step} is presumed unstratified.
  }
}

\value{A matrix of the same dimension as \code{step}.}

\seealso{\code{\link{jump2step}}}

\examples{
fit <- coxdual(Surv(start, stop, status) ~ cluster(id)
               + trans(from, to) + z, data = dualrc, init.coxph = TRUE)
haz <- step2jump(fit$coxph$basehaz, stratum = 3)
Haz <- jump2step(haz, stratum = 3)
all(Haz == fit$coxph$basehaz)
}

\keyword{survival}
