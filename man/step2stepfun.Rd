\name{step2stepfun}
\alias{step2stepfun}

\title{Step function from cumulative increments over time}

\description{
  A utility function for that returns a step function given cumulative
  increments over time.
}

\usage{step2stepfun(step, time = "time", stratum = NULL)}

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

\value{
  A function that returns cumulative increments for given times.  If
  stratified, a list of such functions given in the same order as the strata.
}

\seealso{
  \code{\link{jump2step}}
  \code{\link{step2jump}}
}

\examples{
fit <- coxdual(Surv(start, stop, status) ~ cluster(id)
               + trans(from, to) + z, data = dualrc, init.coxph = TRUE)
head(fit$coxph$basehaz)
Hazfun <- step2stepfun(fit$coxph$basehaz, stratum = 3)
Haz01 <- with(fit$coxph, split(basehaz[, 1:2], basehaz[, 3]))[[1]]
all(Hazfun[[1]](Haz01$time) == Haz01$hazard)
}

\keyword{survival}
