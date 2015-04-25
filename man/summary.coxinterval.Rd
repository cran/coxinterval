\name{summary.coxinterval}
\alias{summary.coxinterval}
\alias{print.summary.coxinterval}

\title{Summary method for Cox model fit}

\description{
  \code{summary} method for the class \code{"coxinterval"}.
}

\usage{
\method{summary}{coxinterval}(object, conf.int = 0.95, scale = 1, ...)
}

\arguments{
  \item{object}{
    an object returned by \code{\link{coxaalen}} or \code{\link{coxdual}}.
  }
  \item{conf.int}{
    level for confidence intervals. If \code{FALSE}, no confidence
    intervals are provided.
  }
  \item{scale}{
    scale factor for the confidence intervals, whose limits represent
    the change in risk associated with one scale unit increase in the
    corresponding variable.
  }
  \item{\dots}{
    further arguments for other methods.
  }
}

\value{
  An object of the class \code{"summary.coxinterval"}.
}

\examples{
# Fit Cox model to dual-right--censored data
fit <- coxdual(Surv(start, stop, status) ~ cluster(id) + trans(from, to)
               + I(z * (to == 1)) + I(z * (from \%in\% 0 & to == 2))
               + I(z* (from \%in\% c(NA, 1) & to == 2)), data = dualrc,
               sieve.rate = 2/5)
fit
}

\keyword{survival}
