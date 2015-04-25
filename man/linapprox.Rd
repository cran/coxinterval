\name{linapprox}
\alias{linapprox}

\title{Linear approximation}

\description{Perform linear interpolation or approximation}

\usage{linapprox(xyin, xout)}

\arguments{
  \item{xyin}{
    matrix whose rows give the coordinate pairs of the points to be
    interpolated or extrapolated.
  }
  \item{xout}{
    a vector of numeric values at which interpolation or extrapolation
    should take place.
  }
}

\value{A vector of interpolated or extrapolated values.}

\seealso{
  \code{\link[stats]{approx}}
  \code{\link{const2lin}}
  \code{\link{lin2const}}
}
\examples{
fit <- coxdual(Surv(start, stop, status) ~ cluster(id)
               + trans(from, to) + z, data = dualrc, init.coxph = TRUE)
head(basehaz)
Haz01 <- with(fit, split(basehaz[, 1:2], basehaz[, 3]))[[1]]
all(Haz01$hazard == with(Haz01, linapprox(cbind(time, hazard), time)))
}
