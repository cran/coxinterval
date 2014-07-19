\name{cosmesis}
\alias{cosmesis}
\docType{data}

\title{Breast cosmesis data}

\description{
  Interval-censored times to cosmetic deterioration for breast cancer
  patients undergoing radiation or radiation plus chemotherapy.
}

\usage{data(cosmesis)}

\format{
  A data frame with 94 observations on the following 3 variables.
  \describe{
    \item{\code{left}}{
      left endpoint of the censoring interval in months
    }
    \item{\code{right}}{
      right endpoint of the censoring interval in months
    }
    \item{\code{treat}}{
      factor with levels \code{RT} and \code{RCT} representing
      radiotherapy-only and radiation plus chemotherapy treatments,
      respectively
    }
  }
}

\source{
  Finkelstein, D. M. and Wolfe, R. A. (1985)
  \href{http://www.jstor.org/stable/2530698}{
    A semiparametric model for regression analysis of interval-censored
    failure time data}.
  \emph{Biometrics} \bold{41}, 933--945.
}

\references{
  Finkelstein, D. M. (1986)
  \href{http://www.jstor.org/stable/2530698}{
    A proportional hazards model for interval-censored failure time data}.
  \emph{Biometrics} \bold{42}, 845--854.
}

\examples{data(cosmesis)}

% RShowDoc("KEYWORDS")
\keyword{datasets}
