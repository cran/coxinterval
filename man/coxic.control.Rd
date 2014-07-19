\name{coxic.control}
\alias{coxic.control}

\title{Control Cox model fit}

\description{
  Set parameters controlling the model fit returned by
  \code{\link{coxic}}.
}

\usage{
coxic.control(eps = 1e-07, iter.max = 50000, coef.typ = 1,
              coef.max = 10, sieve.const = 1, sieve.rate = 1/3)
}

\arguments{
  \item{eps}{
    maximum threshold for the relative change in the model fit required
    to meet the iterative estimation routine's stopping rule.
  }
  \item{iter.max}{
    maximum number of iterations to attempt.  This ensures that
    \code{\link{coxic}} will eventually exit, even when there are
    convergence issues.
  }
  \item{coef.typ}{
    a scalar or vector of typical (absolute) values for the
    regression coefficient.
  }
  \item{coef.max}{
    a scalar or vector of probable upper bounds for the regression
    coefficient.  This and the \code{coef.typ} arguments tune
    variance estimation via the curvature in the profile
    log-likelihood.
  }
  \item{sieve.const}{
    a constant factor that, in part, determines the sieve size.  The
    factor can be made specific to the transition type with
    \code{sieve.const} a vector of length three.  Indexing the states
    from zero, this vector's components correspond to the state 0 to
    state 1, 0 to 2, and 1 to 2 transition types, respectively.
  }
  \item{sieve.rate}{
    a scalar in (1/8, 1/2) determining the rate at which the sieve
    increases with the sample size.
  }
}

\value{
  A list of the above arguments with their final values.
}

\details{
  For a given sample size \emph{n}, the resulting sieve has size
  at most \code{sieve.const*}\emph{n}\code{^sieve.rate}.  Any reduction
  in size from this value is applied to ensure that each subinterval in
  the sieve's time partition captures at least one support point from
  the semiparametric maximum likelihood estimator based on the subsample
  with known progression status (Boruvka and Cook, 2014b).
}

\references{
  Boruvka, A. and Cook, R. J. (2014b)
  Sieve estimation in a Markov illness-death process under dual censoring.
}

\seealso{\code{\link{coxic}}}

\examples{
coxic(Surv(start, stop, status) ~ cluster(id) + trans(from, to)
      + I(z * (to == 1)) + I(z * (from \%in\% 0 & to == 2))
      + I(z * (from \%in\% c(NA, 1) & to == 2)), data = dualrc,
      control = coxic.control(eps = 1e-5, sieve.rate = 2/5))
}

\keyword{survival}
