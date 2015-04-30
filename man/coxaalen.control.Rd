\name{coxaalen.control}
\alias{coxaalen.control}

\title{Control a Cox-Aalen model fit}

\description{
  Set parameters controlling the model fit returned by
  \code{\link{coxaalen}}.
}

\usage{
coxaalen.control(eps = 1e-07, eps.norm = c("max", "grad"),
                 iter.max = 5000, armijo = 1/3, var.coef = TRUE,
                 coef.typ = 1, coef.max = 10, trace = FALSE,
                 thread.max = 1, data = FALSE)
}

\arguments{
  \item{eps}{
    threshold value for the norm used to measure convergence in the
    parameter estimates.
  }
  \item{eps.norm}{
    a character string identifying the norm to use in the convergence
    criteria---either the maximum norm between the current and
    previous parameter values (\code{eps.norm = "max"}) or the absolute
    inner product between the current value and the score
    (\code{eps.norm = "grad"}).
  }
  \item{iter.max}{
    maximum number of iterations to attempt.  This ensures that
    \code{\link{coxaalen}} will eventually exit, even when the
    convergence criteria is not met.  A warning is issued whenever the
    estimation routine has stopped before converging on a final
    parameter value.
  }
  \item{armijo}{
    a scale factor in (0, 1/2) for Armijo's (1966) rule---a line search
    used to ensure that each iteration achieves an adequate increase in
    the log-likelihood.  The model fit is typically not very sensitive
    to this value.
  }
  \item{var.coef}{
    a logical value indicating that standard errors for the
    multiplicative regression coefficients should be estimated.  This is
    done via profile likelihood---an approach that can require an
    inordinate amount of processing time under many regression
    coefficients and larger sample size.
  }
  \item{coef.typ}{
    a scalar or vector of typical (absolute) values for the
    multiplicative regression coefficient.
  }
  \item{coef.max}{
    a scalar or vector of probable upper bounds for the multiplicative
    regression coefficient.  This and the \code{coef.typ} arguments
    tune variance estimation via the curvature in the profile
    log-likelihood.
  }
  \item{trace}{
    a logical value indicating that CPLEX should print its results to
    the screen.
  }
  \item{thread.max}{
    maximum number of CPU threads to allocate to CPLEX.  The default
    value disables multithreading.  A value of zero allows CPLEX to set
    the number of threads automatically.  The actual number of threads
    used is limited by the number of available processors and the CPLEX
    license.
  }
  \item{data}{
    a logical value indicating that the object returned by
    \code{\link{coxaalen}} should contain an element \code{data} that
    gives the maximal intersections and the model matrix split into
    multiplicative and additive terms.
  }
}

\value{
  A list of the above arguments with their final values.
}

\references{
  Boruvka, A. and Cook, R. J. (2015)
  \href{http://dx.doi.org/10.1111/sjos.12113}{
    A Cox-Aalen model for interval-censored data}.
  \emph{Scandinavian Journal of Statistics} \bold{42}, 414--426.

  Armijo, L. (1966)
  \href{http://dx.doi.org/10.2140/pjm.1966.16.1}{
    Minimization of functions having Lipschitz continuous first partial
    derivatives.}
  \emph{Pacific Journal of Mathematics} \bold{16}, 1--3.
}

\seealso{\code{\link{coxaalen}}}

\examples{
if (is.loaded("coxaalen", "coxinterval"))
coxaalen(Surv(left, right, type = "interval2") ~ prop(treat),
         data = cosmesis, control = coxaalen.control(iter.max = 2,
         trace = TRUE))
}

\keyword{survival}
