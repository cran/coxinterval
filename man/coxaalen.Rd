\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{coxaalen}
\alias{coxaalen}

\title{Cox-Aalen model for interval-censored survival data}

\description{
  Fit a Cox-Aalen model to interval-censored survival data with fixed
  covariates.
}

\usage{
coxaalen(formula, data = parent.frame(), subset, init = NULL,
         formula.timereg = NULL, init.timereg = FALSE,
         close.cplex = TRUE, control, ...)
}

\arguments{
  \item{formula}{
    an expression of the form \code{response ~ terms}, where
    \code{response} is a survival object returned by the
    \code{\link[survival]{Surv}} function and \code{terms} contains at
    least one multiplicative term identified by the
    \code{\link[timereg]{prop}} function from the \CRANpkg{timereg}
    package.
  }
  \item{data}{
    an optional data frame in which to interpret the variables named in
    the arguments \code{formula} and \code{formula.timereg}.
  }
  \item{subset}{
    expression specifying which rows of \code{data} should be used
    in the fit.  All observations are included by default.
  }
  \item{init}{
    a list with elements named \code{coef} and \code{basehaz}.  The
    \code{coef} element should be a scalar or vector specifying the
    initial values of the multiplicative regression coefficient.  If
    \code{init = NULL} or \code{coef = NULL}, this coefficient will be
    initialized to zero.  The element \code{basehaz} should be a matrix
    or data frame whose columns represent time and the corresponding
    value for the cumulative baseline hazard function and any remaining
    additive cumulative coefficients (also referred to as
    "regression functions").  The cumulative regression functions are
    then initialized by linearly interpolating or extrapolating these
    values. If \code{init = NULL} or \code{basehaz = NULL}, the
    cumulative baseline hazard function is initialized to a linear
    function of time and any cumulative coefficients are started at
    zero.
  }
  \item{formula.timereg}{
    a list of formula objects specifying models to fit with
    the \CRANpkg{timereg} package's \code{\link[timereg]{cox.aalen}}
    function using right-censored observations in \code{data}.  Here the
    shorthand \code{"."} refers to the same \code{terms} given in
    \code{formula}.
  }
  \item{init.timereg}{
    a logical value indicating that \code{init} should be overrided by
    estimates based on the \code{\link[timereg]{cox.aalen}} fit to the
    first model in \code{formula.timereg}.
  }
  \item{close.cplex}{
    a logical value indicating that CPLEX data structures should be
    freed when \code{coxaalen} exits.
  }
  \item{control}{
    a named list of parameters controlling the model fit, as returned by
    the function \code{\link{coxaalen.control}}.  This defaults to
    \code{coxaalen.control()}.
  }
  \item{\dots}{
    additional arguments to be passed to
    \code{\link{coxaalen.control}}.
  }
}

\details{
  A valid response in the \code{formula} argument can be expressed
  as \preformatted{Surv(<left>, <right>, type = "interval2")}
  where (\code{<left>}, \code{<right>}] is the censoring interval for
  the survival time.  Following the \CRANpkg{survival} package's
  \code{type = "interval2"} censoring for the
  \code{\link[survival]{Surv}} function, we use the convention that any
  right-censoring times are provided in the variable \code{<left>}
  and \code{<right>} is set to the \code{NA} value.

  Terms in the \code{formula} have either time-varying additive effects
  on the survival hazard as in Aalen's additive regression model, or
  fixed multiplicative effects as in the Cox model. Multiplicative terms
  are distinguished by applying \CRANpkg{timereg}'s
  \code{\link[timereg]{prop}} function to each corresponding variable.

  \code{coxaalen} depends on libraries that are loaded only if
  \code{coxinterval} is installed from source on a system with
  \href{http://www-03.ibm.com/software/products/en/ibmilogcpleoptistud/}{
    IBM ILOG CPLEX Optimization Studio}. Refer to the package's INSTALL
  file for detailed instructions.
}

\value{
  An object of the class \code{"coxinterval"} and \code{"coxaalen"},
  which is a list with the following components.
  \item{call}{
    the matched call to \code{coxdual}.
  }
  \item{n}{
    size of the sample used in the model fit.
  }
  \item{p}{
    number of (multiplicative) regression coefficients.
  }
  \item{coef}{
    a named \code{p} vector of regression coefficients.
  }
  \item{var}{
    a named \code{p} by \code{p} covariance matrix of the regression
    coefficients.
  }
  \item{basehaz}{
    a data frame giving the cumulative regression functions evaluated at
    time points given by the maximal intersections of the censoring
    intervals.
  }
  \item{init}{
    list of initial values used in the model fit.
  }
  \item{loglik}{
    a vector giving the initial and final log-likelihood values.
  }
  \item{iter}{
    number of iterations needed to achieve the model fit.
  }
  \item{maxnorm}{
    the maximum norm of the difference between the penultimate and final
    parameter values.
  }
  \item{gradnorm}{
    the inner product between the final parameter value and the score
    function.
  }
  \item{cputime}{
    the processing time for parameter and variance estimation.
  }
  \item{fit.timereg}{
    an unnamed list whose components give the
    \code{\link[timereg]{cox.aalen}} fit to any models specified by the
    \code{formula.timereg} argument.
  }
  \item{na.action}{
    the \code{"na.action"} attribute of the model frame.
  }
  \item{censor.rate}{
    a vector giving the rates of exact, left-censored, interval-censored
    and right-censored observations used in the model fit.
  }
  \item{control}{
    a named list of arguments passed to \code{\link{coxdual.control}}.
  }
}

\references{
  Boruvka, A. and Cook, R. J. (2015)
  \href{http://dx.doi.org/10.1111/sjos.12113}{
    A Cox-Aalen model for interval-censored data}.
  \emph{Scandinavian Journal of Statistics} \bold{42}, 414--426.

  Martinussen, T. and Scheike, T. H. (2006)
  Dynamic Regression Models for Survival Data.
  New York: Springer.

  Scheike, T. H. and Zhang, M.-J. (2002)
  \href{http://dx.doi.org/10.1111/1467-9469.00065}{
    An additive-multiplicative Cox-Aalen regression model}.
  \emph{Scandinavian Journal of Statistics} \bold{29}, 75--88.
}

\seealso{
  \code{\link[timereg]{cox.aalen}}, \code{\link[timereg]{prop}},
  \code{\link[survival]{Surv}}
}

\examples{
# Fit a Cox model to the breast cosmesis dataset
if (is.loaded("coxaalen", "coxinterval")) \{
  fit <- coxaalen(Surv(left, right, type = "interval2") ~ prop(treat),
                  data = cosmesis, init.timereg = TRUE,
                  formula.timereg = list(Surv(pmax(left, right, na.rm = TRUE),
                  !is.na(right)) ~ .))
  summary(fit)
\}
}

\keyword{survival}
