\name{coxdual}
\alias{coxdual}
\alias{na.coxdual}
\alias{na.coxdual.default}
\alias{na.coxdual.data.frame}

\title{Cox model for a Markov illness-death process under dual censoring}

\description{
  Fit a Cox model to a progressive Markov illness-death process observed
  under right-censored survival times and interval- or right-censored
  progression times.
}

\usage{
coxdual(formula, data = parent.frame(), subset, init = NULL,
        formula.coxph = NULL, init.coxph = FALSE, control, ...)
}

\arguments{
  \item{formula}{
    an expression of the form \code{response ~ terms}, where
    \code{response} is a survival object returned by the
    \code{\link[survival]{Surv}} function and \code{terms} contain
    clustering and state-transition variables identified by the
    \code{\link[survival]{cluster}} and \code{\link{trans}} functions,
    respectively.
  }
  \item{data}{
    an optional data frame in which to interpret the variables named in
    the arguments \code{formula} and \code{formula.coxph}.
  }
  \item{subset}{
    expression specifying which rows of \code{data} should be used
    in the fit.  All observations are included by default.
  }
  \item{init}{
    a named list of the vector \code{coef}, specifying the initial
    coefficient values, and matrix or data frame \code{basehaz}, on
    which to base the initial cumulative baseline transition
    intensities.  The \code{basehaz} element should contain columns
    representing the cumulative transition intensity, time and
    transition type. If the columns do not appear in this order, they
    should be indicated by the column names "hazard", "time" and
    "trans", respectively.  The sorted values used to represent the
    transition types should respectively denote the initial to
    intermediate, initial to terminal, and intermediate to terminal
    state transitions.  The \code{NULL} value for \code{init} or its
    components enables default values.  For \code{coef} the default is
    zero.  For \code{basehaz} the default corresponds to linear
    functions of time with an upper bound of one.  Under any
    alternatives arising from this or the \code{init.coxph} arguments,
    \code{basehaz} is interpreted as a step function of time and the
    initial value is its piecewise linear approximation.
  }
  \item{formula.coxph}{
    a list of formula objects specifying models to fit with
    \code{\link[survival]{coxph}} using singly right-censored
    observations in \code{data}.  The shorthand \code{.} indicates
    the same terms given in \code{formula}, with the function
    \code{\link[survival]{strata}} in place by
    \code{\link{trans}}.  Under dual-right--censoring, the first
    component of \code{formula.coxph} defaults to same model in
    \code{formula}, which is fit to observations singly-right--censored
    at the earlier censoring time.
  }
  \item{init.coxph}{
    a logical value indicating that \code{init} should be overrided by
    estimates based on the \code{\link[survival]{coxph}} fit to the
    first model specified in \code{formula.coxph}.
  }
  \item{control}{
    a named list of parameters controlling the model fit, returned by
    the function \code{\link{coxdual.control}}.  This defaults to
    \code{coxdual.control()}.
  }
  \item{\dots}{
    additional arguments to be passed to \code{\link{coxdual.control}}.
  }
}

\details{
  A valid \code{formula} argument can be expressed
  as \preformatted{Surv(<start>, <stop>, <status>)
    ~ cluster(<id>) + trans(<from>, <to>) + <covariate terms>}
  where (\code{<start>}, \code{<stop>}] is largest \emph{known} time
  interval over which individual \code{<id>} is at risk for a transition
  between the states \code{<from>} and \code{<to>}.  The variable
  \code{<status>} indicates whether or not a transition is observed to
  occur at \code{<stop>}.

  Under dual censoring (Boruvka and Cook, 2014), both the originating
  state and the left endpoint of an at-risk interval may be unknown.
  This case is handled with \code{<start> = NA}, \code{<from> = NA},
  \code{<to>} equal to the index of the terminal state and any
  transition-type--specific covariates taking on the values assumed when
  \code{<from>} is equal to the intermediate state index.  Under
  discrete observation of non-terminal events, the right-endpoint of
  some at-risk intervals may be unknown.  For these \code{<start>} is
  the initial observation time (zero, unless left-truncated),
  \code{<stop> = NA} and \code{<from>} is equal to the initial state
  index.  Missing values are retained by the \code{NA} action
  \code{\link{na.coxdual}}.  The default \code{NA} action is used to
  handle any missing values passed to \code{\link[survival]{coxph}}.

  Dual censoring typically arises in two scenarios: (1) dual
  right-censoring, where intermediate events are right-censored before
  terminal events, and (2) interval-censored intermediate events.
  For examples of these refer to \code{\link{dualrc}} and
  \code{\link{dualic}}, respectively.

  A consequence of dual censoring is that any discrete maximum
  likelihood estimator has ambiguous support at any failure times
  associated with these \code{NA} values.  To resolve this, the
  cumulative baseline transition intensities are restricted to piecewise
  linear functions on a sieve with size controlled by arguments passed
  to \code{\link{coxdual.control}}.  This approach requires that both types
  of transitions to the terminal state are, at least for some subjects,
  observed exactly.
}

\value{
  An object of the classes \code{"coxinterval"} and \code{"coxdual"},
  which is a list with the following components.
  \item{call}{
    the matched call to \code{coxdual}.
  }
  \item{censor}{
    a string indicating the dual censoring type.  The value "right"
    corresponds to strictly dual right-censored data.  All other cases
    return "interval".
  }
  \item{n}{
    size of the sample used in the model fit.
  }
  \item{m}{
    number of at-risk intervals used in the model fit.
  }
  \item{p}{
    number of regression coefficients.
  }
  \item{coef}{
    a named \code{p} vector of regression coefficients.
  }
  \item{var}{
    a named \code{p} by \code{p} covariance matrix of the regression
    coefficients.
  }
  \item{basehaz}{
    a data frame giving the cumulative baseline transition intensities
    evaluated over the sieve.
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
  \item{gradnorm}{
    the maximum norm of the score function at the final parameter value.
  }
  \item{maxnorm}{
    the maximum norm of the difference between the penultimate and final
    parameter values.
  }
  \item{cputime}{
    the processing time used for parameter and variance estimation.
  }
  \item{fit.coxph}{
    an unnamed list whose components give the
    \code{\link[survival]{coxph}} fit to any models specified by the
    \code{formula.coxph} argument.
  }
  \item{na.action}{
    the \code{"na.action"} attribute of the model frame.  Here this
    corresponds to the result from the custom \code{NA} action
    \code{\link{na.coxdual}}.
  }
  \item{censor.rate}{
    a named vector of censoring rates.
  }
  \item{control}{
    a list of arguments passed to \code{\link{coxdual.control}}.
  }
  \item{data}{
    a list containing the data relevant to the model fit, if the
    \code{\link{coxdual.control}} argument \code{data} is true.
  }
}

\references{
  Boruvka, A. and Cook, R. J. (2014)
  Sieve estimation in a Markov illness-death process under dual censoring.
}

\seealso{
  \code{\link[survival]{cluster}}, \code{\link{dualic}},
  \code{\link{dualrc}}, \code{\link[survival]{Surv}},
  \code{\link{trans}}
}

\examples{
# Fit Cox model to dual-right--censored data
fit <- coxdual(Surv(start, stop, status) ~ cluster(id) + trans(from, to)
               + I(z * (to == 1)) + I(z * (from \%in\% 0 & to == 2))
               + I(z * (from \%in\% c(NA, 1) & to == 2)), data = dualrc,
               sieve.rate = 2/5)
fit
par(mfrow = c(1, 3))
by(fit$basehaz, fit$basehaz$trans, function(x) plot(x[, 2:1],
   type = "l", main = paste(x[1, 3]), xlim = c(0, 2), ylim = c(0, 4)))

# Fit Cox model to data with interval-censored progression times
fit <- coxdual(Surv(start, stop, status) ~ cluster(id) + trans(from, to)
               + I(z * (to == 1)) + I(z * (from \%in\% 0 & to == 2))
               + I(z * (from \%in\% c(NA, 1) & to == 2)), data = dualic)
fit
par(mfrow=c(1, 3))
by(fit$basehaz, fit$basehaz$trans, function(x) plot(x[, 2:1],
   type = "l", main = paste(x[1, 3]), xlim = c(0, 2), ylim = c(0, 4)))
}

\keyword{survival}
