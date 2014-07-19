\name{trans}
\alias{trans}

\title{Identify transition type in model terms}

\description{
  A special function for Cox models fit with \code{\link{coxic}} that
  identifies formula terms specifying the state-transition type.
}

\usage{trans(from, to)}

\arguments{
  \item{from}{
    a variable representing the originating state.
  }
  \item{to}{
    a variable representing the subsequent state.
  }
}

\value{
  A combination of the \code{from} and \code{to} arguments by column
  with two attributes:
  \item{\code{"states"}}{
    a vector giving the unique non-missing values in the \code{from} and
    \code{to} arguments ordered so that the initial state appears first,
    the intermediate state second, and the terminal state last.
  }
  \item{\code{"types"}}{
    a vector of transition type labels in terms of the values in the
    \code{from} and \code{to} arguments ordered so that the
    intermediate transition appears first, the terminal transition
    directly from the initial state second, and the terminal transition
    from the intermediate state last.
  }
}

\seealso{
  \code{\link{coxic}}
}

\examples{with(dualrc, trans(from, to))}

\keyword{survival}
