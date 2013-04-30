\name{coef.lqm}
\alias{coef.lqm}

\title{
Extract LQM Coefficients
}
\description{
\code{coef} extracts model coefficients from \code{lqm} objects.
}
\usage{
\method{coef}{lqm}(object, ...)
}

\arguments{
  \item{object}{
a fitted object of \code{\link{class}} "lqm".
}
  \item{\dots}{
not used.
}
}

\value{
a vector for single quantiles or a matrix for multiple quantiles. 
}

\author{
Marco Geraci
}

\seealso{
\code{\link{lqm}}
\code{\link{summary.lqm}}
}

\keyword{lqm}
\keyword{coefficients}

