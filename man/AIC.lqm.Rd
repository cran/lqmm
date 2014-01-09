\name{AIC.lqm}
\alias{AIC.lqm}

\title{
Akaike Information Criterion for LQM}

\description{
This function calculates the AIC of a fitted \code{lqm} model.
}

\usage{
\method{AIC}{lqm}(object, ..., k = 2)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} "lqm".
}
  \item{\dots}{
not used.
}
  \item{k}{
numeric, the penalty per parameter to be used; the default \code{k = 2} is the classical AIC.
}
}


\author{
Marco Geraci
}


\seealso{
\code{\link{lqm}}
\code{\link{logLik.lqm}}
}

\keyword{AIC}
