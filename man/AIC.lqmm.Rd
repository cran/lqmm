\name{AIC.lqmm}
\alias{AIC.lqmm}

\title{
Akaike Information Criterion for LQMM}

\description{
This function calculates the AIC of a fitted \code{lqmm} model.
}

\usage{
\method{AIC}{lqmm}(object, ..., k = 2)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} "lqmm".
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
\code{\link{lqmm}}
\code{\link{logLik.lqmm}}
}

\keyword{AIC}
