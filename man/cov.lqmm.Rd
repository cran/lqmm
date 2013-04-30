\name{cov.lqmm}
\alias{cov.lqmm}

\title{
Extract Variance-Covariance Matrix
}
\description{
\code{cov.lqmm} extracts the variance-covariance matrix from a fitted \code{lqmm} object.
}
\usage{
cov.lqmm(object, ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} "lqmm".
}
  \item{\dots}{
not used.
}
}

\details{
This function returns the variance or the variance-covariance matrix of the random effects. It calls \code{covHandling} to manage the output of \code{\link{lqmm.fit.gs}} or \code{\link{lqmm.fit.df}}. A post-fitting approximation to the nearest positive (semi)definite matrix (Higham, 2002) is applied if necessary.
}

\references{
Higham N (2002). Computing the Nearest Correlation Matrix - A Problem from Finance. IMA Journal of Numerical Analysis, 22, 329-343.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{lqmm}}
\code{\link{coef.lqmm}}
}

\keyword{lqmm}
\keyword{covariance}
\keyword{coefficients}
