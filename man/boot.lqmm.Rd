\name{boot.lqmm}
\alias{boot.lqmm}

\title{
Bootstrapping Linear Quantile Mixed Models
}
\description{
This function is used to obtain a bootstrap sample of a fitted LQMM.
}
\usage{
boot.lqmm(object, R = 50, seed = round(runif(1, 1, 10000)), startQR = FALSE)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} "lqmm".
 }
  \item{R}{
number of bootstrap replications.
  }
  \item{seed}{
optional random generator seed.
  }
  \item{startQR}{
logical flag. If \code{TRUE} the estimated parameters in \code{object} are used as starting values in \code{lqmm.fit} applied to bootstrap samples. Otherwise starting values are based on \code{\link{lm}}.}

}
  \value{
an object of class "boot.lqmm".
}

\author{
Marco Geraci
}

\seealso{
\code{\link{lqmm}}, \code{\link{summary.lqmm}}, \code{\link{lqmmControl}}.
}

\keyword{lqmm}
\keyword{bootstrap}
\keyword{standard errors}
