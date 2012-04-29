\name{boot.lqm}
\alias{boot.lqm}

\title{
Bootstrapping Linear Quantile Models
}
\description{
This function is used to obtain a bootstrap sample of a fitted LQM.
}
\usage{
boot.lqm(object, R = 50, seed = round(runif(1, 1, 10000)), startQR = TRUE)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} "lqm".
 }
  \item{R}{
number of bootstrap replications.
  }
  \item{seed}{
optional random generator seed.
  }
  \item{startQR}{
logical flag. If \code{TRUE} (default) the estimated parameters in \code{object} are used as starting values in \code{\link{lqm.fit.gs}} applied to bootstrap samples. Otherwise starting values are based on \code{\link{lm}}.}
}

\value{
an object of class "boot.lqm".
}

\author{
Marco Geraci
}

\seealso{
\code{\link{lqm}}, \code{\link{summary.lqm}}, \code{\link{lqmControl}}.
}

\keyword{lqm}
\keyword{bootstrap}
\keyword{standard errors}
