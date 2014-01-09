\name{boot.lqm}
\alias{boot.lqm}

\title{
Bootstrapping Linear Quantile Models
}
\description{
This function is used to obtain a bootstrap sample of a fitted LQM.
}
\usage{
boot.lqm(object, R = 50, seed = round(runif(1, 1, 10000)), startQR = FALSE)
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
logical flag. If \code{TRUE} the estimated parameters in \code{object} are used as starting values in \code{\link{lqm.fit.gs}} applied to bootstrap samples. Otherwise starting values are based on \code{\link{lm}}.}
}

\value{
An object of class \code{boot.lqm} is a data frame with \code{R} rows and \code{npars} columns containing the bootstrap estimates of \code{theta}. If \code{object} contains results for multiple quantiles, \code{boot.lqm} returns an array of dimension \code{c(R,npars,nt)}, where \code{nt} is the length of \code{tau}.

The following attributes are available:

\item{tau}{index of the quantile(s).}
\item{estimated}{the estimated parameter as given by \code{object}.}
\item{R}{number of bootstrap replications.}
\item{seed}{the random number generator seed used to produce the bootstrap sample.}
\item{npars}{total numer of parameters.}
\item{rdf}{the number of residual degrees of freedom.}
\item{indices}{the bootstrap sample of independent data units.}

}

\author{
Marco Geraci
}

\seealso{
\code{\link{lqm}}, \code{\link{summary.lqm}}, \code{\link{lqmControl}}.
}

\examples{

set.seed(123)
n <- 500
test <- data.frame(x = runif(n,0,1))
test$y <- 30 + test$x + rnorm(n)
fit.lqm <- lqm(y ~ x, data = test, tau = 0.5)
boot.lqm(fit.lqm)

}

\keyword{bootstrap}
\keyword{standard errors}
