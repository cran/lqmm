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
optional random number generator seed.
  }
  \item{startQR}{
logical flag. If \code{TRUE} the estimated parameters in \code{object} are used as starting values in \code{lqmm.fit} applied to bootstrap samples. Otherwise starting values are based on \code{\link{lm}}.}

}
  \value{
An object of class \code{boot.lqmm} is a data frame with \code{R} rows and \code{npars} columns containing the bootstrap estimates of \code{theta_x}, \code{theta_z}, and \code{scale}. If \code{object} contains results for multiple quantiles, \code{boot.lqmm} returns an array of dimension \code{c(R,npars,nt)}, where \code{nt} is the length of \code{tau}.

The elements of \code{theta_z} are labelled with \code{reStruct}. See function \code{\link{covHandling}} and the example below on how to derive the variance-covariance matrix of the random effects starting from \code{theta_z}. The following attributes are available:

\item{tau}{index of the quantile(s).}
\item{estimated}{the estimated parameter as given by \code{object}.}
\item{R}{number of bootstrap replications.}
\item{seed}{the random number generator seed used to produce the bootstrap sample.}
\item{nn}{labels of the fixed effects.}
\item{npars}{total numer of parameters.}
\item{indices}{the bootstrap sample of independent data units.}

}

\references{
Geraci M and Bottai M (2013). Linear quantile mixed models. Statistics and Computing, doi: 10.1007/s11222-013-9381-9.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{lqmm}}, \code{\link{summary.lqmm}}, \code{\link{lqmmControl}}.
}


\examples{

## Orthodont data
data(Orthodont)

# Random intercept model
fit <- lqmm(distance ~ age, random = ~ 1, group = Subject,
	tau = 0.5, data = Orthodont)
fit.boot <- boot.lqmm(fit)


}


\keyword{bootstrap}
\keyword{standard errors}
