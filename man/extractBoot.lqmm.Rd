\name{extractBoot.lqmm}
\alias{extractBoot.lqmm}

\title{
Extract Fixed and Random Bootstrapped Parameters
}
\description{
This function extracts the fixed and random components of bootstrapped estimates of an \code{lqmm} object.
}
\usage{
extractBoot.lqmm(object, which = "fixed")
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{boot.lqmm}.
}
  \item{which}{
character indicating whether \code{"fixed"} or \code{"random"} parameters.
}
}
\details{
The \code{"random"} parameters refer to the "raw" parameters of the variance-covariance matrix of the random effects as returned by \code{\link{lqmm.fit.gs}} and \code{\link{lqmm.fit.df}}.
}

\value{
a matrix of bootstrapped estimates.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{boot.lqmm}}, \code{\link{lqmm.fit.gs}}, \code{\link{lqmm.fit.df}}
}


\examples{

## Orthodont data
data(Orthodont)

# Random intercept model
fit <- lqmm(distance ~ age, random = ~ 1, group = Subject,
	tau = 0.5, data = Orthodont)
fit.boot <- boot.lqmm(fit)

# extract fixed effects
B <- extractBoot.lqmm(fit.boot, "fixed")

# covariance matrix estimated fixed parameters
cov(B)

}


\keyword{bootstrap}
