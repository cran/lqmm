\name{lqm.counts}
\alias{lqm.counts}
\title{
Quantile Regression for Counts}
\description{
This function is used to fit a quantile regression model when the response is a count variable.
}
\usage{
lqm.counts(formula, data, weights = NULL, offset = NULL, contrasts = NULL,
	alpha = 0.5, M = 50, zeta = 1e-05, B = 0.999, cn = NULL)
}

\arguments{
  \item{formula}{
an object of class \code{\link{formula}}: a symbolic description of the model to be fitted.
}
  \item{data}{
an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lqm is called.
}
  \item{weights}{
an optional vector of weights to be used in the fitting process.
}

  \item{offset}{
an optional offset to be included in the model frame.  
}

  \item{contrasts}{
an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
}
  \item{alpha}{
quantile to be estimated.
}
  \item{M}{
number of dithered samples.
}
  \item{zeta}{
small constant (see References).
}

  \item{B}{
right boundary for uniform random noise U[0,B] to be added to the response variable (see References).
}

 \item{cn}{
small constant to be passed to \code{\link{F.lqm}} (see Reference).
}

}

\details{A linear quantile regression model if fitted to the log--transformed response. Additional tranformation functions will be implemented. The notation used here follows closely that of Machado and Santos Silva (2005).}

\value{
an object of class "list" containing the following components

\item{estimate_log_scale}{regression quantile on the log--scale}
\item{std.err}{standard error}
\item{fitted}{predicted quantile}
\item{effective}{number of dithered samples used for standard error estimation that gave an invertible D matrix (Machado and Santos Silva, 2005).}
}

\references{
Machado JAF and Santos Silva JMC (2005). Quantiles for counts. Journal of the American Statistical Association, 100(472), 1226--1237. 
}
\author{
Marco Geraci
}



\examples{

n <- 100
x <- runif(n)
data <- data.frame(x = x, y = rpois(n, 2*x))
lqm.counts(y ~ x, data = data, M = 20)


}
\keyword{quantiles for counts}
