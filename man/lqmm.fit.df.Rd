\name{lqmm.fit.df}
\alias{lqmm.fit.df}

\title{
Linear Quantile Mixed Models Fitting by Derivative-Free Optimization
}
\description{
This function controls the arguments to be passed to \code{\link{optim}} and \code{\link{optimize}} for LQMM estimation.
}

\usage{
lqmm.fit.df(theta_0, x, y, z, weights, cov_name, V, W, sigma_0,	
	iota, group, control)
}

\arguments{
  \item{theta_0}{
starting values for the linear predictor.
}
  \item{x}{
the model matrix for fixed effects.
}
  \item{y}{
the model response.
}
  \item{z}{
the model matrix for random effects.
}
  \item{weights}{
the weights used in the fitting process.
}

  \item{cov_name}{
variance--covariance matrix of the random effects. Default is \code{pdIdent}. See details.
}

  \item{V}{
nodes of the quadrature.
}
  \item{W}{
weights of the quadrature.
}
  \item{sigma_0}{
starting value for the scale parameter.
  }
  \item{iota}{
the quantile(s) to be estimated.
}
  \item{group}{
the grouping factor.
}
  \item{control}{
list of control parameters used for optimization (see \code{\link{lqmmControl}}).
}
}
\details{
In \code{\link{lqmm}}, see argument \code{fit} for generating a list of arguments to be called by this function; see argument \code{covariance} for alternative variance--covariance matrices.
}

\value{
An object of class "list" containing the following components:

\item{theta}{a vector of coefficients, including the "raw" variance--covariance parameters (see \code{\link{cov.lqmm}}).}
\item{scale}{the scale parameter.}
\item{logLik}{the log--likelihood.}
\item{opt}{number of iterations when the estimation algorithm stopped for lower (theta) and upper (scale) loop.}.
}

\author{
Marco Geraci}


\seealso{
\code{\link{lqmm}}
}

\examples{
set.seed(123)

M <- 50
n <- 10
test <- data.frame(x = runif(n*M,0,1), group = rep(1:M,each=n))
test$y <- 10*test$x + rep(rnorm(M, 0, 2), each = n) + rchisq(n*M, 3)
lqmm.ls <- lqmm(fixed = y ~ x, random = ~ 1, group = group, data = test,
	fit = FALSE)

do.call("lqmm.fit.df", lqmm.ls)

}

\keyword{fitting}
