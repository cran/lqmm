\name{summary.lqmm}
\alias{summary.lqmm}

\title{Summary for an \code{lqmm} Object}

\description{
Summary method for class \code{lqmm}.}

\usage{
\method{summary}{lqmm}(object, method = "boot", alpha = 0.05, covariance = FALSE, ...)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{lqmm}.
	}

  \item{method}{
specifies the method used to compute standard errors. Currently, only the bootstrap method ("boot") is available.
	}
  \item{alpha}{
significance level.
	}
  \item{covariance}{
logical flag. If \code{TRUE} the bootstrapped variance-covariance matrix is provided.
	}
  \item{\dots}{
see \code{\link{boot.lqmm}} for additional arguments.
	}
}

\details{
\code{print.summary.lqmm} formats the coefficients, standard errors, etc. and additionally gives `significance stars'.}

\value{
an object of class \code{summary.lqmm}. The function \code{summary.lqmm} computes and returns a list of summary statistics of the fitted linear quantile mixed model given in \code{object}, using the components (list elements) from its argument, plus

\item{Cov}{the covariance matrix obtained from the bootstrapped estimates (if \code{covariance = TRUE}).}
\item{tTable}{a matrix with estimates, standard errors, etc.}
}

\author{
Marco Geraci}

\seealso{
\code{\link{print.summary.lqmm}}
\code{\link{lqmm}}
}

\examples{
data(Orthodont)
fitOi.lqmm <- lqmm(distance ~ age, random = ~ 1, group = Subject,
	iota = c(0.1,0.5,0.9), data = Orthodont)
summary(fitOi.lqmm)
}

\keyword{lqmm}
\keyword{bootstrap}
\keyword{standard errors}
