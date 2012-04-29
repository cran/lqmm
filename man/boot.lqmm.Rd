\name{boot.lqmm}
\alias{boot.lqmm}

\title{
Bootstrapping Linear Quantile Mixed Models
}
\description{
This function is used to obtain a bootstrap sample of a fitted LQMM.
}
\usage{
boot.lqmm(object, R = 50, seed = round(runif(1, 1, 10000)))
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
