\name{raneff.lqmm}
\alias{raneff.lqmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random Effects
}
\description{
This function computes random effects for a linear quantile mixed model.
}
\usage{
raneff.lqmm(object)
}

\arguments{
  \item{object}{
an object of \code{\link{class}} \code{lqmm}.
}
}

\details{
The prediction of the random effects is done via estimated best linear prediction. See references for details.
}

\value{
a data frame or a list of data frames of predicted random effects.
}

\references{
Geraci M and Bottai M (1 June 2011). Linear Quantile Mixed Models. Unpublished manuscript.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{lqmm}}, \code{\link{coef.lqmm}}
}

\keyword{random effects}
\keyword{coefficients}
\keyword{lqmm}
