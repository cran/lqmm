\name{predict.lqmm}
\alias{predict.lqmm}

\title{
Predictions from an \code{lqmm} Object
}

\description{
The predictions at level 0 correspond to population predictions (based only on the fixed effects estimates). The predictions at level 1 are obtained by adding the best linear predictions of the random effects to the predictions at level 0.
}

\usage{
\method{predict}{lqmm}(object, level = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an \code{lqmm} object.
}
  \item{level}{
an optional integer vector giving the level of grouping to be used in obtaining the predictions. Level zero corresponds to the population predictions.
}
  \item{\dots}{
not used.
}
}

\value{
a matrix of predictions.
}

\references{
Geraci M and Bottai M (1 June 2011). Linear Quantile Mixed Models. Unpublished manuscript.
}

\author{
Marco Geraci
}


\seealso{
\code{\link{lqmm}}, \code{\link{raneff.lqmm}}, \code{\link{coef.lqmm}}
}

\keyword{prediction}
\keyword{lqmm}
