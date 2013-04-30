\name{predict.lqm}
\alias{predict.lqm}

\title{
Predictions from an \code{lqm} Object
}

\description{
This function computes predictions based on fitted linear quantile model.
}

\usage{
\method{predict}{lqm}(object, newdata, interval = FALSE,
	level = 0.95, na.action = na.pass, ...)
}

\arguments{
  \item{object}{
an \code{lqm} object.
}
  \item{newdata}{
an optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.
}
  \item{interval}{
logical flag. If \code{TRUE}, bootstrap percentile intervals for predictions are provided.
}
  \item{level}{
confidence level.
}
  \item{na.action}{
function determining what should be done with missing values in \code{newdata}. The default is to predict \code{NA}.}

  \item{\dots}{
further arguments passed to \code{boot.lqm}.
}

}

\value{
a vector or a matrix or an array of predictions.
}


\author{
Marco Geraci
}


\seealso{
\code{\link{residuals.lqm}}, \code{\link{lqm}}, \code{\link{coef.lqmm}}, \code{\link{boot.lqm}}
}

\keyword{prediction}
\keyword{lqm}
