\name{addnoise}
\alias{addnoise}

\title{
Jittering
}
\description{
This function adds a uniformly generated noise to a vector.
}
\usage{
addnoise(x, centered = TRUE, B = 0.999)
}

\arguments{
  \item{x}{
a numeric vector.
}
  \item{centered}{
logical flag. If \code{TRUE} (default) noise is drawn from U[0,B] and if FALSE from U[-B/2,B/2].
}
  \item{B}{
right boundary for uniform random noise U[0,B].
}
}

\value{
a dithered vector.
}

\author{
Marco Geraci
}

\seealso{
\code{\link{lqm.counts}}
}


\keyword{jittering}

