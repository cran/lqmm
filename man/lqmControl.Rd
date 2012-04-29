\name{lqmControl}
\alias{lqmControl}

\title{
Control parameters for lqm estimation
}
\description{
A list of parameters for controlling the fitting process.
}
\usage{
lqmControl(loop_tol_ll = 1e-6, loop_tol_theta = 1e-6, check_theta = FALSE,
	loop_step = NULL, beta = 0.5, gamma = 1, reset_step = FALSE,
	loop_max_iter = 500, verbose = FALSE)
}

\arguments{
  \item{loop_tol_ll}{
tolerance expressed as absolute change of the log-likelihood.
}
  \item{loop_tol_theta}{
tolerance expressed as absolute change of the estimates.
}
  \item{check_theta}{
logical flag. If \code{TRUE} the algorithm performs an additional check on the change in the estimates.
}

  \item{loop_step}{
step size (default standard deviation of response).
}

  \item{beta}{
decreasing step factor for line search (0,1).
}

  \item{gamma}{
nondecreasing step factor for line search (>= 1).
}

  \item{reset_step}{
logical flag. If \code{TRUE} the step size is resetted to the initial value at each iteration.
}

  \item{loop_max_iter}{
maximum number of iterations.
}
  \item{verbose}{
logical flag.
}
}

\value{
a list of control parameters.

}

\author{
Marco Geraci}

\seealso{
\code{\link{lqm}}
}

\keyword{fitting}
\keyword{control}
