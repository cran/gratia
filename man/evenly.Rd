% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evenly.R
\name{evenly}
\alias{evenly}
\alias{seq_min_max}
\title{Create a sequence of evenly-spaced values}
\usage{
evenly(x, n = 100, by = NULL, lower = NULL, upper = NULL)

seq_min_max(x, n, by = NULL, lower = NULL, upper = NULL)
}
\arguments{
\item{x}{numeric; vector over which evenly-spaced values are returned}

\item{n}{numeric; the number of evenly-spaced values to return. A default of
\code{100} is used for convenience as that what is typically used when
evaluating a smooth.}

\item{by}{numeric; the increment of the sequence. If specified, argument \code{n}
is ignored and the sequence returned will be from \code{min(x)} to \code{max(x)} in
increments of \code{by}.}

\item{lower}{numeric; the lower bound of the interval.}

\item{upper}{numeric; the upper bound of the interval.}
}
\value{
A numeric vector of length \code{n}.
}
\description{
For a continuous vector \code{x}, \code{evenly} and \code{seq_min_max()}
create a sequence of \code{n} evenly-spaced values over the range \code{lower}
-- \code{upper}. By default, \code{lower} is defined as \code{min(x)} and \code{upper} as
\code{max(x)}, excluding \code{NA}s. For a factor \code{x}, the function returns
\code{levels(x)}.
}
\examples{
\dontshow{
set.seed(1)
}
x <- rnorm(10)
n <- 10L

# 10 values evenly over the range of `x`
evenly(x, n = n)

# evenly spaced values, incrementing by 0.2
evenly(x, by = 0.2)

# evenly spaced values, incrementing by 0.2, starting at -2
evenly(x, by = 0.2, lower = -2)
}
\seealso{
See \code{\link[base:seq]{base::seq()}} for details of the behaviour of \code{evenly()} when
using \code{by}.
}
