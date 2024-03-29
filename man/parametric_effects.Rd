% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parametric-effects.R
\name{parametric_effects}
\alias{parametric_effects}
\alias{parametric_effects.gam}
\title{Estimated values for parametric model terms}
\usage{
parametric_effects(object, ...)

\method{parametric_effects}{gam}(
  object,
  terms = NULL,
  data = NULL,
  unconditional = FALSE,
  unnest = TRUE,
  ci_level = 0.95,
  envir = environment(formula(object)),
  transform = FALSE,
  ...
)
}
\arguments{
\item{object}{a fitted model object.}

\item{...}{arguments passed to other methods.}

\item{terms}{character; which model parametric terms should be drawn? The
Default of \code{NULL} will plot all parametric terms that can be drawn.}

\item{data}{a optional data frame that may or may not be used? FIXME!}

\item{unconditional}{logical; should confidence intervals include the
uncertainty due to smoothness selection? If \code{TRUE}, the corrected Bayesian
covariance matrix will be used.}

\item{unnest}{logical; unnest the parametric effect objects?}

\item{ci_level}{numeric; the coverage required for the confidence interval.
Currently ignored.}

\item{envir}{an environment to look up the data within.}

\item{transform}{logical; if \code{TRUE}, the parametric effect will be plotted on
its transformed scale which will result in the effect being a straight
line. If FALSE, the effect will be plotted against the raw data (i.e. for
\code{log10(x)}, or \code{poly(z)}, the x-axis of the plot will be \code{x} or \code{z}
respectively.)}
}
\description{
Estimated values for parametric model terms
}
