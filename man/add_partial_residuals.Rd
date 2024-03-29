% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add-methods.R
\name{add_partial_residuals}
\alias{add_partial_residuals}
\alias{add_partial_residuals.gam}
\title{Add partial residuals}
\usage{
add_partial_residuals(data, model, ...)

\method{add_partial_residuals}{gam}(data, model, select = NULL, partial_match = FALSE, ...)
}
\arguments{
\item{data}{a data frame containing values for the variables used to fit the
model. Passed to \code{\link[stats:residuals]{stats::residuals()}} as \code{newdata}.}

\item{model}{a fitted model for which a \code{\link[stats:residuals]{stats::residuals()}} method is
available. S3 method dispatch is performed on the \code{model} argument.}

\item{...}{arguments passed to other methods.}

\item{select}{character, logical, or numeric; which smooths to plot. If
\code{NULL}, the default, then all model smooths are drawn. Numeric \code{select}
indexes the smooths in the order they are specified in the formula and
stored in \code{object}. Character \code{select} matches the labels for smooths
as shown for example in the output from \code{summary(object)}. Logical
\code{select} operates as per numeric \code{select} in the order that smooths are
stored.}

\item{partial_match}{logical; should smooths be selected by partial matches
with \code{select}? If \code{TRUE}, \code{select} can only be a single string to match
against.}
}
\description{
Add partial residuals
}
\examples{
load_mgcv()
\dontshow{
op <- options(pillar.sigfig = 3, cli.unicode = FALSE)
}
df <- data_sim("eg1", seed = 1)
df <- df[, c("y", "x0", "x1", "x2", "x3")]
m <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = df, method = "REML")

## add partial residuals
add_partial_residuals(df, m)

## add partial residuals for selected smooths
add_partial_residuals(df, m, select = "s(x0)")
\dontshow{
options(op)
}
}
