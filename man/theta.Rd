% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution-helpers.R
\name{theta}
\alias{theta}
\alias{theta.gam}
\title{General extractor for additional parameters in mgcv models}
\usage{
theta(object, ...)

\method{theta}{gam}(object, transform = TRUE, ...)
}
\arguments{
\item{object}{a fitted model}

\item{...}{arguments passed to other methods.}

\item{transform}{logical; transform to the natural scale of the parameter}
}
\value{
Returns a numeric vector of additional parameters
}
\description{
General extractor for additional parameters in mgcv models
}
\examples{
load_mgcv()
df <- data_sim("eg1", dist = "poisson", seed = 42, scale = 1 / 5)
m <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3),
  data = df, method = "REML",
  family = nb()
)
p <- theta(m)
}
