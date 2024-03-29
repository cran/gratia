% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/concurvity.R
\name{draw.pairwise_concurvity}
\alias{draw.pairwise_concurvity}
\alias{draw.overall_concurvity}
\title{Plot concurvity measures}
\usage{
\method{draw}{pairwise_concurvity}(
  object,
  title = "Smooth-wise concurvity",
  subtitle = NULL,
  caption = NULL,
  x_lab = "Term",
  y_lab = "With",
  fill_lab = "Concurvity",
  continuous_colour = NULL,
  ...
)

\method{draw}{overall_concurvity}(
  object,
  title = "Overall concurvity",
  subtitle = NULL,
  caption = NULL,
  y_lab = "Concurvity",
  x_lab = NULL,
  bar_col = "steelblue",
  bar_fill = "steelblue",
  ...
)
}
\arguments{
\item{object}{An object inheriting from class \code{"concurvity"}, usually the
result of a call to \code{\link[=model_concurvity]{model_concurvity()}} or its abbreviated form
\code{\link[=concrvity]{concrvity()}}.}

\item{title}{character; the plot title.}

\item{subtitle}{character; the plot subtitle.}

\item{caption}{character; the plot caption}

\item{x_lab}{character; the label for the x axis.}

\item{y_lab}{character; the label for the y axis.}

\item{fill_lab}{character; the label to use for the fill guide.}

\item{continuous_colour}{function; continuous colour (fill) scale to use.}

\item{...}{arguments passed to other methods.}

\item{bar_col}{colour specification for the bar colour.}

\item{bar_fill}{colour specification for the bar fill}
}
\description{
Plot concurvity measures
}
