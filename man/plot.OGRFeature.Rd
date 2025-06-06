% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3_methods.R
\name{plot.OGRFeature}
\alias{plot.OGRFeature}
\title{Plot the geometry of an \code{OGRFeature} object}
\usage{
\method{plot}{OGRFeature}(x, xlab = "x", ylab = "y", main = "", ...)
}
\arguments{
\item{x}{An \code{OGRFeature} object.}

\item{xlab}{Title for the x axis.}

\item{ylab}{Title for the y axis.}

\item{main}{The main title (on top).}

\item{...}{Optional arguments passed to \code{wk::wk_plot()}.}
}
\value{
The input, invisibly.
}
\description{
Plot the geometry of an \code{OGRFeature} object
}
