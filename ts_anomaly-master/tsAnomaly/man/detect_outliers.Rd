% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tsAnomaly.R
\name{detect_outliers}
\alias{detect_outliers}
\title{Performs outlier detection}
\usage{
detect_outliers(x, STL = FALSE, sigma = 0.01, ccF = 2.5)
}
\arguments{
\item{x}{Integer. A vector of integers.}

\item{STL}{Logical. Should the series be seasonally detrended before isolation
forest is run.}

\item{sigma}{Numeric. The significance level for outlier detection.}

\item{ccF}{Numeric. The significance level for outlier detection. See Details.}
}
\value{
A list containing both outlier indicies from \code{\link{iso_forest}}
and \code{\link{count_dist_outlier}}.
}
\description{
Performs outlier detection
}

