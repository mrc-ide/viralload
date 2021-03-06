% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viralload.R
\name{prepare_observed}
\alias{prepare_observed}
\title{Prepare observed data}
\usage{
prepare_observed(observed)
}
\arguments{
\item{observed}{A list of \code{data.frame}s, each with columns \code{vl}
and \code{count}}
}
\value{
An object of class \code{observed} to pass through to
\link{log_likelihood}. \strong{DO NOT ALTER THIS OBJECT}
}
\description{
Prepare observed data
}
