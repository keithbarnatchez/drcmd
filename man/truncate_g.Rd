% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{truncate_g}
\alias{truncate_g}
\title{Truncate treatment propensity scores}
\usage{
truncate_g(x, cutoff = 0.025)
}
\arguments{
\item{x}{A vector of treatment propensity scores}
}
\value{
A vector of treatment propensity scores truncated to interval \link{c, 1-c}
}
\description{
Truncate propensity scores to interval \link{c, 1-c}
}
