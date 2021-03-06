% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/common.R
\name{genmultscore}
\alias{genmultscore}
\title{Calculate the general multi-group score test.}
\usage{
genmultscore(x, g, s = NULL)
}
\arguments{
\item{x}{response vector}

\item{g}{group indicator.}

\item{s}{scores to be averaged over ties.}
}
\value{
an object of class htest, with the chi-square approximate p-valu.e
}
\description{
Calculate the general multi-group score test.
}
