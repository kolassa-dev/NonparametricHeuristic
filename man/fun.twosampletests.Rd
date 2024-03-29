% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twosampletests.R
\name{fun.twosampletests}
\alias{fun.twosampletests}
\title{Calculate a vector of two-sample tests of the null hypothesis of equal center}
\usage{
fun.twosampletests(x, y)
}
\arguments{
\item{x}{the data vector for the first sample.}

\item{y}{the data vector for the second sample.}
}
\value{
A vector of two-sided p-values
}
\description{
Calculates t-test, exact Wilcoxon, approximate Wilcoxon, normal scores, Savage scores, and Mood tests.
}
\details{
The exact Wilcoxon is performed only for small samples.
}
\examples{
fun.twosampletests(rnorm(10),rnorm(15))

}
