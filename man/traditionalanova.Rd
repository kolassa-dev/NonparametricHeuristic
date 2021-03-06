% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/traditionalanova.R
\name{traditionalanova}
\alias{traditionalanova}
\title{Construct a analysis of variance table for multiple regression, highlighting
the null hypothesis of no effect from any of the covariates.  This table
is common in statistical text books.}
\usage{
traditionalanova(object)
}
\arguments{
\item{object}{A fit as created by lm.}
}
\value{
A matrix with three rows, corresponing to model, error and total,
  and five columns, representing sum of squares, degrees of freedom, mean
  squares, the F statistic, and the p-value.
}
\description{
Construct a analysis of variance table for multiple regression, highlighting
the null hypothesis of no effect from any of the covariates.  This table
is common in statistical text books.
}
