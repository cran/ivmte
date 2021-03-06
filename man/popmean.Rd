% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testfunctions_covariates.R
\name{popmean}
\alias{popmean}
\title{Calulating population mean}
\usage{
popmean(formula, distribution, density = "f")
}
\arguments{
\item{formula}{formula, each term of which will have its mean
calculated.}

\item{distribution}{data.table, characterizing the distribution of
the variables entering into \code{formula}.}

\item{density}{string, name of the variable \code{data}
characterizing the density.}
}
\value{
vector, the means for each term in \code{formula}.
}
\description{
Given a distribution, this function calculates the population mean
for each term in a formula.
}
