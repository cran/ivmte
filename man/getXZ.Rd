% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/callcheck.R
\name{getXZ}
\alias{getXZ}
\title{Auxiliary function: extract X and Z covariates from a formula}
\usage{
getXZ(fm, inst = FALSE, terms = FALSE, components = FALSE)
}
\arguments{
\item{fm}{the formula.}

\item{inst}{boolean expression, set to TRUE if the instrument names
are to be extracted. Otherwise, the covariate names are
extracted.}

\item{terms}{boolean expression, set to TRUE if the terms in the
formula \code{fm} should be returned instead of the variable
names.}

\item{components}{boolean expression, set to FALSE by
default. Indicates that the formula being considered is
constructed from a list of components, and thus the term
'intercept' is permitted.}
}
\value{
vector of variable names.
}
\description{
Auxiliary function that takes in a two-sided formula, and extracts
the variable names of either the covariates or instruments. The
function returns an error if the formula includes a variable called
'intercept'.
}
