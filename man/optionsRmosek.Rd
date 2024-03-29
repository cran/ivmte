% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lp.R
\name{optionsRmosek}
\alias{optionsRmosek}
\title{Function to parse options for Gurobi}
\usage{
optionsRmosek(options, debug)
}
\arguments{
\item{options}{list. Each set of options should be passed as a
list, with the name of each entry being the name of the class
of options. For example, options for double parameters should
be contained in the entry\code{dparam = list(BASIS_TOL_X = 1e-06)}.}

\item{debug}{boolean, indicates whether or not the function should
provide output when obtaining bounds. The output provided is
the same as what Mosek would send to the console.}
}
\value{
list, the set of options declared by the user, including
    some additional default values.
}
\description{
This function constructs a list of options to be parsed when
\code{solver} is set to \code{Rmosek}. This function really
implements the default feasibility tolerances.
}
