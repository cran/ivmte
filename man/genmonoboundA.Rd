% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/monobound.R
\name{genmonoboundA}
\alias{genmonoboundA}
\title{Generating monotonicity and boundedness constraints}
\usage{
genmonoboundA(
  pm0,
  pm1,
  support,
  grid_index,
  uvec,
  splinesobj,
  monov,
  uname,
  m0,
  m1,
  sset,
  gstar0,
  gstar1,
  m0.lb,
  m0.ub,
  m1.lb,
  m1.ub,
  mte.lb,
  mte.ub,
  m0.dec,
  m0.inc,
  m1.dec,
  m1.inc,
  mte.dec,
  mte.inc,
  solution.m0.min = NULL,
  solution.m1.min = NULL,
  solution.m0.max = NULL,
  solution.m1.max = NULL,
  audit.tol,
  direct
)
}
\arguments{
\item{pm0}{A list of the monomials in the MTR for d = 0.}

\item{pm1}{A list of the monomials in the MTR for d = 1.}

\item{support}{a matrix for the support of all variables that enter
into the MTRs.}

\item{grid_index}{a vector, the row numbers of \code{support} used
to generate the grid preceding the audit.}

\item{uvec}{a vector, the points in the interval [0, 1] that the
unobservable takes on.}

\item{splinesobj}{a list of lists. Each of the inner lists contains
details on the splines declared in the MTRs.}

\item{monov}{name of variable for which the monotonicity conditions
applies to.}

\item{uname}{name declared by user to represent the unobservable
term in the MTRs.}

\item{m0}{one-sided formula for marginal treatment response
function for the control group. The formula may differ from
what the user originally input in \code{\link{ivmte}}, as the
spline components should have been removed. This formula is
simply a linear combination of all covariates that enter into
the original \code{m0} declared by the user in
\code{\link{ivmte}}.}

\item{m1}{one-sided formula for marginal treatment response
function for the treated group. The formula may differ from
what the user originally input in \code{\link{ivmte}}, as the
spline components should have been removed. This formula is
simply a linear combination of all covariates that enter into
the original \code{m1} declared by the user in
\code{\link{ivmte}}.}

\item{sset}{a list containing the point estimates and gamma
components associated with each element in the S-set.}

\item{gstar0}{set of expectations for each terms of the MTR for the
control group.}

\item{gstar1}{set of expectations for each terms of the MTR for the
control group.}

\item{m0.lb}{scalar, lower bound on MTR for control group.}

\item{m0.ub}{scalar, upper bound on MTR for control group.}

\item{m1.lb}{scalar, lower bound on MTR for treated group.}

\item{m1.ub}{scalar, upper bound on MTR for treated group.}

\item{mte.lb}{scalar, lower bound on MTE.}

\item{mte.ub}{scalar, upper bound on MTE.}

\item{m0.dec}{boolean, indicating whether the MTR for the control
group is monotone decreasing.}

\item{m0.inc}{boolean, indicating whether the MTR for the control
group is monotone increasing.}

\item{m1.dec}{boolean, indicating whether the MTR for the treated
group is monotone decreasing.}

\item{m1.inc}{boolean, indicating whether the MTR for the treated
group is monotone increasing.}

\item{mte.dec}{boolean, indicating whether the MTE is monotone
decreasing.}

\item{mte.inc}{boolean, indicating whether the MTE is monotone
increasing.}

\item{solution.m0.min}{vector, the coefficients for the MTR for
\code{D = 0} corresponding to the lower bound of the target
parameter. If passed, this will initiate checks of shape
constraints.}

\item{solution.m1.min}{vector, the coefficients for the MTR for
\code{D = 1} corresponding to the lower bound of the target
parameter. If passed, this will initiate checks of shape
constraints.}

\item{solution.m0.max}{vector, the coefficients for the MTR for
\code{D = 0} corresponding to the upper bound of the target
parameter. If passed, this will initiate checks of shape
constraints.}

\item{solution.m1.max}{vector, the coefficients for the MTR for
\code{D = 1} corresponding to the upper bound of the target
parameter. If passed, this will initiate checks of shape
constraints.}

\item{audit.tol}{feasibility tolerance when performing the
audit. By default to set to be equal \code{1e-06}. This
parameter should only be changed if the feasibility tolerance
of the solver is changed, or if numerical issues result in
discrepancies between the solver's feasibility check and the
audit.}

\item{direct}{boolean, set to \code{TRUE} if the direct MTR
regression is used.}
}
\value{
a list containing a unified constraint matrix, unified
    vector of inequalities, and unified RHS vector for the
    boundedness and monotonicity constraints of an LP/QCQP problem.
}
\description{
This is a wrapper function generating the matrices and vectors
associated with the monotonicity and boundedness constraints
declared by the user. Since this function generates all the
components required for the shape constraints, it is also the
function that performs the audit. That is, MTR coefficients are
passed, then this function will verify whether they satisfy the
shape constraints.
}
