% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mst.R
\name{gmmEstimate}
\alias{gmmEstimate}
\title{GMM estimate of TE under point identification}
\usage{
gmmEstimate(
  sset,
  gstar0,
  gstar1,
  center = NULL,
  subsetList = NULL,
  n = NULL,
  redundant = NULL,
  identity = FALSE,
  nMoments,
  splines,
  noisy = TRUE
)
}
\arguments{
\item{sset}{a list of lists constructed from the function
\link{genSSet}. Each inner list should include a coefficient
corresponding to a term in an IV specification, a matrix of the
estimates of the gamma moments conditional on (X, Z) for the
control group, and a matrix of the estimates of the gamma
moments conditional on (X, Z) for the treated group. The column
means of the last two matrices is what is used to generate the
gamma moments.}

\item{gstar0}{vector, the target gamma moments for the control
group.}

\item{gstar1}{vector, the target gamma moments for the treated
group.}

\item{center}{numeric, the GMM moment equations from the original
sample. When bootstrapping, the solution to the point
identified case obtained from the original sample can be passed
through this argument to recenter the bootstrap distribution of
the J-statistic.}

\item{subsetList}{list of subset indexes, one for each IV-like
specification.}

\item{n}{number of observations in the data. This option is only
used when subsetting is involved.}

\item{redundant}{vector of integers indicating which components in
the S-set are redundant.}

\item{identity}{boolean, default set to \code{FALSE}. Set to
\code{TRUE} if GMM point estimate should use the identity
weighting matrix (i.e. one-step GMM).}

\item{nMoments}{number of linearly independent moments. This option
is used to determine the cause of underidentified cases.}

\item{splines}{boolean, set to \code{TRUE} if the MTRs involve
splines. This option is used to determine the cause of
underidentified cases.}

\item{noisy}{boolean, default set to \code{TRUE}. If \code{TRUE},
then messages are provided throughout the estimation
procedure. Set to \code{FALSE} to suppress all messages,
e.g. when performing the bootstrap.}
}
\value{
a list containing the point estimate of the treatment
    effects, and the MTR coefficient estimates. The moment
    conditions evaluated at the solution are also returned, along
    with the J-test results. However, if the option \code{center}
    is passed, then the moment conditions and J-test are centered
    (this is to perform the J-test via bootstrap).
}
\description{
If the user sets the argument \code{point = TRUE} in the function
\code{ivmte}, then it is assumed that the treatment effect
parameter is point identified. The observational equivalence
condition is then set up as a two-step GMM problem. Solving this
GMM problem recovers the coefficients on the MTR functions m0 and
m1. Combining these coefficients with the target gamma moments
allows one to estimate the target treatment effect.
}
\examples{
dtm <- ivmte:::gendistMosquito()

## Declare empty list to be updated (in the event multiple IV like
## specifications are provided
sSet <- list()

## Declare MTR formulas
formula1 = ~ 0 + u
formula0 = ~ 0 + u

## Construct object that separates out non-spline components of MTR
## formulas from the spline components. The MTR functions are
## obtained from this object by the function 'genSSet'.
splinesList = list(removeSplines(formula0), removeSplines(formula1))

## Construct MTR polynomials
polynomials0 <- polyparse(formula = formula0,
                 data = dtm,
                 uname = u,
                 as.function = FALSE)
polynomials1 <- polyparse(formula = formula0,
                 data = dtm,
                 uname = u,
                 as.function = FALSE)

## Generate propensity score model
propensityObj <- propensity(formula = d ~ z,
                            data = dtm,
                            link = "linear")

## Generate IV estimates
ivEstimates <- ivEstimate(formula = ey ~ d | z,
                          data = dtm,
                          components = l(intercept, d),
                          treat = d,
                          list = FALSE)

## Generate target gamma moments
targetGamma <- genTarget(treat = "d",
                         m0 = ~ 1 + u,
                         m1 = ~ 1 + u,
                         target = "atu",
                         data = dtm,
                         splinesobj = splinesList,
                         pmodobj = propensityObj,
                         pm0 = polynomials0,
                         pm1 = polynomials1)

## Construct S-set. which contains the coefficients and weights
## corresponding to various IV-like estimands
sSet <- genSSet(data = dtm,
                sset = sSet,
                sest = ivEstimates,
                splinesobj = splinesList,
                pmodobj = propensityObj$phat,
                pm0 = polynomials0,
                pm1 = polynomials1,
                ncomponents = 2,
                scount = 1,
                yvar = "ey",
                dvar = "d",
                means = FALSE)

## Obtain point estimates using GMM
gmmEstimate(sset = sSet$sset,
            gstar0 = targetGamma$gstar0,
            gstar1 = targetGamma$gstar1)

}
