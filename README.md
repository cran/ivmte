ivmte: An R Package for Marginal Treatment Effect Methods
================
Joshua Shea and Alexander Torgovitsky

  - [Introduction](#introduction)
  - [Scope of this Vignette](#scope-of-this-vignette)
  - [Installation and Requirements](#installation-and-requirements)
  - [Usage Demonstration](#usage-demonstration)
      - [Data and Background](#data-and-background)
      - [Syntax and Output Overview](#syntax-and-output-overview)
      - [Specifying the MTR Functions](#specifying-the-mtr-functions)
      - [Specifying the Target
        Parameter](#specifying-the-target-parameter)
      - [Specifying the IV–Like
        Estimands](#specifying-the-ivlike-estimands)
      - [Specifying the Propensity
        Score](#specifying-the-propensity-score)
      - [Point Identified Models](#point-identified-models)
      - [Confidence Intervals](#confidence-intervals)
      - [Specification Tests](#specification-tests)
  - [Help, Feature Requests and Bug
    Reports](#help-feature-requests-and-bug-reports)
  - [References](#references)

## Introduction

Heckman and Vytlacil (2005) introduced the marginal treatment effect
(MTE) to provide a choice-theoretic interpretation for the widely used
instrumental variables model of Imbens and Angrist (1994). The MTE can
be used to formally extrapolate from the compliers to estimate treatment
effects for other subpopulations.

The **ivmte** package provides a flexible set of methods for conducting
this extrapolation. The package uses the moment-based framework
developed by Mogstad, Santos, and Torgovitsky (2018), which accommodates
both point identification and partial identification (bounds), both
parametric and nonparametric models, and allows the user to maintain
additional shape constraints.

## Scope of this Vignette

This vignette is intended as a guide to using **ivmte** for users
already familiar with MTE methods. The key defintions and concepts from
that literature are used without further explanation. We have written a
paper (Shea and Torgovitsky 2019) that discusses both the MTE
methodology and the usage of **ivmte**, which should be helpful for
users unfamiliar with MTE methods. The survey article by Mogstad and
Torgovitsky (2018) provides additional theoretical background on MTE
methods, including the moment-based implementation used in this module.

## Installation and Requirements

**ivmte** can be installed from CRAN via

``` r
install.packages("ivmte")
```

If you have the **devtools** package, you can install the latest version
of the module directly from our GitHub repository via

``` r
devtools::install_github("jkcshea/ivmte")
```

Two additional packages are also required for **ivmte**:

1.  **splines2** for specifying models with polynomial splines. This
    package is available on CRAN.
2.  A package for solving linear programs. There are three options here:
    1.  Gurobi and the Gurobi R package **gurobi**, which can be
        obtained from [Gurobi
        Optimization](http://www.gurobi.com/index). This option requires
        a Gurobi software license, which Gurobi Optimization offers at
        no cost to academic researchers.
    2.  CPLEX and the package **cplexAPI**, which is available on CRAN.
        CPLEX can be obtained from
        [IBM](https://www.ibm.com/analytics/cplex-optimizer). This
        option requires a software license, which IBM offers at no cost
        to academic researchers.
    3.  The **lpSolveAPI** package, which is free and open-source, and
        available from CRAN. Note that **lpSolveAPI** is a wrapper for
        [lp\_solve](http://lpsolve.sourceforge.net/5.5/), which is no
        longer actively developed.

**ivmte** tries to automatically choose a solver from those available,
with preference being given to Gurobi and CPLEX. We have provided the
option to use **lpSolveAPI** because it appears to be the only interface
for solving linear programs that can be installed solely through
`install.packages`. However, we *strongly* recommend using Gurobi or
CPLEX, since these are actively developed, much more stable, and
typically an order of magnitude faster than **lpSolveAPI**. A very clear
installation guide for Gurobi can be found
[here](https://cran.r-project.org/package=prioritizr/vignettes/gurobi_installation.html)

## Usage Demonstration

### Data and Background

We will use a subsample of 209133 women from the data used in Angrist
and Evans (1998). The data is included with **ivmte** and has the
following structure.

``` r
library(ivmte)
knitr::kable(head(AE, n = 10))
```

| worked | hours | morekids | samesex | yob |
| -----: | ----: | -------: | ------: | --: |
|      0 |     0 |        0 |       0 |  52 |
|      0 |     0 |        0 |       0 |  45 |
|      1 |    35 |        0 |       1 |  49 |
|      0 |     0 |        0 |       0 |  50 |
|      0 |     0 |        0 |       0 |  50 |
|      0 |     0 |        1 |       1 |  51 |
|      1 |    40 |        0 |       0 |  51 |
|      1 |    40 |        0 |       1 |  46 |
|      0 |     0 |        0 |       1 |  45 |
|      1 |    30 |        1 |       0 |  44 |

We will use these variables as follows:

  - `worked` is a binary outcome variable that indicates whether the
    woman was working.
  - `hours` is a multivalued outcome variable that measures the number
    of hours the woman worked per week.
  - `morekids` is a binary treatment variable that indicates whether a
    woman had more than two children or exactly two children.
      - Note that Angrist and Evans (1998) removed women with fewer than
        two children when constructing their sample.
      - Also note that the treatment variable must be binary.
  - `samesex` is an indicator for whether the woman’s first two children
    had the same sex (male-male or female-female). This is used as an
    instrument for `morekids`.
  - `yob` is the woman’s year of birth (i.e. her age), which will be
    used as a conditioning variable.

A goal of Angrist and Evans (1998) was to estimate the causal effect of
`morekids` on either `worked` or `hours`. A linear regression of
`worked` on `morekids` yields:

``` r
lm(data = AE, worked ~ morekids)
#> 
#> Call:
#> lm(formula = worked ~ morekids, data = AE)
#> 
#> Coefficients:
#> (Intercept)     morekids  
#>      0.5822      -0.1423
```

The regression shows that women with three or more children are about
14–15% less likely to be working than women with exactly two children.
Is this because children have a negative impact on a woman’s labor
supply? Or is it because women with a weaker attachment to the labor
market choose to have more children?

Angrist and Evans (1998) address this question by using `samesex` as an
instrument for `morekids`. We expect that `samesex` is randomly
assigned, so that among women with two or more children, those with weak
labor market attachment are just as likely as those with strong labor
market attachment to have to have two boys or two girls for their first
two children. A regression shows that women whose first two children had
the same sex are also more likely to go on to have a third child:

``` r
lm(data = AE, morekids ~ samesex)
#> 
#> Call:
#> lm(formula = morekids ~ samesex, data = AE)
#> 
#> Coefficients:
#> (Intercept)      samesex  
#>     0.30214      0.05887
```

Thus, `samesex` constitutes a potential instrumental variable for
`morekids`. We can run a simple instrumental variable regression using
the `ivreg` command from the **AER** package.

``` r
library("AER")
ivreg(data = AE, worked ~ morekids | samesex )
#> 
#> Call:
#> ivreg(formula = worked ~ morekids | samesex, data = AE)
#> 
#> Coefficients:
#> (Intercept)     morekids  
#>     0.56315     -0.08484
```

The coefficient on `morekids` is smaller in magnitude than it was in the
linear regression of `worked` on `morekids`. This suggests that the
linear regression overstates the negative impact that children have on a
woman’s labor supply. The likely explanation is that women who have more
children were less likely to work anyway.

An important caveat to this reasoning, first discussed by Imbens and
Angrist (1994), is that it applies only to the group of *compliers* who
would have had a third child if and only if their first two children
were same sex. (This interpretation requires the so-called
*monotonicity* condition.) The first stage regression of `morekids` on
`samesex` shows that this group comprises less than 6% of the entire
population. Thus, the complier subpopulation is a small and potentially
unrepresentative group of individuals.

Is the relationship between fertility and labor supply for the compliers
the same as for other groups? The answer is important if we want to use
an instrumental variable estimator to inform policy questions. The
purpose of **ivmte** is to provide a formal framework for answering this
type of extrapolation question.

-----

For demonstrating some of the features of **ivmte**, it will also be
useful to use a simulated dataset. The following code, which is
contained in `./extdata/ivmteSimData.R`, generates some simulated data
from a simple DGP. The simulated data is also included with the package
as `./data/ivmteSimData.rda`.

``` r
set.seed(1)
n <- 5000
u <- runif(n)
z <- rbinom(n, 3, .5)
x <- as.numeric(cut(rnorm(n), 10)) # normal discretized into 10 bins
d <- as.numeric(u < z*.25 + .01*x)

v0 <- rnorm(n) + .2*u
m0 <- 0
y0 <- as.numeric(m0 + v0 + .1*x > 0)

v1 <- rnorm(n) - .2*u
m1 <- .5
y1 <- as.numeric(m1 + v1 - .3*x > 0)

y <- d*y1 + (1-d)*y0

ivmteSimData <- data.frame(y,d,z,x)
knitr::kable(head(ivmteSimData, n = 10))
```

| y | d | z | x |
| -: | -: | -: | -: |
| 0 | 0 | 0 | 3 |
| 0 | 0 | 1 | 6 |
| 0 | 0 | 1 | 3 |
| 0 | 0 | 3 | 7 |
| 0 | 1 | 1 | 7 |
| 1 | 0 | 1 | 4 |
| 1 | 0 | 2 | 4 |
| 1 | 0 | 1 | 6 |
| 1 | 0 | 1 | 6 |
| 1 | 1 | 3 | 5 |

### Syntax and Output Overview

The main command of the **ivmte** package is also called `ivmte`. It has
the following basic syntax:

``` r
library("ivmte")
results <- ivmte(data = AE,
                 target = "att",
                 m0 = ~ u + yob,
                 m1 = ~ u + yob,
                 ivlike = worked ~ morekids + samesex + morekids*samesex,
                 propensity = morekids ~ samesex)
```

Here’s what these required parameters refer to:

  - `data = AE` is the usual reference to the data to be used.
  - `target = "att"` specifies the target parameter to be the average
    treatment on the treated (ATT).
  - `m0` and `m1` are formulas specifying the MTR functions for the
    untreated and treated arms, respectively. The symbol `u` in the
    formula refers to the unobservable latent variable in the selection
    equation.
  - `ivlike` indicates the regressions to be run to create moments to
    which the model is fit.
  - `propensity` specifies a model for the propensity score.

This is what happens when the code above is run:

``` r
results <- ivmte(data = AE,
                 target = "att",
                 m0 = ~ u + yob,
                 m1 = ~ u + yob,
                 ivlike = worked ~ morekids + samesex + morekids*samesex,
                 propensity = morekids ~ samesex)
#> 
#> LP solver: Gurobi ('gurobi')
#> 
#> Obtaining propensity scores...
#> 
#> Generating target moments...
#>     Integrating terms for control group...
#>     Integrating terms for treated group...
#> 
#> Generating IV-like moments...
#>     Moment 1...
#>     Moment 2...
#>     Moment 3...
#>     Moment 4...
#>     Independent moments: 4 
#> 
#> Performing audit procedure...
#>     Generating audit grid...
#>     Generating initial constraint grid...
#> 
#>     Audit count: 1
#>     Minimum criterion: 0
#>     Obtaining bounds...
#>     Violations: 0
#>     Audit finished.
#> 
#> Bounds on the target parameter: [-0.09835591, -0.08285941]
```

The `ivmte` function indicates its progress in performing a sequence of
intermediate operations. It then runs through an audit procedure to
enforce shape constraints. The audit procedure is described in more
detail [ahead](#audit). After the audit procedure terminates, `ivmte`
returns bounds on the target parameter, which in this case is the
average treatment on the treated (`target = "att"`). These bounds are
the primary output of interest. The detailed output can be suppressed by
passing `noisy = FALSE` as an additional parameter:

``` r
results <- ivmte(data = AE,
                 target = "att",
                 m0 = ~ u + yob,
                 m1 = ~ u + yob,
                 ivlike = worked ~ morekids + samesex + morekids*samesex,
                 propensity = morekids ~ samesex,
                 noisy = FALSE)
#> 
#> Bounds on the target parameter: [-0.09835591, -0.08285941]
#> Audit terminated successfully after 1 round.
#> 
```

### Specifying the MTR Functions

#### Basics

The required `m0` and `m1` arguments use the standard R formula syntax
familiar from functions like `lm` or `glm`. However, the formulas
involve an unobservable variable, `u`, which corresponds to the latent
propensity to take treatment in the selection equation. This variable
can be included in formulas in the same way as other observable
variables in the data. For example,

``` r
args <- list(data = AE,
             ivlike =  worked ~ morekids + samesex + morekids*samesex,
             target = "att",
             m0 = ~ u + I(u^2) + yob + u*yob,
             m1 = ~ u + I(u^2) + I(u^3) + yob + u*yob,
             propensity = morekids ~ samesex,
             noisy = FALSE)
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.2987098, 0.1289499]
#> Audit terminated successfully after 1 round.
#> 
```

A restriction that we make for computational purposes is that `u` can
only appear in polynomial terms (perhaps interacted with other
variables). Thus, the following raises an error

``` r
args[["m0"]] <- ~ log(u) + yob
r <- do.call(ivmte, args)
#> Error: The following terms are not declared properly.
#>   m0: log(u) 
#> The unobserved variable 'u' must be declared as a monomial, e.g. u, I(u^3). The monomial can be interacted with other variables, e.g. x:u, x:I(u^3). Expressions where the unobservable term is not a monomial are either not permissable or will not be parsed correctly, e.g. exp(u), I((x * u)^2). Try to rewrite the expression so that 'u' is only included in monomials.
```

Names other than `u` can be used for the selection equation
unobservable, but one must remember to pass the option `uname` to
indicate the new name.

``` r
args[["m0"]] <- ~ v + I(v^2) + yob + v*yob
args[["m1"]] <- ~ v + I(v^2) + I(v^3) + yob + v*yob
args[["uname"]] <- "v"
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.2987098, 0.1289499]
#> Audit terminated successfully after 1 round.
#> 
```

There are some limitations regarding the use of factor variables. For
example, the following formula for `m1` will trigger an error.

``` r
args[["uname"]] <- ~ "u"
args[["m0"]] <- ~ u + yob
args[["m1"]] <- ~ u + factor(yob)55 + factor(yob)60
```

However, one can work around this in a natural way.

``` r
args[["m1"]] <- ~ u + (yob == 55) + (yob == 60)
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.09835591, -0.08285941]
#> Audit terminated successfully after 1 round.
#> 
```

#### Polynomial Splines

In addition to global polynomials in `u`, `ivmte` also allows for
polynomial splines in `u` using the keyword `uSplines`. For example,

``` r
args <- list(data = AE,
             ivlike =  worked ~ morekids + samesex + morekids*samesex,
             target = "att",
             m0 = ~ u + uSplines(degree = 2, knots = c(.2, .4, .6, .8)) + yob,
             m1 = ~ uSplines(degree = 3, knots = c(.1, .3, .5, .7))*yob,
             propensity = morekids ~ samesex,
             noisy = FALSE)
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.4479412, 0.2875447]
#> Audit terminated successfully after 2 rounds.
#> 
```

The `degree` refers to the polynomial degree, so that `degree = 2` is a
quadratic spline, and `degree = 3` is a cubic spline. Constant splines,
which have an important role in some of the theory developed by Mogstad,
Santos, and Torgovitsky (2018), can be implemented with `degree = 0`.
The vector `knots` indicates the piecewise regions for the spline. Note
that `0` and `1` are always automatically included as `knots`, so that
only the interior knots need to be specified.

#### Shape Restrictions

One can also require the MTR functions to satisfy the following
nonparametric shape restrictions.

  - Boundedness of the MTR functions, via the parameters `m0.lb`,
    `m0.ub`, `m1.lb`, and `m1.ub`, which all take scalar values. In
    order to produce non-trivial bounds in cases of partial
    identification, `m0.lb` and `m1.lb` are by default set to the
    smallest value of the response variable (which is inferred from
    `ivlike`) that is observed in the data, while `m0.ub` and `m1.ub`
    are set to the largest value.
  - Boundedness of the MTE function, that is, `m1 - m0`, by setting
    `mte.lb` and `mte.ub`. By default, these are set to the values
    logically implied by the values for `m0.lb`, `m0.ub`, `m1.lb`, and
    `m1.ub`.
  - Monotonicity of the MTR functions in `u` for each value of the other
    observed variables, via the parameters `m0.dec`, `m0.inc`, `m1.dec`,
    `m1.inc`, which all take boolean values.
  - Monotonicity of the MTE function in `u`, via the parameters
    `mte.dec` and `mte.inc`.

Here is an example with monotonicity restrictions:

``` r
args <- list(data = AE,
             ivlike =  worked ~ morekids + samesex + morekids*samesex,
             target = "att",
             m0 = ~ u + uSplines(degree = 2, knots = c(.2, .4, .6, .8)) + yob,
             m1 = ~ uSplines(degree = 3, knots = c(.1, .3, .5, .7))*yob,
             m1.inc = TRUE,
             m0.inc = TRUE,
             mte.dec = TRUE,
             propensity = morekids ~ samesex,
             noisy = FALSE)
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.08827673, 0.2006589]
#> Audit terminated successfully after 1 round.
#> 
```

#### The Audit Procedure

The shape restrictions in the previous section are enforced by adding
constraints to a linear program. It is difficult in general to ensure
that a function satisfies restrictions such as boundedness or
monotonicity, even if that function is known to be a polynomial. This
difficulty is addressed in **ivmte** through an auditing procedure.

The auditing procedure is based on two grids: A relatively small
*constraint grid*, and a relatively large *audit grid*. The MTR
functions (and implied MTE function) are made to satisfy the specified
shape restrictions on the constraint grid by adding constraints to the
linear programs. While this ensures that the MTR functions satisfy the
shape restrictions on the constraint grid, it does not ensure that they
satisfy the restrictions “everywhere.” Thus, after solving the programs,
we evaluate the solution MTR functions on the large audit grid, which is
relatively inexpensive computationally. If there are points in the audit
grid at which the solution MTRs violate the desired restrictions, then
we add some of these points to the constraint grid and repeat the
procedure. The procedure ends (the audit is passed), when the solution
MTRs satisfy the constraints over the entire audit grid.

There are certain parameters that can be used to tune the audit
procedure: The number of initial points in the constraint grid
(`initgrid.nu` and `initgrid.nx`), the number of points in the audit
grid (`audit.nu` and `audit.nx`), the maximum number of points that are
added from the audit grid to the constraint grid (`audit.add`) for each
constraint, and the maximum number of times the audit procedure is
repeated before giving up (`audit.max`). We have tried to choose
defaults for these parameters that should be suitable for most
applications. However, if `ivmte` takes a very long time to run, one
might want to try adjusting these parameters. Also, if `audit.max` is
hit, which should be unlikely given the default settings, one should
either adjust the settings or examine the `audit.grid$violations` field
of the returned results to see the extent to which the shape
restrictions are not satisfied.

### Specifying the Target Parameter

The target parameter is the object the researcher wants to know.
**ivmte** has built-in support for conventional target parameters and a
class of generalized LATEs. It also has a system that allows users to
construct their own target parameters by defining the associated
weights.

#### Conventional Target Parameters

The target parameter can be set to the average treatment effect (ATE),
average treatment on the treated (ATT), or the average treatment on the
untreated by setting `target` to `ate`, `att`, or `atu`, respectively.
For example:

``` r
args <- list(data = AE,
             ivlike =  worked ~ morekids + samesex + morekids*samesex,
             target = "att",
             m0 = ~ u + I(u^2) + yob + u*yob,
             m1 = ~ u + I(u^2) + I(u^3) + yob + u*yob,
             propensity = morekids ~ samesex,
             noisy = FALSE)
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.2987098, 0.1289499]
#> Audit terminated successfully after 1 round.
#> 
args[["target"]] <- "ate"
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.3800976, 0.2003928]
#> Audit terminated successfully after 1 round.
#> 
```

#### LATEs and Generalized LATEs

LATEs can be set as target parameters by passing `target = late` and
including named lists called `late.from` and `late.to`. The named lists
should contain variable name and value pairs, where the variable names
must also appear in the propensity score formula. Typically, one would
choose instruments for these variables, although `ivmte` will let you
include control variables as well. We demonstrate this using the
simulated data.

``` r
args <- list(data = ivmteSimData,
             ivlike =  y ~ d + z + d*z,
             target = "late",
             late.from = c(z = 1),
             late.to = c(z = 3),
             m0 = ~ u + I(u^2) + I(u^3) + x,
             m1 = ~ u + I(u^2) + I(u^3) + x,
             propensity = d ~ z + x,
             noisy = FALSE)
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.6931532, -0.4397993]
#> Audit terminated successfully after 2 rounds.
#> 
```

We can condition on covariates in the definition of the LATE by adding
the named list `late.X` as follows.

``` r
args[["late.X"]] = c(x = 2)
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.8419396, -0.2913049]
#> Audit terminated successfully after 2 rounds.
#> 
args[["late.X"]] = c(x = 8)
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.7721625, -0.3209851]
#> Audit terminated successfully after 2 rounds.
#> 
```

**ivmte** also provides support for generalized LATEs, i.e. LATEs where
the intervals of `u` that are being considered may or may not correspond
to points in the support of the propensity score. These objects are
useful for a number of extrapolation purposes, see e.g. Mogstad, Santos,
and Torgovitsky (2018). They are set with `target = "genlate"` and the
additional scalars `genlate.lb` and `genlate.ub`. For example,

``` r
args <- list(data = ivmteSimData,
             ivlike =  y ~ d + z + d*z,
             target = "genlate",
             genlate.lb = .2,
             genlate.ub = .4,
             m0 = ~ u + I(u^2) + I(u^3) + x,
             m1 = ~ u + I(u^2) + I(u^3) + x,
             propensity = d ~ z + x,
             noisy = FALSE)
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.7597831, -0.2984248]
#> Audit terminated successfully after 2 rounds.
#> 
args[["genlate.ub"]] <- .41
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.7551905, -0.3083199]
#> Audit terminated successfully after 2 rounds.
#> 
args[["genlate.ub"]] <- .42
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.7504255, -0.3182317]
#> Audit terminated successfully after 2 rounds.
#> 
```

Generalized LATEs can also be made conditional-on-covariates by passing
`late.X`:

``` r
args[["late.X"]] <- c(x = 2)
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.867551, -0.2750135]
#> Audit terminated successfully after 2 rounds.
#> 
```

#### Custom Target Parameters

Since there are potentially a large and diverse array of target
parameters that a researcher might be interested in across various
applications, the **ivmte** package also allows the specification of
custom target parameters. This is done by specifying the two weighting
functions for the target parameters.

To facilitate computation, these functions must be expressible as
constant splines in `u`. For the untreated weights, the knots of the
spline are specified with `target.knots0` and the weights to place in
between each knot is given by `target.weight0`. Since the weighting
functions might depend on both the instrument and other covariates, both
`target.knots0` and `target.weight0` should be lists consisting of
functions or scalars. (A scalar is interpreted as a constant function.)
Specifying the treated weights works the same way but with
`target.knots1` and `target.weight1`. Note that the `target` option is
ignored when a set of knots and weights are passed.

Here is an example of using custom weights to replicate the conditional
LATE from above.

``` r
pmodel <- r$propensity$model # returned from the previous run of ivmte
xeval = 2 # x = xeval is the group that is conditioned on
px <- mean(ivmteSimData$x == xeval) # probability that x = xeval
z.from = 1
z.to = 3
weight1  <- function(x) {
    if (x != xeval) {
        return(0)
    } else {
        xz.from <- data.frame(x = xeval, z = z.from)
        xz.to <- data.frame(x = xeval, z = z.to)
        p.from <- predict(pmodel, newdata = xz.from, type = "response")
        p.to <- predict(pmodel, newdata = xz.to, type = "response")
        return(1 / ((p.to - p.from) * px))
    }
}
weight0 <- function(x) {
    return(-weight1(x))
}
## Define knots (same for treated and control)
knot1 <- function(x) {
    xz.from <- data.frame(x = x, z = z.from)
    p.from <- predict(pmodel, newdata = xz.from, type = "response")
    return(p.from)
}
knot2 <- function(x) {
    xz.to <- data.frame(x = x, z = z.to)
    p.to <- predict.glm(pmodel, newdata = xz.to, type = "response")
    return(p.to)
}

args <- list(data = ivmteSimData,
             ivlike =  y ~ d + z + d*z,
             target.knots0 = c(knot1, knot2),
             target.knots1 = c(knot1, knot2),
             target.weight0 = c(0, weight0, 0),
             target.weight1 = c(0, weight1, 0),
             m0 = ~ u + I(u^2) + I(u^3) + x,
             m1 = ~ u + I(u^2) + I(u^3) + x,
             propensity = d ~ z + x,
             noisy = FALSE)
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.8419396, -0.2913049]
#> Audit terminated successfully after 2 rounds.
#> 
```

The knot specification here is the same for both the treated and
untreated weights. It specifies two knots that depend on whether `x
== 2`, so that there are four knots total when accounting for 0 and 1,
which are always included automatically. These four knots divide the
interval between 0 and 1 into three regions. The first region is from 0
to the value of the propensity score when evaluated at `x` and `z.from`.
The second region is between this point and the propensity score
evaluated at `x` and `z.to`. The third and final region is from this
point up to 1.

Since the knot specification creates three regions, three weight
functions must be passed. Here, the weights in the first and third
regions are set to 0 regardless of the value of `x` by just passing a
scalar `0`. The weights in the second region are only non-zero when `x
== 2`, in which case they are equal to the inverse of the distance
between the second and third knot points. Thus, the weighting scheme
applies equal weight to compliers with `x == 2`, and zero weight to all
others. As expected, the resulting bounds match those that we computed
above using `target == late`.

### Specifying the IV–Like Estimands

The IV–like estimands refer to the moments in the data that are used to
fit the model. In **ivmte** these are restricted to be only moments that
can be expressed as the coefficient in either a linear regression (LR)
or two stage least squares (TSLS) regression of the outcome variable on
the other variables in the dataset.

There are three aspects that can be changed, which we discuss in order
subsequently. First, one can specify one or multiple LR and TSLS
formulas from which moments are used. By default, all moments from each
formula are used. Second, one can specify that only certain moments from
a formula are used. Third, one can specify a subset of the data on which
the formula is run, which provides an easy way to nonparametrically
condition on observables.

#### List of Linear Model Formulas

The required `ivlike` option expects a list of R formulas. All of our
examples up to now have had a single LR. However, one can include
multiple LRs, as well as TSLS regressions using the `|` syntax familiar
from the **AER** package. Each formula must have the same variable on
the left-hand side, which is how the outcome variable gets inferred. For
example,

``` r
args <- list(data = ivmteSimData,
             ivlike =  c(y ~ (z == 1) + (z == 2) + (z == 3) + x,
                         y ~ d + x,
                         y ~ d | z),
             target = "ate",
             m0 = ~ uSplines(degree = 1, knots = c(.25, .5, .75)) + x,
             m1 = ~ uSplines(degree = 1, knots = c(.25, .5, .75)) + x,
             propensity = d ~ z + x,
             noisy = FALSE)
r <- do.call(ivmte, args)
#> Warning: The following IV-like specifications do not include the treatment
#> variable: 1. This may result in fewer independent moment conditions than
#> expected.
#> 
#> Bounds on the target parameter: [-0.6427017, -0.3727193]
#> Audit terminated successfully after 1 round.
#> 
```

There are 10 moments being fit in this specification. Five of these
moments correspond to the constant term, the coefficients on the three
dummies `(z == 1)`, `(z == 2)` and `(z == 3)`, and the coefficient on
`x` from the first linear regression. There are three more coefficients
from the second regression of `y` on `d`, `x`, and a constant. Finally,
there are two moments from the TSLS (simple IV in this case) regression
of `y` on `d` and a constant, using `z` and a constant as instruments.

#### Using Only Subcomponents

The default behavior of `ivmte` is to use all of the moments from each
formula included in `ivlike`. There are reasons one might not want to do
this, such as if certain moments are estimated poorly (i.e. with
substantial statistical error). For more information, see the
discussions in Mogstad, Santos, and Torgovitsky (2018) and Mogstad and
Torgovitsky (2018).

One can tell `ivmte` to only include certain moments and not others by
passing the `components` option. This option expects a list of the same
length as the list passed to `ivlike`. Each component of the list is
itself a list that contains the variable names for the coefficients to
be included from that formula. The list should be declared using the `l`
function, which is a generalization of the `list` function. The `l`
function allows the user to list variables and expressions without
having to enclose them by quotation marks. For example, the following
includes the coefficients on the intercept and `x` from the first
formula, the coefficient on `d` from the second formula, and all
coefficients in the third formula.

``` r
args[["components"]] <- l(c(intercept, x), c(d), )
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.6865291, -0.2573525]
#> Audit terminated successfully after 1 round.
#> 
```

Note that `intercept` is a reserved word that is used to specify the
coefficient on the constant term. If the function `list` is used to pass
the `components` option, an error will follow.

``` r
args[["components"]] <- list(c(intercept, x), c(d), )
#> Error in eval(expr, envir, enclos): object 'intercept' not found
args[["components"]] <- list(c("intercept", "x"), c("d"), "")
r <- do.call(ivmte, args)
#> Error in terms.formula(fi, ...): invalid model formula in ExtractVars
```

#### Subsetting

The formulas can be run conditional on certain subgroups by adding the
`subset` option. This option expects a list of logical statements with
the same length as `ivlike` declared using the `l` function. One can use
the entire data by leaving the statement blank, or inserting a tautology
such as `1 == 1`. For example, the following would run the first
regression only on observations with `x` less than or equal to 9, the
second regression on the entire sample, and the third (TSLS) formula
only on those observations that have `z` equal to 1 or 3.

``` r
args <- list(data = ivmteSimData,
             ivlike =  c(y ~ z + x,
                         y ~ d + x,
                         y ~ d | z),
             subset = l(x <= 9, 1 == 1, z %in% c(1,3)),
             target = "ate",
             m0 = ~ uSplines(degree = 3, knots = c(.25, .5, .75)) + x,
             m1 = ~ uSplines(degree = 3, knots = c(.25, .5, .75)) + x,
             propensity = d ~ z + x,
             noisy = FALSE)
r <- do.call(ivmte, args)
#> Warning: The following IV-like specifications do not include the treatment
#> variable: 1. This may result in fewer independent moment conditions than
#> expected.
#> 
#> Bounds on the target parameter: [-0.6697228, -0.3331582]
#> Audit terminated successfully after 2 rounds.
#> 
```

### Specifying the Propensity Score

The procedure implemented by `ivmte` requires first estimating the
propensity score, that is, the probability that the binary treatment
variable is 1, conditional on the instrument and other covariates. This
must be specified with the `propensity` option, which expects a formula.
The treatment variable is inferred to be the variable on the left-hand
side of the `propensity`. The default is to estimate the formula as a
logit model via the `glm` command, but this can be changed to `"linear"`
or `"probit"` with the `link` option.

``` r
results <- ivmte(data = AE,
                 target = "att",
                 m0 = ~ u + yob,
                 m1 = ~ u + yob,
                 ivlike = worked ~ morekids + samesex + morekids*samesex,
                 propensity = morekids ~ samesex + yob,
                 link = "probit",
                 noisy = FALSE)
#> 
#> Bounds on the target parameter: [-0.100781, -0.0825274]
#> Audit terminated successfully after 1 round.
#> 
```

### Point Identified Models

The default of the `ivmte` command is to estimate bounds on the target
parameter. If one believes the model is in fact point identified, then
`point == TRUE` can be included as an option. Point identification will
typically occur when `ivlike`, `components` and `subset` are such that
they have more (non-redundant) components than there are free parameters
in `m0` and `m1`. Bounds in this case will typically shrink to a point.

``` r
args <- list(data = ivmteSimData,
             ivlike =  y ~ d + factor(z),
             target = "ate",
             m0 = ~ u,
             m1 = ~ u,
             propensity = d ~ factor(z),
             noisy = FALSE)
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.5349027, -0.5349027]
#> Audit terminated successfully after 1 round.
#> 
args[["point"]] <- TRUE
r <- do.call(ivmte, args)
#> 
#> Point estimate of the target parameter: -0.5389508
```

Notice that even though `point = FALSE` and `point = TRUE` provide a
single number for the target parameter in this scenario, that number is
slightly different. The reason is that, for computational reasons,
`point = TRUE` uses the usual quadratic loss, while `point = FALSE` uses
absolute deviations loss. One should use `point = TRUE` if the model is
point identified, since it computes confidence intervals and
specification tests (discussed ahead) in a well-understood way. The
default when `point = TRUE` is to use the optimal two-step GMM
weighting, however the identity weighting can be computed by passing
`point.eyeweight = TRUE`.

### Confidence Intervals

The `ivmte` command does provide functionality for constructing
confidence regions, although this is turned off by default, since it can
be computationally expensive. To turn it on, set `bootstraps` to be a
positive integer.

``` r
results <- ivmte(data = AE,
                 target = "att",
                 m0 = ~ u + yob,
                 m1 = ~ u + yob,
                 ivlike = worked ~ morekids + samesex + morekids*samesex,
                 propensity = morekids ~ samesex,
                 bootstraps = 50,
                 noisy = FALSE)
#> 
#> Bounds on the target parameter: [-0.09835591, -0.08285941]
#> Audit terminated successfully after 1 round. 
#> 
#> Bootstrapped confidence intervals (backward):
#>     90%: [-0.2013432, -0.01624443]
#>     95%: [-0.2185989, -0.003729633]
#>     99%: [-0.2300391, 0.0006812143]
#> p-value: 0.04
#> Number of bootstraps: 50
```

Other options relevent to confidence region construction are
`bootstraps.m`, which indicates the number of observations to draw from
the sample on each bootstrap replication, and `bootstraps.replace` to
indicate whether these observations should be drawn with or without
replacement. The default is to set `bootstraps.m` to the sample size of
the observed data with `bootstraps.replace = TRUE`. This corresponds to
the usual nonparametric bootstrap. Choosing `bootstraps.m` to be smaller
than the sample size and setting `bootstraps.replace` to be `TRUE` or
`FALSE` corresponds to the “m out of n” bootstrap or subsampling,
respectively. Regardless of these settings, two types of confidence
regions are constructed following the terminology of Andrews and Han
(2009); see Shea and Torgovitsky (2019) for more detail, references, and
important theoretical caveats to the procedures.

Confidence regions can also be constructed when `point == TRUE` in a
similar way.

``` r
args <- list(data = AE,
             target = "att",
             m0 = ~ u,
             m1 = ~ u,
             ivlike = worked ~ morekids + samesex + morekids*samesex,
             propensity = morekids ~ samesex,
             point = TRUE,
             bootstraps = 50,
             noisy = FALSE)
r <- do.call(ivmte, args)
#> 
#> Point estimate of the target parameter: -0.09160436
#> 
#> Bootstrapped confidence intervals (nonparametric):
#>     90%: [-0.1869028, -0.02510774]
#>     95%: [-0.2077314, -0.02507778]
#>     99%: [-0.227466, 0.0001054838]
#> p-value: 0.08
#> Number of bootstraps: 50
```

The p-value reported in both cases is for the null hypothesis that the
target parameter is equal to 0, which we infer here by finding the
largest confidence region that does not contain 0. By default, `ivmte`
returns 99%, 95%, and 90% confidence intervals. This can be changed with
the `levels` option.

### Specification Tests

The moment-based framework implemented by `ivmte` is amenable to
specification tests. These tests are based on whether the minimum value
of the criterion function is statistically different from 0. In the
point identified case, this is the well known GMM overidentification
test (Hansen 1982). Here, we implement it via bootstrapping as discussed
by Hall and Horowitz (1996), because the moments depend on the estimated
propensity score which is estimated in a first stage. In the partially
identified case, we implement the misspecification test discussed by
Bugni, Canay, and Shi (2015). The specification tests are automatically
conducted when `bootstraps` is positive and the minimum criterion value
in the sample problem is larger than 0. However, it can be turned off by
setting `specification.test = FALSE`.

``` r
args <- list(data = ivmteSimData,
             ivlike =  y ~ d + factor(z),
             target = "ate",
             m0 = ~ u,
             m1 = ~ u,
             m0.dec = TRUE,
             m1.dec = TRUE,
             propensity = d ~ factor(z),
             bootstraps = 50,
             noisy = FALSE)
r <- do.call(ivmte, args)
#> 
#> Bounds on the target parameter: [-0.6189484, -0.6189484]
#> Audit terminated successfully after 1 round. 
#> 
#> Bootstrapped confidence intervals (backward):
#>     90%: [-0.6411833, -0.6011875]
#>     95%: [-0.64164, -0.5966297]
#>     99%: [-0.6465837, -0.5944238]
#> p-value: 0
#> Bootstrapped specification test p-value: 0.36
#> Number of bootstraps: 50

args[["ivlike"]] <- y ~ d + factor(z) + d*factor(z) # many more moments
args[["point"]] <- TRUE
r <- do.call(ivmte, args)
#> Warning: If argument 'point' is set to TRUE, then shape restrictions on m0
#> and m1 are ignored, and the audit procedure is not implemented.
#> 
#> Point estimate of the target parameter: -0.5559325
#> 
#> Bootstrapped confidence intervals (nonparametric):
#>     90%: [-0.6066197, -0.5112431]
#>     95%: [-0.6101133, -0.5086833]
#>     99%: [-0.6147004, -0.4743232]
#> p-value: 0
#> Bootstrapped J-test p-value: 0.44
#> Number of bootstraps: 50
```

## Help, Feature Requests and Bug Reports

Please post an issue on the [GitHub
repository](https://github.com/jkcshea/ivmte/issues).

## References

<div id="refs" class="references">

<div id="ref-andrewshan2009ej">

Andrews, Donald W. K., and Sukjin Han. 2009. “Invalidity of the
Bootstrap and the M Out of N Bootstrap for Confidence Interval Endpoints
Defined by Moment Inequalities.” *Econometrics Journal* 12: S172–S199.
<http://dx.doi.org/10.1111/j.1368-423X.2008.00265.x>.

</div>

<div id="ref-angristevans1998taer">

Angrist, Joshua D., and William N. Evans. 1998. “Children and Their
Parents’ Labor Supply: Evidence from Exogenous Variation in Family
Size.” *The American Economic Review* 88 (3): 450–77.
<http://www.jstor.org/stable/116844>.

</div>

<div id="ref-bugnicanayshi2015joe">

Bugni, Federico A., Ivan A. Canay, and Xiaoxia Shi. 2015. “Specification
Tests for Partially Identified Models Defined by Moment Inequalities.”
*Journal of Econometrics* 185 (1): 259–82.
<http://www.sciencedirect.com/science/article/pii/S0304407614002577>.

</div>

<div id="ref-hallhorowitz1996e">

Hall, Peter, and Joel L. Horowitz. 1996. “Bootstrap Critical Values for
Tests Based on Generalized-Method-of-Moments Estimators.” *Econometrica*
64 (4): 891–916. <http://www.jstor.org/stable/2171849>.

</div>

<div id="ref-hansen1982e">

Hansen, Lars Peter. 1982. “Large Sample Properties of Generalized Method
of Moments Estimators.” *Econometrica* 50 (4): 1029–54.
<http://www.jstor.org/stable/1912775>.

</div>

<div id="ref-heckmanvytlacil2005e">

Heckman, James J., and Edward Vytlacil. 2005. “Structural Equations,
Treatment Effects, and Econometric Policy Evaluation.” *Econometrica* 73
(3): 669–738. <http://dx.doi.org/10.1111/j.1468-0262.2005.00594.x>.

</div>

<div id="ref-imbensangrist1994e">

Imbens, Guido W., and Joshua D. Angrist. 1994. “Identification and
Estimation of Local Average Treatment Effects.” *Econometrica* 62 (2):
467–75. <http://www.jstor.org/stable/2951620>.

</div>

<div id="ref-mogstadsantostorgovitsky2018e">

Mogstad, Magne, Andres Santos, and Alexander Torgovitsky. 2018. “Using
Instrumental Variables for Inference About Policy Relevant Treatment
Parameters.” *Econometrica* 86 (5): 1589–1619.
<https://dx.doi.org/10.3982/ecta15463>.

</div>

<div id="ref-mogstadtorgovitsky2018aroe">

Mogstad, Magne, and Alexander Torgovitsky. 2018. “Identification and
Extrapolation of Causal Effects with Instrumental Variables.” *Annual
Review of Economics* 10 (1).
<https://dx.doi.org/10.1146/annurev-economics-101617-041813>.

</div>

<div id="ref-sheatorgovitsky2019wp">

Shea, Joshua Ka Chun, and Alexander Torgovitsky. 2019. “Ivmte: An R
Package for Implementing Marginal Treatment Effect Methods.” *Working
Paper*.

</div>

</div>
