## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE
)
require(pander)
require(data.table)
require(AER)
require(ivmte)

## ----eval = FALSE--------------------------------------------------------
#  install.packages("ivmte")

## ----eval = FALSE--------------------------------------------------------
#  devtools::install_github("jkcshea/ivmte")

## ---- drawData-ae--------------------------------------------------------
library(ivmte)
knitr::kable(head(AE, n = 10))

## ---- ols, results='markdown'--------------------------------------------
lm(data = AE, worked ~ morekids)

## ---- fs, results='markdown'---------------------------------------------
lm(data = AE, morekids ~ samesex)

## ---- ivreg, results='markdown'------------------------------------------
library("AER")
ivreg(data = AE, worked ~ morekids | samesex )

## ---- drawData-sim, results='markdown'-----------------------------------
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

## ---- syntax, eval = FALSE-----------------------------------------------
#  library("ivmte")
#  results <- ivmte(data = AE,
#                   target = "att",
#                   m0 = ~ u + yob,
#                   m1 = ~ u + yob,
#                   ivlike = worked ~ morekids + samesex + morekids*samesex,
#                   propensity = morekids ~ samesex)

## ---- syntaxrun----------------------------------------------------------
results <- ivmte(data = AE,
                 target = "att",
                 m0 = ~ u + yob,
                 m1 = ~ u + yob,
                 ivlike = worked ~ morekids + samesex + morekids*samesex,
                 propensity = morekids ~ samesex)

## ---- syntaxrun.quiet----------------------------------------------------
results <- ivmte(data = AE,
                 target = "att",
                 m0 = ~ u + yob,
                 m1 = ~ u + yob,
                 ivlike = worked ~ morekids + samesex + morekids*samesex,
                 propensity = morekids ~ samesex,
                 noisy = FALSE)

## ---- mtrbasics----------------------------------------------------------
args <- list(data = AE,
             ivlike =  worked ~ morekids + samesex + morekids*samesex,
             target = "att",
             m0 = ~ u + I(u^2) + yob + u*yob,
             m1 = ~ u + I(u^2) + I(u^3) + yob + u*yob,
             propensity = morekids ~ samesex,
             noisy = FALSE)
r <- do.call(ivmte, args)

## ---- non.polynomial.u.error, error = TRUE-------------------------------
args[["m0"]] <- ~ log(u) + yob
r <- do.call(ivmte, args)

## ---- uname--------------------------------------------------------------
args[["m0"]] <- ~ v + I(v^2) + yob + v*yob
args[["m1"]] <- ~ v + I(v^2) + I(v^3) + yob + v*yob
args[["uname"]] <- "v"
r <- do.call(ivmte, args)

## ---- factor.error, eval = FALSE-----------------------------------------
#  args[["uname"]] <- ~ "u"
#  args[["m0"]] <- ~ u + yob
#  args[["m1"]] <- ~ u + factor(yob)55 + factor(yob)60

## ---- factor.ok.prepare, echo = FALSE------------------------------------
args[["uname"]] <- ~ "u"
args[["m0"]] <- ~ u + yob

## ---- factor.ok----------------------------------------------------------
args[["m1"]] <- ~ u + (yob == 55) + (yob == 60)
r <- do.call(ivmte, args)

## ---- spline-------------------------------------------------------------
args <- list(data = AE,
             ivlike =  worked ~ morekids + samesex + morekids*samesex,
             target = "att",
             m0 = ~ u + uSplines(degree = 2, knots = c(.2, .4, .6, .8)) + yob,
             m1 = ~ uSplines(degree = 3, knots = c(.1, .3, .5, .7))*yob,
             propensity = morekids ~ samesex,
             noisy = FALSE)
r <- do.call(ivmte, args)

## ---- bounds-------------------------------------------------------------
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

## ---- conventional.tgt.params--------------------------------------------
args <- list(data = AE,
             ivlike =  worked ~ morekids + samesex + morekids*samesex,
             target = "att",
             m0 = ~ u + I(u^2) + yob + u*yob,
             m1 = ~ u + I(u^2) + I(u^3) + yob + u*yob,
             propensity = morekids ~ samesex,
             noisy = FALSE)
r <- do.call(ivmte, args)
args[["target"]] <- "ate"
r <- do.call(ivmte, args)

## ---- late---------------------------------------------------------------
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

## ---- late.cond----------------------------------------------------------
args[["late.X"]] = c(x = 2)
r <- do.call(ivmte, args)
args[["late.X"]] = c(x = 8)
r <- do.call(ivmte, args)

## ---- genlate------------------------------------------------------------
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
args[["genlate.ub"]] <- .41
r <- do.call(ivmte, args)
args[["genlate.ub"]] <- .42
r <- do.call(ivmte, args)

## ---- genlate.cond-------------------------------------------------------
args[["late.X"]] <- c(x = 2)
r <- do.call(ivmte, args)

## ---- custom.weights-----------------------------------------------------
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

## ---- multi.ivlike-------------------------------------------------------
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

## ---- ivlike.components--------------------------------------------------
args[["components"]] <- l(c(intercept, x), c(d), )
r <- do.call(ivmte, args)

## ---- l.example, error = TRUE--------------------------------------------
args[["components"]] <- list(c(intercept, x), c(d), )
args[["components"]] <- list(c("intercept", "x"), c("d"), "")
r <- do.call(ivmte, args)

## ---- ivlike.subsets-----------------------------------------------------
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

## ---- pscore-------------------------------------------------------------
results <- ivmte(data = AE,
                 target = "att",
                 m0 = ~ u + yob,
                 m1 = ~ u + yob,
                 ivlike = worked ~ morekids + samesex + morekids*samesex,
                 propensity = morekids ~ samesex + yob,
                 link = "probit",
                 noisy = FALSE)

## ---- point.id-----------------------------------------------------------
args <- list(data = ivmteSimData,
             ivlike =  y ~ d + factor(z),
             target = "ate",
             m0 = ~ u,
             m1 = ~ u,
             propensity = d ~ factor(z),
             noisy = FALSE)
r <- do.call(ivmte, args)
args[["point"]] <- TRUE
r <- do.call(ivmte, args)

## ---- partial.ci---------------------------------------------------------
results <- ivmte(data = AE,
                 target = "att",
                 m0 = ~ u + yob,
                 m1 = ~ u + yob,
                 ivlike = worked ~ morekids + samesex + morekids*samesex,
                 propensity = morekids ~ samesex,
                 bootstraps = 50,
                 noisy = FALSE)

## ---- point.id.bootstrap-------------------------------------------------
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

## ---- misspecification.test----------------------------------------------
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

args[["ivlike"]] <- y ~ d + factor(z) + d*factor(z) # many more moments
args[["point"]] <- TRUE
r <- do.call(ivmte, args)

