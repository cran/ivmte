## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE,
  fig.path = "vignettes/ivmte_files/figure-gfm/"
)
require(pander)
require(data.table)
require(AER)
require(ivmte)
require(ggplot2)
require(gridExtra)
require(splines2)

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
#                   propensity = morekids ~ samesex,
#                   noisy = TRUE)

## ---- syntaxrun----------------------------------------------------------
results <- ivmte(data = AE,
                 target = "att",
                 m0 = ~ u + yob,
                 m1 = ~ u + yob,
                 ivlike = worked ~ morekids + samesex + morekids*samesex,
                 propensity = morekids ~ samesex, 
                 noisy = TRUE)

## ---- syntaxrun.quiet----------------------------------------------------
results <- ivmte(data = AE,
                 target = "att",
                 m0 = ~ u + yob,
                 m1 = ~ u + yob,
                 ivlike = worked ~ morekids + samesex + morekids*samesex,
                 propensity = morekids ~ samesex,
                 noisy = FALSE)
results
cat(results$messages, sep = "\n")

## ---- mtrbasics----------------------------------------------------------
args <- list(data = AE,
             ivlike =  worked ~ morekids + samesex + morekids*samesex,
             target = "att",
             m0 = ~ u + I(u^2) + yob + u*yob,
             m1 = ~ u + I(u^2) + I(u^3) + yob + u*yob,
             propensity = morekids ~ samesex)
r <- do.call(ivmte, args)
r

## ---- non.polynomial.u.error, error = TRUE-------------------------------
args[["m0"]] <- ~ log(u) + yob
r <- do.call(ivmte, args)

## ---- uname--------------------------------------------------------------
args[["m0"]] <- ~ v + I(v^2) + yob + v*yob
args[["m1"]] <- ~ v + I(v^2) + I(v^3) + yob + v*yob
args[["uname"]] <- "v"
r <- do.call(ivmte, args)
r

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
r

## ---- spline-------------------------------------------------------------
args <- list(data = AE,
             ivlike =  worked ~ morekids + samesex + morekids*samesex,
             target = "att",
             m0 = ~ u + uSplines(degree = 1, knots = c(.2, .4, .6, .8)) + yob,
             m1 = ~ uSplines(degree = 2, knots = c(.1, .3, .5, .7))*yob,
             propensity = morekids ~ samesex)
r <- do.call(ivmte, args)
r

## ---- bounds-------------------------------------------------------------
args <- list(data = AE,
             ivlike =  worked ~ morekids + samesex + morekids*samesex,
             target = "att",
             m0 = ~ u + uSplines(degree = 1, knots = c(.2, .4, .6, .8)) + yob,
             m1 = ~ uSplines(degree = 2, knots = c(.1, .3, .5, .7))*yob,
             m1.inc = TRUE,
             m0.inc = TRUE,
             mte.dec = TRUE,
             propensity = morekids ~ samesex)
r <- do.call(ivmte, args)
r

## ---- conventional.tgt.params--------------------------------------------
args <- list(data = AE,
             ivlike =  worked ~ morekids + samesex + morekids*samesex,
             target = "att",
             m0 = ~ u + I(u^2) + yob + u*yob,
             m1 = ~ u + I(u^2) + I(u^3) + yob + u*yob,
             propensity = morekids ~ samesex)
r <- do.call(ivmte, args)
r
args[["target"]] <- "ate"
r <- do.call(ivmte, args)
r

## ---- late---------------------------------------------------------------
args <- list(data = ivmteSimData,
             ivlike =  y ~ d + z + d*z,
             target = "late",
             late.from = c(z = 1),
             late.to = c(z = 3),
             m0 = ~ u + I(u^2) + I(u^3) + x,
             m1 = ~ u + I(u^2) + I(u^3) + x,
             propensity = d ~ z + x)
r <- do.call(ivmte, args)
r

## ---- late.cond----------------------------------------------------------
args[["late.X"]] = c(x = 2)
r <- do.call(ivmte, args)
r
args[["late.X"]] = c(x = 8)
r <- do.call(ivmte, args)
r

## ---- genlate------------------------------------------------------------
args <- list(data = ivmteSimData,
             ivlike =  y ~ d + z + d*z,
             target = "genlate",
             genlate.lb = .2,
             genlate.ub = .4,
             m0 = ~ u + I(u^2) + I(u^3) + x,
             m1 = ~ u + I(u^2) + I(u^3) + x,
             propensity = d ~ z + x)
r <- do.call(ivmte, args)
args[["genlate.ub"]] <- .41
r <- do.call(ivmte, args)
args[["genlate.ub"]] <- .42
r <- do.call(ivmte, args)
r

## ---- genlate.cond-------------------------------------------------------
args[["late.X"]] <- c(x = 2)
r <- do.call(ivmte, args)
r

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
             propensity = d ~ z + x)
r <- do.call(ivmte, args)
r

## ---- multi.ivlike-------------------------------------------------------
args <- list(data = ivmteSimData,
             ivlike =  c(y ~ (z == 1) + (z == 2) + (z == 3) + x,
                         y ~ d + x,
                         y ~ d | z),
             target = "ate",
             m0 = ~ uSplines(degree = 1, knots = c(.25, .5, .75)) + x,
             m1 = ~ uSplines(degree = 1, knots = c(.25, .5, .75)) + x,
             propensity = d ~ z + x)
r <- do.call(ivmte, args)
r

## ---- ivlike.components--------------------------------------------------
args[["components"]] <- l(c(intercept, x), c(d), )
r <- do.call(ivmte, args)
r

## ---- l.example, error = TRUE--------------------------------------------
args[["components"]] <- list(c(intercept, x), c(d), )
args[["components"]] <- list(c("intercept", "x"), c("d"), "")
r <- do.call(ivmte, args)
r

## ---- ivlike.subsets-----------------------------------------------------
args <- list(data = ivmteSimData,
             ivlike =  c(y ~ z + x,
                         y ~ d + x,
                         y ~ d | z),
             subset = l(x <= 9, 1 == 1, z %in% c(1,3)),
             target = "ate",
             m0 = ~ uSplines(degree = 3, knots = c(.25, .5, .75)) + x,
             m1 = ~ uSplines(degree = 3, knots = c(.25, .5, .75)) + x,
             propensity = d ~ z + x)
r <- do.call(ivmte, args)
r

## ---- pscore-------------------------------------------------------------
results <- ivmte(data = AE,
                 target = "att",
                 m0 = ~ u + yob,
                 m1 = ~ u + yob,
                 ivlike = worked ~ morekids + samesex + morekids*samesex,
                 propensity = morekids ~ samesex + yob,
                 link = "probit")
results

## ---- point.id-----------------------------------------------------------
args <- list(data = ivmteSimData,
             ivlike =  y ~ d + factor(z),
             target = "ate",
             m0 = ~ u,
             m1 = ~ u,
             propensity = d ~ factor(z))
r <- do.call(ivmte, args)
args[["point"]] <- TRUE
r <- do.call(ivmte, args)
r

## ---- partial.ci---------------------------------------------------------
r <- ivmte(data = AE,
           target = "att",
           m0 = ~ u + yob,
           m1 = ~ u + yob,
           ivlike = worked ~ morekids + samesex + morekids*samesex,
           propensity = morekids ~ samesex,
           bootstraps = 50)
summary(r)

## ---- partial.ci.show----------------------------------------------------
r$bounds.ci

## ---- bootstrap.data, eval = TRUE----------------------------------------
head(r$bounds.bootstrap)

## ---- bootstrap.plot, eval = TRUE, echo = FALSE--------------------------
bootstraps <- data.frame(type = rep(c('lb', 'ub'), each = 50),
                         value = c(r$bounds.bootstraps[, 1],
                                   r$bounds.bootstraps[, 2]))
ggplot(data = bootstraps,
       aes(x = value, color = type, fill = type)) + 
    geom_histogram(position = "dodge", alpha = 0.5, binwidth = 0.05) + 
    geom_vline(xintercept = r$bounds[1],
               linetype = "dashed", size = 1.5, color = "orange") +
    geom_vline(xintercept = r$bounds[2], 
               linetype = "dashed", size = 1.5, color = "lightskyblue") +
    ## Labeling options
    labs(x = "Bound",
         y = "Frequency") +
    ## Presentation options
    theme(text = element_text(size = 10),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size=10,
                                   angle = 0),
          panel.background = element_blank(),
          legend.key = element_blank(),
          legend.position = "bottom",
          legend.box = "vertical",
          legend.title = element_blank(),
          legend.text = element_text(size = 10)) +
    scale_color_manual("", values = c("orange", "lightskyblue"),
                       labels = c(" Lower bound   ", " Upper bound ")) +
    scale_fill_manual("", values = c("orange", "lightskyblue"),
                      labels = c(" Lower bound   ", " Upper bound "))

## ---- point.id.bootstrap-------------------------------------------------
args <- list(data = AE,
             target = "att",
             m0 = ~ u,
             m1 = ~ u,
             ivlike = worked ~ morekids + samesex + morekids*samesex,
             propensity = morekids ~ samesex,
             point = TRUE,
             bootstraps = 50)
r <- do.call(ivmte, args)
summary(r)

## ---- misspecification.test----------------------------------------------
args <- list(data = ivmteSimData,
             ivlike =  y ~ d + factor(z),
             target = "ate",
             m0 = ~ u,
             m1 = ~ u,
             m0.dec = TRUE,
             m1.dec = TRUE,
             propensity = d ~ factor(z),
             bootstraps = 50)
r <- do.call(ivmte, args)
r

args[["ivlike"]] <- y ~ d + factor(z) + d*factor(z) # many more moments
args[["point"]] <- TRUE
r <- do.call(ivmte, args)
r

## ---- draw.specification, eval = TRUE------------------------------------
args <- list(data = AE,
             ivlike =  worked ~ morekids + samesex + morekids*samesex,
             target = "att",
             m0 = ~ 0 + uSplines(degree = 2, knots = c(1/3, 2/3)),
             m1 = ~ 0 + uSplines(degree = 2, knots = c(1/3, 2/3)),
             m1.inc = TRUE,
             m0.inc = TRUE,
             mte.dec = TRUE,
             propensity = morekids ~ samesex)
r <- do.call(ivmte, args)
r

## ---- draw.m0.splines, eval = TRUE---------------------------------------
specs0 <- r$splines.dict$m0[[1]]
specs0

## ---- draw.m0.coef, eval = TRUE------------------------------------------
r$gstar.coef$min.g0

## ---- draw.design.m0, eval = TRUE----------------------------------------
uSeq <- seq(0, 1, by = 0.01)
dmat0 <- bSpline(x = uSeq,
                      degree = specs0$degree,
                      intercept = specs0$intercept,
                      knots = specs0$knots,
                      Boundary.knots = c(0, 1))
m0.min <- dmat0 %*% r$gstar.coef$min.g0
m0.max <- dmat0 %*% r$gstar.coef$max.g0

## ---- draw.design.m1, eval = TRUE, echo = FALSE--------------------------
specs1 <- r$splines.dict$m1[[1]]
dmat1 <- bSpline(x = uSeq,
                 degree = specs1$degree,
                 intercept = specs1$intercept,
                 knots = specs1$knots,
                 Boundary.knots = c(0, 1))
m1.min <- dmat1 %*% r$gstar.coef$min.g1
m1.max <- dmat1 %*% r$gstar.coef$max.g1

## ---- eval = TRUE--------------------------------------------------------
mte.min <- m1.min - m0.min
mte.max <- m1.max - m0.max

## ---- draw.plots, eval = TRUE, echo = FALSE, fig.width = 8, fig.height = 10----
lim0 <- c(0, 0.8)
lim1 <- c(0.0, 0.8)
limte <- c(-0.2, 0.5)
dt1 <- data.frame(x = uSeq, y = m0.min)
dt2 <- data.frame(x = uSeq, y = m0.max)
dt3 <- data.frame(x = uSeq, y = m1.min)
dt4 <- data.frame(x = uSeq, y = m1.max)
dt5 <- data.frame(x = uSeq, y = mte.min)
dt6 <- data.frame(x = uSeq, y = mte.max)
for (i in 1:6) {
    if (i %in% c(1, 2)) {
        ylim <- lim0
        if (i == 1) ylab <- expression(paste("Min. ", m[0]))
        if (i == 2) ylab <- expression(paste("Max. ", m[0]))
    }
    if (i %in% c(3, 4)) {
        ylim <- lim1
        if (i == 3) ylab <- expression(paste("Min. ", m[1]))
        if (i == 4) ylab <- expression(paste("Max. ", m[1]))
    }
    if (i %in% c(5, 6)) {
        ylim <- limte
        if (i == 5) ylab <- expression(paste("Min. ", MTE(u)))
        if (i == 6) ylab <- expression(paste("Max. ", MTE(u)))
    }
    assign("dt", get(paste0("dt", i)))
    figure <- ggplot() +
        geom_line(data = dt,
                  aes(x = x,
                      y = y),
                  size = 1.0) +
        ## Labeling options
        labs(x = "u",
             y = ylab) +
        scale_x_continuous(limits = c(0, 1),
                           breaks = seq(0, 1, 0.25)) +
        scale_y_continuous(limits = ylim) +
        ## Presentation options
        theme(text = element_text(size = 10),
              axis.line = element_line(color = "black"),
              axis.text = element_text(size=10,
                                       angle = 0),
              panel.background = element_blank(),
              legend.key = element_blank(),
              legend.position = "bottom",
              legend.box = "horizontal",
              legend.title = element_blank(),
              legend.text = element_text(size = 10))
    assign(paste0("figure", i), figure)
}
grid.arrange(figure1, figure2, figure3, figure4, figure5, figure6, ncol=2)

## ---- weights.matrix, eval = TRUE, echo = TRUE---------------------------
## Target weights
w1 <- cbind(r$gstar.weights$w1$lb,
                  r$gstar.weights$w1$ub,
                  r$gstar.weights$w1$multiplier)
## IV-like estimand weights
sw1 <- cbind(r$s.set$s1$w1$lb,
             r$s.set$s1$w1$ub,
             r$s.set$s1$w1$multiplier)
sw2 <- cbind(r$s.set$s2$w1$lb,
             r$s.set$s2$w1$ub,
             r$s.set$s2$w1$multiplier)
sw3 <- cbind(r$s.set$s3$w1$lb,
             r$s.set$s3$w1$ub,
             r$s.set$s3$w1$multiplier)
sw4 <- cbind(r$s.set$s4$w1$lb,
             r$s.set$s4$w1$ub,
             r$s.set$s4$w1$multiplier)
## Assign column names
colnames(w1) <-
    colnames(sw1) <-
    colnames(sw2) <-
    colnames(sw3) <-
    colnames(sw4) <- c("lb", "ub", "mp")

## ---- weights.prop, eval = TRUE, echo = TRUE-----------------------------
pscore <- sort(unique(w1[, "ub"]))
pscore

## ---- weights.data, eval = TRUE, echo = TRUE-----------------------------
avg1 <- NULL ## The data.frame that will contain the average weights
i <- 0 ## An index for the type of weight
for (j in c('w1', 'sw1', 'sw2', 'sw3', 'sw4')) {
    dt <- data.frame(get(j))
    avg <- rbind(
        ## Average for u in [0, pscore[1])
        c(s = i, d = 1, lb = 0, ub = pscore[1],
          avgWeight = mean(dt$mp)),
        ## Average for u in [pscore[1], pscore[2])
        c(s = i, d = 1, lb = pscore[1], ub = pscore[2],
          avgWeight = mean(as.integer(dt$ub ==
                                      pscore[2]) * dt$mp)),
        ## Average for u in [pscore[2], 1]
        c(s = i, d = 1, lb = pscore[2], ub = 1, avgWeight = 0))
    avg1 <- rbind(avg1, avg)
    i <- i + 1
}
avg1 <- data.frame(avg1)

## ---- evaluate = TRUE, echo = FALSE--------------------------------------
knitr::kable(avg1)

## ---- weights.plot1, eval = TRUE, echo = FALSE---------------------------
## Account for overlap
avg1[avg1$lb > 0 & avg1$lb < pscore[2] & avg1$s == 2, ]$avgWeight <- -0.07
avg1[avg1$lb > pscore[1] & avg1$s == 2, ]$avgWeight <- -0.07
avg1[avg1$lb > pscore[1] & avg1$s == 0, ]$avgWeight <- +0.07
avg1[avg1$lb > pscore[1] & avg1$s == 4, ]$avgWeight <- +0.14
## Draw plot
wFigure1 <- ggplot() +
    geom_segment(data = avg1[avg1$s == 0, ],
                 mapping = aes(x=lb, xend=ub, y = avgWeight, yend = avgWeight,
                               color= "Target"), size = 2, linetype = "dashed") +
    geom_segment(data = avg1[avg1$s == 2, ],
                 mapping = aes(x=lb, xend=ub, y = avgWeight, yend = avgWeight,
                               color= "More kids", ), size = 2) +
    geom_segment(data = avg1[avg1$s == 3, ],
                 mapping = aes(x=lb, xend=ub, y = avgWeight, yend = avgWeight,
                               color= "Same sex"), size = 2) +
    geom_segment(data = avg1[avg1$s == 4, ],
                 mapping = aes(x=lb, xend=ub, y = avgWeight, yend = avgWeight,
                               color= "More kids x Same sex"), size = 2) +
    scale_x_continuous(breaks=seq(0, 1, 0.2)) +
    theme(text = element_text(size = 10),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size=10,
                                   angle = 0),
          panel.background = element_blank(),
          legend.key = element_blank(),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.title = element_blank(),
          legend.text = element_text(size = 10)) +
    labs(x="u", y="Average weights, d = 1") +
    scale_color_manual("",
                       breaks = c("Target", "More kids", "Same sex", "More kids x Same sex"),
                       values = c("lightskyblue", "orange", "olivedrab", "gray60"))

## ---- weights.plot0, eval = TRUE, echo = FALSE---------------------------
w0 <- cbind(r$gstar.weights$w0$lb,
                  r$gstar.weights$w0$ub,
                  r$gstar.weights$w0$multiplier)
sw1 <- cbind(r$s.set$s1$w0$lb,
             r$s.set$s1$w0$ub,
             r$s.set$s1$w0$multiplier)
sw2 <- cbind(r$s.set$s2$w0$lb,
             r$s.set$s2$w0$ub,
             r$s.set$s2$w0$multiplier)
sw3 <- cbind(r$s.set$s3$w0$lb,
             r$s.set$s3$w0$ub,
             r$s.set$s3$w0$multiplier)
sw4 <- cbind(r$s.set$s4$w0$lb,
             r$s.set$s4$w0$ub,
             r$s.set$s4$w0$multiplier)
colnames(w0) <-
    colnames(sw1) <-
    colnames(sw2) <-
    colnames(sw3) <-
    colnames(sw4) <- c("lb", "ub", "mp")
## Construct data frame
avg0 <- NULL
i <- 0
for (j in c('w0', 'sw1', 'sw2', 'sw3', 'sw4')) {
    dt <- data.frame(get(j))
    if (j == 'w0') {
        avg <- rbind(c(s = i, d = 1, lb = 0, ub = pscore[1],
                       avgWeight = mean(dt$mp)),
                     c(s = i, d = 1, lb = pscore[1], ub = pscore[2],
                       avgWeight = mean(as.integer(dt$ub ==
                                                   pscore[2]) * dt$mp)),
                     c(s = i, d = 1, lb = pscore[2], ub = 1, avgWeight = 0))
    } else {
        avg <- rbind(c(s = i, d = 0, lb = 0, ub = pscore[1],
                       avgWeight = 0),
                     c(s = i, d = 0, lb = pscore[1], ub = pscore[2],
                       avgWeight = mean(as.integer(dt$ub ==
                                                   pscore[1]) * dt$mp)),
                     c(s = i, d = 0, lb = pscore[2], ub = 1,
                       avgWeight = mean(dt$mp)))
    }
    avg0 <- rbind(avg0, avg)
    i <- i + 1
}
avg0 <- data.frame(avg0)
## Deal with overlaps
avg0[avg0$lb > 0 & avg0$lb < pscore[2] & avg0$s == 2, ]$avgWeight <- -0.055
avg0[avg0$lb > 0 & avg0$lb < pscore[2] & avg0$s == 4, ]$avgWeight <- +0.055
avg0[avg0$lb == 0 & avg0$s == 2, ]$avgWeight <- -0.055
avg0[avg0$lb == 0 & avg0$s == 4, ]$avgWeight <- +0.055
## Draw plot
wFigure2 <- ggplot() +
    geom_segment(data = avg0[avg0$s == 0, ],
                 mapping = aes(x=lb, xend=ub, y = avgWeight, yend = avgWeight,
                               color= "Target"), size = 2, linetype = "dashed") +
    geom_segment(data = avg0[avg0$s == 2, ],
                 mapping = aes(x=lb, xend=ub, y = avgWeight, yend = avgWeight,
                               color= "More kids", ), size = 2) +
    geom_segment(data = avg0[avg0$s == 3, ],
                 mapping = aes(x=lb, xend=ub, y = avgWeight, yend = avgWeight,
                               color= "Same sex"), size = 2) +
    geom_segment(data = avg0[avg0$s == 4, ],
                 mapping = aes(x=lb, xend=ub, y = avgWeight, yend = avgWeight,
                               color= "More kids x Same sex"), size = 2) +
    scale_x_continuous(breaks=seq(0, 1, 0.2)) +
    theme(text = element_text(size = 10),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size=10,
                                   angle = 0),
          panel.background = element_blank(),
          legend.key = element_blank(),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.title = element_blank(),
          legend.text = element_text(size = 10)) +
    labs(x="u", y="Average weights, d = 0") +
    scale_color_manual("",
                       breaks = c("Target", "More kids", "Same sex", "More kids x Same sex"),
                       values = c("lightskyblue", "orange", "olivedrab", "gray60"))

## ---- weights.plot.combine, eval = TRUE, echo = FALSE, fig.width = 8, fig.height = 5----
## Combine plots
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
jointLegend <- g_legend(wFigure2)
grid.arrange(arrangeGrob(wFigure2 + theme(legend.position = "none"),
                         wFigure1 + theme(legend.position = "none"),
                         nrow = 1),
             jointLegend, heights = c(10, 1))

