# Fitting a distributional model

An ordinary regression models the mean and treats everything else about
the distribution as a nuisance held constant. A distributional model
gives every parameter its own equation, so the scale, the shape and a
zero-inflation probability may each depend on covariates, and each may
be a smooth, a random effect or a penalized block.

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
reads one formula that carries all of them, in which a bar separates the
equations, and fits it. This vignette walks through a fit and then
through the functions that read it.

``` r

n <- 300
da <- data.frame(x = sort(runif(n, min = -3, max = 3)))
da$y <- rnorm(n, mean = 1 + 0.8 * da$x, sd = exp(-0.6 + 0.35 * da$x))
```

The scale of these data grows with `x`, from about $`0.19`$ at
$`x = -3`$ to about $`1.57`$ at $`x = 3`$.

## A first fit

The first equation belongs to the first parameter of the distribution,
so a formula with no bar models the mean and leaves the other parameters
as intercepts:

``` r

f1 <- statmod(y ~ x, distrib = gaussian1_distrib(), data = da)
coef(f1)
#> $mu
#> (Intercept)           x 
#>   0.9823413   0.7605631 
#> 
#> $sigma
#> (Intercept) 
#>  -0.2493397
```

[`coef()`](https://rdrr.io/r/stats/coef.html) returns one compartment
per parameter, because the coefficients of the mean and of the scale are
different quantities on different scales, and a single flat vector would
run them together. For this model the compartment of the mean holds the
coefficients of an ordinary regression, which the following check
confirms:

``` r

all.equal(unname(coef(f1)$mu), unname(coef(lm(y ~ x, data = da))),
          tolerance = 1e-6)
#> [1] TRUE
```

## Modeling the scale

A parameter is named on the left of its own equation, after a bar:

``` r

f2 <- statmod(y ~ x | sigma ~ x, distrib = gaussian1_distrib(), data = da)
lapply(coef(f2), round, digits = 3)
#> $mu
#> (Intercept)           x 
#>       1.008       0.795 
#> 
#> $sigma
#> (Intercept)           x 
#>      -0.575       0.356
```

The truth is $`\mu = 1 + 0.8x`$ and $`\log\sigma = -0.6 + 0.35x`$. The
equation of the scale is on the log scale, which is the link of that
parameter: an equation is fitted on the scale that keeps its parameter
inside its domain, so a combination of covariates cannot produce a
negative standard deviation.

The log-likelihoods of the two fits show the gain from the second
equation:

``` r

c(constant_scale = loglik(f1), modeled_scale = loglik(f2))
#> constant_scale  modeled_scale 
#>      -349.8797      -245.7139
```

## Terms

Any `modelterms7` term goes in any equation. A smooth carries a penalty,
and its smoothing parameter is estimated from the data:

``` r

db <- data.frame(z = sort(runif(n, min = -3, max = 3)),
                 g = factor(rep(1:20, length.out = n)))
db$y <- rnorm(n, mean = 1 + sin(1.4 * db$z) + rnorm(20, sd = 0.6)[db$g], sd = 0.4)

f3 <- statmod(y ~ s(z, bspline_smooth(k = 8)),
              distrib = gaussian1_distrib(), data = db)
f3@edf
#>   parameter                        term coefficients      edf
#> 1        mu                      linpar            1 1.000000
#> 2        mu s(z, bspline_smooth(k = 8))            7 5.264483
#> 3     sigma                      linpar            1 1.000000
```

The smooth has 7 coefficients and uses about 5.3 of them. This effective
number of degrees of freedom is the one that a criterion counts, in
place of the number of coefficients.
[`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
reports the hyperparameters and the source of each:

``` r

hyper(f3)
#>   parameter                        term   name estimate  held source   id
#> 1        mu s(z, bspline_smooth(k = 8)) lambda 1.064316 FALSE   reml <NA>
```

A random effect is written the same way and is fitted by the same
machinery, and its penalty is a Gaussian prior over the group effects:

``` r

f4 <- statmod(y ~ s(z, bspline_smooth(k = 8)) + random(~ 1 | g),
              distrib = gaussian1_distrib(), data = db)
f4@edf
#>   parameter                        term coefficients       edf
#> 1        mu                      linpar            1  1.000000
#> 2        mu s(z, bspline_smooth(k = 8))            7  6.083231
#> 3        mu              random(~1 | g)           20 18.540240
#> 4     sigma                      linpar            1  1.000000
hyper(f4)
#>   parameter                        term   name  estimate  held source   id
#> 1        mu s(z, bspline_smooth(k = 8)) lambda 1.0345354 FALSE   reml <NA>
#> 2        mu              random(~1 | g)  sigma 0.7070243 FALSE   reml <NA>
```

The 20 group effects spend about 18.5 degrees of freedom here, and the
estimated prior scale is 0.71, against the 0.6 used to draw the groups.
A random intercept is a ridge penalty on an indicator block, so it is
fitted like a smooth.

## Who chooses the hyperparameters

Every hyperparameter is estimated unless the term holds it. The
criterion is set by the arguments `outer_criterion` and
`sparse_criterion` of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md).
The argument `outer_criterion` applies to the hyperparameters of a
penalty that is twice differentiable, such as a smooth or a random
effect:

- [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  and
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  are marginal criteria. They integrate out the penalized coefficients
  and are differentiable, so the search uses an exact gradient and,
  where it exists, an exact Hessian.
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  is the default.
- [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  and
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  score an estimate of prediction error instead.

A penalty with a kink (`lasso()`, `enet()`, `scad()` and `mcp()`) is
chosen along a path of its own values, from the value that empties the
block downwards, and the path is scored by `sparse_criterion`. This
argument is
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) by
default and may be
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md), and
a marginal criterion is rejected there. A marginal criterion expands the
penalized objective around the penalized mode, and for a lasso that mode
lies *at* the kink, where the second derivative does not exist. The
function
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
serves only as `sparse_criterion`.

A term holds a hyperparameter when the call names its value, for example
`lasso(x, lambda = 3)` or
`s(z, bspline_smooth(k = 8), hyper = c(lambda = 2))`, and the `held`
column of
[`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
then reports `TRUE`:

``` r

f5 <- statmod(y ~ s(z, bspline_smooth(k = 8), hyper = c(lambda = 2)),
              distrib = gaussian1_distrib(), data = db)
hyper(f5)
#>   parameter                                               term   name estimate
#> 1        mu s(z, bspline_smooth(k = 8), hyper = c(lambda = 2)) lambda        2
#>   held source   id
#> 1 TRUE  fixed <NA>
```

## Reading a fit

``` r

summary(f3)
#> A statmod fit
#> 
#> Call:  statmod(formula = y ~ s(z, bspline_smooth(k = 8)), distrib = gaussian1_distrib(), 
#>             data = db)
#> 
#> Distribution: gaussian1     Observations: 300
#> 
#> === mu   [identity link]
#> 
#> Parametric terms
#>                estimate      se     z       p  lower upper
#>   (Intercept)    0.9152 0.04728 19.36 < 1e-16 0.8225 1.008
#> 
#> s(z, bspline_smooth(k = 8))   [7 coefficients, edf 5.26]
#>                  estimate      se     z        p   lower  upper
#>   lambda [reml]     1.064  0.7107                 0.2875   3.94
#>   lin              0.1389 0.04736 2.933 0.003361 0.04606 0.2317
#> 
#> === sigma   [log link]
#> 
#> Parametric terms
#>                estimate      se      z         p   lower   upper
#>   (Intercept)   -0.1998 0.04124 -4.844 1.273e-06 -0.2806 -0.1189
#> 
#> 95% intervals, bayesian variance
#> conditional log-likelihood -362.616546    effective df 7.26
#> cAIC 739.762    cBIC 766.668
#> fitted in 258 ms
#> inner   max |grad|/se 3.7e-16   min eigen 1
#> outer   max |grad|/se 2.1e-08   min eigen 0.99
#> 1 note: print(summary(fit), notes = TRUE)
```

The log-likelihood is the conditional one, read at the fitted
coefficients, and the effective degrees of freedom go with it: `cAIC`
and `cBIC` are built from that pair, and comparing a conditional
log-likelihood with a marginal one from another package compares two
different quantities. The interval of a hyperparameter, such as a
smoothing parameter or a prior standard deviation, is built on the scale
that keeps it positive and mapped back, so it cannot contain a negative
number, while the interval of a coefficient is on the link scale of its
equation. The last lines, headed `inner` and `outer`, give the largest
gradient in standard errors and the smallest eigenvalue of the curvature
for the coefficients and for the hyperparameters, and both are read at
the point that the fit reached.

## The certificate

The convergence flag of an optimizer records whether a stopping rule
fired.
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
tests a different property: whether the point reported is a mode of the
penalized objective and a stationary point of the criterion in the
hyperparameters.

``` r

cert <- statmod_certificate(f4)
c(state = cert$state, mode_error = signif(cert$mode_error, digits = 3))
#>       state  mode_error 
#> "converged"  "4.93e-30"
```

`mode_error` is how far above its own penalized mode the inner fit
stopped, in log-likelihood units, so one limit serves every model
whatever its scale. The `state` is `converged`; `boundary`, where a
hyperparameter or an estimated coefficient lies at an edge of its range
with the criterion flat in it, or the outer curvature has a flat
direction; `not converged`; or `unknown`, where the derivatives that the
certificate needs cannot be computed.

## Prediction, and the ordinary generics

[`predict()`](https://rdrr.io/r/stats/predict.html) returns one
compartment per parameter, with a standard error on request:

``` r

nd <- data.frame(z = c(-1, 0, 1),
                 g = factor(c(1, 2, 3), levels = levels(db$g)))
predict(f4, newdata = nd, se = TRUE)$mu
#>          fit        se      lower       upper
#> 1 -0.2837629 0.1188419 -0.5166888 -0.05083693
#> 2  1.8468308 0.1199200  1.6117920  2.08186963
#> 3  2.4646126 0.1219466  2.2256016  2.70362367
```

Prediction reapplies each term’s blueprint instead of rebuilding it, so
the knots of a smooth and the levels of a factor are those of the fitted
model. As a consequence, predicting at rows that were used to fit the
model returns the fitted values there.

The generics of `stats` apply to a fit:

``` r

c(nobs = nobs(f4), df.residual = df.residual(f4))
#>        nobs df.residual 
#>    300.0000    273.3765
length(sigma(f4))
#> [1] 300
```

[`sigma()`](https://rdrr.io/r/stats/sigma.html) returns the standard
deviation of the response at each observation, so it is a vector
wherever the scale is modeled or the variance depends on the mean.
[`df.residual()`](https://rdrr.io/r/stats/df.residual.html) subtracts
the *effective* count of degrees of freedom from the number of
observations, so it is not an integer.

Three generics signal an error, and each error message names the
function to use instead. [`terms()`](https://rdrr.io/r/stats/terms.html)
fails because a fit has one set of terms per distribution parameter,
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) fails
because a fit keeps the blueprint of each term and not the data, and
[`anova()`](https://rdrr.io/r/stats/anova.html) fails because a
penalized fit whose hyperparameters were chosen from the same data has
no null distribution against which to test.

## Simulating from a model

[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
draws from a model at coefficients supplied by the caller, as a
simulation study needs:

``` r

sim <- rstatmod(y ~ x, distrib = gaussian1_distrib(), data = da,
                par = list(mu = c(1, 0.8), sigma = -0.6))
names(sim)
#> [1] "data"       "par"        "theta"      "latent"     "structural"
#> [6] "hyper"      "n_sim"      "call"
head(sim$data, n = 3)
#>           x          y
#> 1 -2.921535 -1.6210728
#> 2 -2.919658 -0.3713365
#> 3 -2.860013 -1.1992527
```

It returns the true parameters beside the data and not as an attribute
of the data. An attribute survives a row subset without being subset
itself, so if the parameters were stored in `attr(sim$data, "theta")`,
the result of `attr(sim$data[1:10, ], "theta")` would still hold every
row’s parameters and would be misaligned with the rows kept.

## Summary

- One formula, in which the equations are separated by a bar, one per
  distribution parameter; an equation that is omitted is an intercept.
- Each equation is fitted on the link scale of its parameter, so a
  combination of covariates cannot leave the domain of the parameter.
- [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) return one
  compartment per parameter.
- Any `modelterms7` term goes in any equation, and its penalty travels
  with it.
- Hyperparameters are estimated unless a term holds one: by
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  or
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  (or
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  or
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md))
  as `outer_criterion` where the penalty is twice differentiable, and
  along a path scored by `sparse_criterion`
  ([`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  by default, or
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  or
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md))
  where it has a kink.
- The log-likelihood reported is the conditional one, and `cAIC` and
  `cBIC` use the effective degrees of freedom and not a parameter count.
- [`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
  tests whether the point is a mode of the model, which the flag of the
  optimizer does not do.
