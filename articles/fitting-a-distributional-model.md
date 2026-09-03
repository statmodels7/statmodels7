# Fitting a distributional model

An ordinary regression models the mean and treats everything else about
the distribution as a nuisance held constant. A distributional model
gives **every parameter its own equation**, so the scale, the shape and
a zero-inflation probability may each depend on covariates, and each may
be a smooth, a random effect or a penalized block.

[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
reads one formula carrying all of them, the equations separated by a
bar, and fits it. This vignette walks through a fit and the questions a
reader asks of one.

``` r

n <- 300
da <- data.frame(x = sort(runif(n, -3, 3)))
da$y <- rnorm(n, 1 + 0.8 * da$x, exp(-0.6 + 0.35 * da$x))
```

The scale of these data grows with `x`, from about $`0.19`$ at
$`x = -3`$ to about $`1.57`$ at $`x = 3`$.

## A first fit

The first equation belongs to the distribution’s first parameter, so a
formula with no bar models the mean and leaves the rest as intercepts:

``` r

f1 <- statmod(y ~ x, gaussian1_distrib(), da)
coef(f1)
#> $mu
#> (Intercept)           x 
#>   0.9823413   0.7605631 
#> 
#> $sigma
#> (Intercept) 
#>  -0.2526842
```

[`coef()`](https://rdrr.io/r/stats/coef.html) answers with **one
compartment per parameter**, because the coefficients of the mean and of
the scale are different quantities on different scales and a single flat
vector would run them together. Where the model is an ordinary
regression, it is one:

``` r

all.equal(unname(coef(f1)$mu), unname(coef(lm(y ~ x, da))), tolerance = 1e-6)
#> [1] TRUE
```

## Modeling the scale

Name a parameter on the left of its own equation, after a bar:

``` r

f2 <- statmod(y ~ x | sigma ~ x, gaussian1_distrib(), da)
lapply(coef(f2), round, 3)
#> $mu
#> (Intercept)           x 
#>       1.008       0.795 
#> 
#> $sigma
#> (Intercept)           x 
#>      -0.579       0.356
```

against a truth of $`\mu = 1 + 0.8x`$ and $`\log\sigma = -0.6 + 0.35x`$.
The scale’s equation is on the **log scale**, which is that parameter’s
own link: an equation is fitted on the scale that keeps its parameter
inside its domain, so no combination of covariates can produce a
negative standard deviation.

Comparing the two fits shows what the second equation bought:

``` r

c(constant_scale = loglik(f1), modeled_scale = loglik(f2))
#> constant_scale  modeled_scale 
#>      -349.8763      -245.7099
```

## Terms

Any `modelterms7` term goes in any equation. A smooth carries a penalty,
and its smoothing parameter is estimated, not set:

``` r

db <- data.frame(z = sort(runif(n, -3, 3)),
                 g = factor(rep(1:20, length.out = n)))
db$y <- rnorm(n, 1 + sin(1.4 * db$z) + rnorm(20, sd = 0.6)[db$g], 0.4)

f3 <- statmod(y ~ s(z, k = 8), gaussian1_distrib(), db)
f3@edf
#>   parameter        term coefficients     edf
#> 1        mu      linpar            1 1.00000
#> 2        mu s(z, k = 8)            7 5.27609
#> 3     sigma      linpar            1 1.00000
```

The smooth carries 7 coefficients and **spends** about 5.3 of them: that
is what the penalty buys, and it is the number a criterion should count,
and never the 7.
[`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
reports the hyperparameters and who chose each:

``` r

hyper(f3)
#>   parameter        term   name estimate  held source   id
#> 1        mu s(z, k = 8) lambda 1.067872 FALSE   reml <NA>
```

A random effect is written the same way and is the same machinery, its
penalty being a Gaussian prior over the group effects:

``` r

f4 <- statmod(y ~ s(z, k = 8) + random(~ 1 | g), gaussian1_distrib(), db)
f4@edf
#>   parameter           term coefficients      edf
#> 1        mu         linpar            1  1.00000
#> 2        mu    s(z, k = 8)            7  6.12777
#> 3        mu random(~1 | g)           20 18.57777
#> 4     sigma         linpar            1  1.00000
hyper(f4)
#>   parameter           term   name  estimate  held source   id
#> 1        mu    s(z, k = 8) lambda 1.0454567 FALSE   reml <NA>
#> 2        mu random(~1 | g)  sigma 0.7061116 FALSE   reml <NA>
```

20 group effects spend about 18.6 degrees of freedom here, and the
estimated prior scale is close to the 0.6 the groups were drawn at. A
random intercept is a ridge penalty on an indicator block, which is why
it needs no separate machinery.

## Who chooses the hyperparameters

Every hyperparameter is estimated unless the term holds it. Which
criterion does the choosing is `outer_criterion`:

- [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  and
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  are marginal criteria: they integrate out what the penalty shrinks and
  are differentiable, so the search uses an exact gradient and, where it
  exists, an exact Hessian.
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  is the default.
- [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  and
  [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  score an estimate of prediction error instead, and are what a penalty
  with a **kink** needs: a marginal criterion expands around the
  penalized mode, and for a lasso that mode sits *at* the kink, where
  the second derivative it wants does not exist.

A term holds a hyperparameter by naming it, and
[`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
then says so in its `held` column.

## Reading a fit

``` r

summary(f3)
#> A statmod fit
#> 
#> Call:  statmod(formula = y ~ s(z, k = 8), distrib = gaussian1_distrib(), 
#>             data = db)
#> 
#> Distribution: gaussian1     Observations: 300
#> 
#> === mu   [identity link]
#> 
#> Parametric terms
#>                estimate      se     z       p  lower upper
#>   (Intercept)    0.9152 0.04678 19.56 < 1e-16 0.8235 1.007
#> 
#> s(z, k = 8)   [7 coefficients, edf 5.28]
#>                  estimate      se     z        p   lower  upper
#>   lambda [reml]    1.0680 0.71220                0.28890 3.9470
#>   lin              0.1389 0.04686 2.964 0.003039 0.04704 0.2307
#> 
#> === sigma   [log link]
#> 
#> Parametric terms
#>                estimate      se      z         p   lower   upper
#>   (Intercept)   -0.2103 0.04082 -5.152 2.574e-07 -0.2904 -0.1303
#> 
#> 95% intervals, bayesian variance
#> conditional log-likelihood -362.579490    effective df 7.28
#> cAIC 739.711    cBIC 766.660
#> fitted in 765 ms   search: converged
#> certificate: CONVERGED   outer gradient 1.87e-05   5.77e-13 above the mode
#> 1 note: print(summary(fit), notes = TRUE)
```

Three things in that output are worth naming. The log-likelihood is the
**conditional** one, read at the fitted coefficients, and the effective
degrees of freedom go with it: `cAIC` and `cBIC` are built from that
pair, and comparing a conditional log-likelihood against a marginal one
from another package is comparing two different quantities. Every
interval is built on the scale that keeps its quantity in its own set
and mapped back, so a scale’s interval cannot contain a negative number.
And the last line is the certificate.

## The certificate

The optimizer’s convergence flag answers whether a stopping rule fired.
[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
answers a different question: is the point reported the one the model
asks for?

``` r

cert <- statmod_certificate(f4)
c(state = cert$state, mode_error = signif(cert$mode_error, 3))
#>       state  mode_error 
#> "converged"  "5.91e-11"
```

`mode_error` is how far above its own penalized mode the inner fit
stopped, in log-likelihood units, so one limit serves every model
whatever its scale. The `state` is `converged`, `boundary` where a
hyperparameter has run to an edge and the criterion is genuinely flat
there, or `not converged`.

## Prediction, and the ordinary generics

[`predict()`](https://rdrr.io/r/stats/predict.html) answers per
parameter, with a standard error on request:

``` r

nd <- data.frame(z = c(-1, 0, 1),
                 g = factor(c(1, 2, 3), levels = levels(db$g)))
predict(f4, newdata = nd, se = TRUE)$mu
#>          fit        se      lower       upper
#> 1 -0.2844637 0.1137421 -0.5073941 -0.06153321
#> 2  1.8482711 0.1147791  1.6233081  2.07323404
#> 3  2.4667234 0.1167143  2.2379676  2.69547925
```

Prediction **reapplies** each term’s blueprint rather than rebuilding
it, so the knots of a smooth and the levels of a factor are the ones the
model was fitted with. The identity that pins it is that predicting on
rows the model saw returns the fitted values there.

The `stats` generics behave as they read:

``` r

c(nobs = nobs(f4), df.residual = df.residual(f4))
#>        nobs df.residual 
#>    300.0000    273.2945
length(sigma(f4))
#> [1] 300
```

[`sigma()`](https://rdrr.io/r/stats/sigma.html) returns a **vector**,
because the whole point of the layer is that a scale may be modeled: a
single residual standard deviation exists only where that equation is an
intercept, and returning its first value everywhere else would answer a
different question.
[`df.residual()`](https://rdrr.io/r/stats/df.residual.html) subtracts
the *effective* count, so it is not an integer.

Three generics signal an error instead of guessing.
[`terms()`](https://rdrr.io/r/stats/terms.html) would have to choose one
equation’s terms of several;
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) has no frame
to return, a fit keeping each term’s blueprint and never the data; and
[`anova()`](https://rdrr.io/r/stats/anova.html) has no null distribution
to test against when the hyperparameters were chosen from the same data.

## Simulating from a model

[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
draws from a model at coefficients you supply, as a simulation study
needs:

``` r

sim <- rstatmod(y ~ x, gaussian1_distrib(), da,
                par = list(mu = c(1, 0.8), sigma = -0.6))
names(sim)
#> [1] "data"       "par"        "theta"      "latent"     "structural"
#> [6] "n_sim"      "call"
head(sim$data, 3)
#>           x          y
#> 1 -2.921535 -1.6210728
#> 2 -2.919658 -0.3713365
#> 3 -2.860013 -1.1992527
```

It returns the truth **beside** the data and never as an attribute of
it: an attribute survives a row subset without being subset itself, so
`attr(sim[1:10, ], "theta")` would still hold every row’s parameters and
be silently misaligned, which for a simulation study is the worse
failure.

## Summary

- One formula, the equations separated by a bar, one per distribution
  parameter; an equation omitted is an intercept.
- Each equation is fitted on its parameter’s own link scale, so no
  combination of covariates leaves the parameter’s domain.
- [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) answer with one
  compartment per parameter.
- Any `modelterms7` term goes in any equation, and its penalty travels
  with it.
- Hyperparameters are estimated unless a term holds one:
  [`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  and
  [`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
  where the penalty is twice differentiable,
  [`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md),
  [`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md)
  or [`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
  where it has a kink.
- The log-likelihood reported is the conditional one, and `cAIC` and
  `cBIC` go with the effective degrees of freedom, never a parameter
  count.
- [`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
  asks whether the point is the model’s mode, which the optimizer’s flag
  does not.
