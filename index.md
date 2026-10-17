# statmodels7

[statmodels7](https://statmodels7.github.io/statmodels7/) is the
modeling layer of the [statmodels7](https://statmodels7.github.io)
toolkit, and the package that installs the rest of it.
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
reads one formula carrying an equation for every parameter of a
distribution and fits it, with any term of
[modelterms7](https://statmodels7.github.io/modelterms7/) in any of
those equations.

## Installation

``` r

# install.packages("pak")
pak::pak("statmodels7/statmodels7")
```

Attaching it attaches the eight member packages.

``` r

library(statmodels7)
#> -- Attaching the statmodels7 toolkit 0.203.3
#> v basis7         0.15.1   v numericals7    0.23.3
#> v distributions7 0.80.2   v optimizers7    0.15.1
#> v linkfunctions7 0.5.4    v parameters7    0.26.0
#> v modelterms7    0.99.3   v penalties7     0.33.0
```

## Fitting a model

The equations are separated by a bar, one for each parameter of the
distribution, and a parameter with no equation of its own is given an
intercept. Below, the mean carries a penalized smooth of `x` and a
random intercept over `g`, and the standard deviation is linear in `x`
on its log scale.

``` r

set.seed(1)
n <- 400
d <- data.frame(x = runif(n, min = -2, max = 2),
                g = factor(rep(1:20, length.out = n)))
u <- rnorm(20, mean = 0, sd = 0.5)
d$y <- sin(2 * d$x) + u[as.integer(d$g)] + rnorm(n, sd = exp(-0.7 + 0.4 * d$x))

fit <- statmod(y ~ s(x) + random(~ 1 | g) | sigma ~ x,
               distrib = gaussian1_distrib(), data = d)
```

The smoothing parameter of the smooth and the scale of the random effect
are hyperparameters. They are estimated by REML, whose gradient and
Hessian are computed exactly, and
[`hyper()`](https://statmodels7.github.io/statmodels7/reference/hyper.md)
reports each value with its source: the criterion that estimated it, or
`fixed` for a value held in the term.

``` r

hyper(fit)
#>   parameter           term   name  estimate  held source   id
#> 1        mu           s(x) lambda 1.3774223 FALSE   reml <NA>
#> 2        mu random(~1 | g)  sigma 0.5804077 FALSE   reml <NA>
```

[`summary()`](https://rdrr.io/r/base/summary.html) prints the
coefficients and the hyperparameters of each equation, with the
effective degrees of freedom of each penalized term.

``` r

summary(fit)
#> A statmod fit
#> 
#> Call:  statmod(formula = y ~ s(x) + random(~1 | g) | sigma ~ x, distrib = gaussian1_distrib(), 
#>             data = d)
#> 
#> Distribution: gaussian1     Observations: 400
#> 
#> === mu   [identity link]
#> 
#> Parametric terms
#>                estimate     se      z      p   lower  upper
#>   (Intercept)   0.02975 0.1332 0.2234 0.8233 -0.2313 0.2908
#> 
#> s(x)   [9 coefficients, edf 6.81]
#>                  estimate     se    z        p  lower  upper
#>   lambda [reml]     1.377 0.8096               0.4353  4.359
#>   lin              0.2184 0.0354 6.17 6.83e-10 0.1491 0.2878
#> 
#> random(~1 | g)   [20 coefficients, edf 18.45]
#>                 estimate      se  lower  upper
#>   sigma [reml]    0.5804 0.09696 0.4184 0.8052
#>   predicted effects (20 levels of g)
#>                  Min. 1st Qu.   Median Mean 3rd Qu.  Max.
#>     (Intercept)  -1.1 -0.3411 0.008096    0  0.2206 1.151
#> 
#> === sigma   [log link]
#> 
#> Parametric terms
#>                estimate      se      z       p  lower   upper
#>   (Intercept)   -0.6693 0.03661 -18.28 < 1e-16 -0.741 -0.5976
#>   x              0.3736 0.03385  11.04 < 1e-16 0.3072  0.4399
#> 
#> 95% intervals, bayesian variance
#> conditional log-likelihood -282.744059    effective df 28.26
#> cAIC 622.006    cBIC 734.801
#> fitted in 1.99 s
#> inner   max |grad|/se 1.9e-15   min eigen 0.014
#> outer   max |grad|/se 2.9e-12   min eigen 0.97
#> 1 note: print(summary(fit), notes = TRUE)
```

After the conditional log-likelihood and the information criteria, the
table ends with two lines headed `inner` and `outer`. Each gives the
largest gradient in standard errors and the smallest eigenvalue of the
curvature, for the coefficients and for the hyperparameters
respectively, and both are read at the point that the fit reached.

[`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md)
is computed at the same point and not read off the stopping rule of the
optimizer. Its `mode_error` is how far above its own penalized mode the
coefficients stopped, in log-likelihood units, and its `state` is
`converged`, `boundary`, `not converged` or `unknown`.

``` r

statmod_certificate(fit)[c("state", "mode_error")]
#> $state
#> [1] "converged"
#> 
#> $mode_error
#> [1] 4.831658e-30
```

[`predict()`](https://rdrr.io/r/stats/predict.html) returns the fitted
value of any parameter, or of a moment of the response, with a standard
error that follows from every term in the model.

``` r

predict(fit, what = "mu", newdata = head(d, n = 3), se = TRUE)
#>           fit        se      lower      upper
#> 1 -0.81384092 0.1046663 -1.0189830 -0.6086988
#> 2  0.06690214 0.1012942 -0.1316309  0.2654351
#> 3  1.41682665 0.1327520  1.1566376  1.6770157
```

## What the layer fits

Any term of `modelterms7` goes in any equation: the parametric block,
penalized smooths of one and several covariates, grouped random
intercepts and slopes, ridge, lasso, elastic net, SCAD and MCP, a
nonlinear parametric term, effects that change at estimated
break-points, score-driven dynamics and latent Markov regimes. A formula
carries at most one score-driven term, regime term or break-point term
built with `marginal = TRUE`. A score-driven term is fitted in one
system with the coefficients, and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) is built from the observed
information of that joint model. A regime term, or a marginal
break-point term, is fitted by alternating with the coefficients.

The hyperparameters of a smooth penalty are estimated by
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
or
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md),
or by
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) on
an estimate of prediction error, all passed as `outer_criterion`. A
penalty with a kink (`lasso()`, `enet()`, `scad()` and `mcp()`) has no
second derivative at the mode, so its hyperparameter is chosen along a
path of its own values, from the value that empties the block downwards.
The path is scored by `sparse_criterion`, which is
[`bic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) by
default and may also be
[`aic()`](https://statmodels7.github.io/statmodels7/reference/aic.md) or
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md). The
function
[`cv()`](https://statmodels7.github.io/statmodels7/reference/cv.md)
serves only as `sparse_criterion`.

Random effects written with the same label and the same grouping share
one covariance block, in whichever equation they appear, as with the
middle bar of `(1 | u | id)` in `brms`. The formula
`y ~ random(~ 1 | u | id) | sigma ~ random(~ 1 | u | id)` estimates the
correlation between a subject’s level and the logarithm of its scale.
The label reaches a subformula too, so
`seg(x, psi ~ random(~ 1 | u | id))` in the same formula correlates a
random break-point with a random scale. The block is one multivariate
prior over the effects of one level, parametrized by `dr_prod()` by
default, whose free coordinates are the logarithms of the standard
deviations followed by the angles of the correlation.

[`rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.md)
simulates from a model written the same way and returns the parameters
used beside the data.

``` r

sim <- rstatmod(y ~ x, distrib = gaussian1_distrib(), data = d[1:20, ],
                par = list(mu = c(1, 2), sigma = log(0.3)))
sim$par
#> $mu
#> (Intercept)           x 
#>           1           2 
#> 
#> $sigma
#> (Intercept) 
#>   -1.203973
```

The vignette [fitting a distributional
model](https://statmodels7.github.io/statmodels7/articles/fitting-a-distributional-model.html)
works through a fit from the formula to the fitted object.

## The toolkit in one package

Installing this package installs the eight members, and attaching it
attaches them. The members are listed in `Imports` and attached by the
package itself, which prints one message with the version of each;
`Depends` would attach them without a report.

``` r

statmodels7_packages()
#> [1] "basis7"         "distributions7" "linkfunctions7" "modelterms7"   
#> [5] "numericals7"    "optimizers7"    "parameters7"    "penalties7"
```

| package | what it provides |
|----|----|
| [numericals7](https://statmodels7.github.io/numericals7/) | the numerical layer at the root: stencils, batched quadrature and series, the combinatorial enumerations, special functions |
| [linkfunctions7](https://statmodels7.github.io/linkfunctions7/) | link functions as objects, with exact derivatives to fifth order in both directions |
| [distributions7](https://statmodels7.github.io/distributions7/) | univariate and multivariate distributions with exact score, information and higher derivatives |
| [optimizers7](https://statmodels7.github.io/optimizers7/) | optimization algorithms as objects, with composable stopping rules |
| [basis7](https://statmodels7.github.io/basis7/) | basis expansions with derivatives, anchored integrals and exact Gram matrices |
| [parameters7](https://statmodels7.github.io/parameters7/) | constrained parameters as maps from an unconstrained vector, exact to fourth order |
| [penalties7](https://statmodels7.github.io/penalties7/) | penalties as objects, with exact derivatives, proximal operators and the pieces that a marginal criterion consumes |
| [modelterms7](https://statmodels7.github.io/modelterms7/) | model terms as objects, each carrying its design block, its penalty and its blueprint for new data |

[`statmodels7_versions()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_versions.md)
reports the installed version of each member,
[`statmodels7_conflicts()`](https://statmodels7.github.io/statmodels7/reference/statmodels7_conflicts.md)
the names that members mask from one another, and
`statmodels7_update(action = "install")` installs or updates every
member.

## Related

The mathematical companion to the toolkit is the
[book](https://statmodels7.github.io/book/), which presents the formulas
that the packages implement, explains where each quantity comes from and
cites the source of the full argument.
