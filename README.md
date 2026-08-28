
<!-- README.md is generated from README.Rmd. Please edit that file, then
     regenerate with devtools::build_readme(). Do not use knitr::knit(): it
     processes the code but leaves this YAML header in the output as literal
     text, which GitHub and pkgdown both render verbatim. -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/statmodels7/statmodels7/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/statmodels7/statmodels7/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/statmodels7/statmodels7/graph/badge.svg)](https://app.codecov.io/gh/statmodels7/statmodels7)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# statmodels7 <img src="man/figures/logo.png" align="right" height="139" alt="" />

`{statmodels7}` is the modeling layer of the
[statmodels7](https://statmodels7.github.io) toolkit, and the package
that installs the rest of it. `statmod()` reads one formula carrying an
equation for every parameter of a distribution and fits it, with any
term of [modelterms7](https://statmodels7.github.io/modelterms7/) in any
of those equations.

## Installation

``` r
# install.packages("pak")
pak::pak("statmodels7/statmodels7")
```

Attaching it attaches the eight member packages.

``` r
library(statmodels7)
#> -- Attaching the statmodels7 toolkit 0.90.0
#> v basis7         0.6.0    v numericals7    0.12.0
#> v distributions7 0.41.0   v optimizers7    0.8.0
#> v linkfunctions7 0.3.0    v parameters7    0.18.0
#> v modelterms7    0.64.0   v penalties7     0.19.0
```

## Fitting a model

The equations are separated by a bar, one for each parameter of the
distribution, and a parameter with no equation of its own is given an
intercept. Below, the mean carries a penalized smooth of `x` and a
random intercept over `g`, while the standard deviation is linear in `x`
on its log scale.

``` r
set.seed(1)
n <- 400
d <- data.frame(x = runif(n, -2, 2), g = factor(rep(1:20, length.out = n)))
u <- rnorm(20, 0, 0.5)
d$y <- sin(2 * d$x) + u[as.integer(d$g)] + rnorm(n, sd = exp(-0.7 + 0.4 * d$x))

fit <- statmod(y ~ s(x) + random(~ 1 | g) | sigma ~ x,
               gaussian1_distrib(), d)
```

The smoothing parameter of the smooth and the scale of the random effect
are hyperparameters rather than coefficients. They are estimated by
REML, whose gradient and Hessian are computed exactly, and `hyper()`
reports each value with what put it there.

``` r
hyper(fit)
#>   parameter           term   name  estimate  held source
#> 1        mu           s(x) lambda 1.3902202 FALSE   reml
#> 2        mu random(~1 | g)  sigma 0.5796875 FALSE   reml
```

`summary()` prints the two together, one block per equation, with the
effective degrees of freedom each penalized term spends.

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
#>   (Intercept)   0.02981 0.1329 0.2243 0.8225 -0.2307 0.2903
#> 
#> s(x)   [9 coefficients, edf 6.85]
#>                  estimate      se     z         p  lower  upper
#>   lambda [reml]    1.3900 0.81460                 0.4409 4.3840
#>   lin              0.2187 0.03498 6.251 4.076e-10 0.1501 0.2872
#> 
#> random(~1 | g)   [20 coefficients, edf 18.50]
#>                 estimate      se  lower  upper
#>   sigma [reml]    0.5797 0.09684 0.4178 0.8043
#> 
#> === sigma   [log link]
#> 
#> Parametric terms
#>                estimate      se      z       p   lower   upper
#>   (Intercept)   -0.7032 0.03537 -19.89 < 1e-16 -0.7726 -0.6339
#>   x              0.3903 0.03181  12.27 < 1e-16  0.3280  0.4527
#> 
#> 95% intervals, bayesian variance
#> conditional log-likelihood -282.115367    effective df 28.35
#> cAIC 620.936    cBIC 734.104
#> fitted in 2.55 s   search: converged
#> certificate: CONVERGED   outer gradient 6.79e-07   2.19e-11 above the mode
#> 1 note: print(summary(fit), notes = TRUE)
```

The last line of that table comes from `statmod_certificate()`, which is
computed at the point the fit reached rather than read off the
optimizer’s stopping rule: it reports how far above its own mode the
coefficients stopped and whether the hyperparameters are stationary
there.

``` r
statmod_certificate(fit)[c("state", "gradient", "mode_error")]
#> $state
#> [1] "converged"
#> 
#> $gradient
#> [1] 6.791092e-07
#> 
#> $mode_error
#> [1] 2.192239e-11
```

`predict()` answers for any parameter, or for a moment of the response,
and carries a standard error through whatever the terms are.

``` r
predict(fit, "mu", head(d, 3), se = TRUE)
#>           fit         se      lower      upper
#> 1 -0.81662977 0.09999405 -1.0126145 -0.6206450
#> 2  0.06755757 0.09687091 -0.1223059  0.2574211
#> 3  1.41839312 0.12838002  1.1667729  1.6700133
```

## What the layer fits

Every term `modelterms7` defines goes in any equation: the parametric
block, penalized smooths of one and several covariates, grouped random
intercepts and slopes, ridge, lasso, elastic net, SCAD and MCP, a
nonlinear parametric term, effects that change at estimated
break-points, score-driven dynamics and latent Markov regimes. A term
that rewrites the likelihood rather than contributing a design block is
fitted in the same system as the coefficients, with the observed
information of the joint model behind `vcov()`.

The hyperparameters of a smooth penalty come from `reml()` or `ml()`, or
from `aic()`, `bic()` or `cv()` on an estimate of prediction error. A
penalty with a kink has no second derivative at the mode, so its
hyperparameter is chosen instead along a path of its own values, from
the value that empties the block downwards.

`rstatmod()` simulates from a model written the same way, which is what
a simulation study needs: the truth is returned beside the data.

``` r
sim <- rstatmod(y ~ x, gaussian1_distrib(), d[1:20, ],
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
attaches them. They sit in `Imports` rather than `Depends`, which is
what lets the attaching be reported rather than merely happen: `Depends`
would attach them in the order the field lists them, with no message and
no way for a caller to see what arrived at which version.

``` r
statmodels7_packages()
#> [1] "basis7"         "distributions7" "linkfunctions7" "modelterms7"   
#> [5] "numericals7"    "optimizers7"    "parameters7"    "penalties7"
```

| package | what it provides |
|----|----|
| [numericals7](https://statmodels7.github.io/numericals7/) | the numerical layer at the root: stencils, batched quadrature and series, the combinatorial enumerations, special functions |
| [linkfunctions7](https://statmodels7.github.io/linkfunctions7/) | link functions as objects, with exact derivatives to fourth order in both directions |
| [distributions7](https://statmodels7.github.io/distributions7/) | univariate and multivariate distributions with exact score, information and higher derivatives |
| [optimizers7](https://statmodels7.github.io/optimizers7/) | optimization algorithms as objects, with composable stopping rules |
| [basis7](https://statmodels7.github.io/basis7/) | basis expansions with derivatives, anchored integrals and exact Gram matrices |
| [parameters7](https://statmodels7.github.io/parameters7/) | constrained parameters as maps from an unconstrained vector, exact to fourth order |
| [penalties7](https://statmodels7.github.io/penalties7/) | penalties as objects, with their normalizing constants and the pieces a marginal criterion consumes |
| [modelterms7](https://statmodels7.github.io/modelterms7/) | model terms as objects, each carrying its design block, its penalty and its blueprint for new data |

`statmodels7_versions()` reports what is installed,
`statmodels7_conflicts()` what masks what, and
`statmodels7_update("install")` installs or updates every member.

## Related

The mathematical companion to the whole toolkit is the
[book](https://statmodels7.github.io/book/), which gives every formula
the packages implement with its derivation and its citation.
