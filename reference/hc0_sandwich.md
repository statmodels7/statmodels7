# Sandwich Variance of the Coefficients

Select the meat of the frequentist variance \\V = V_b B V_b\\, with
\\V_b = (H + S)^{-1}\\ the inverse of the penalized information and
\\B\\ built from the scores of the observations. Pass the object as the
`sandwich` argument of
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md).

## Usage

``` r
hc0_sandwich(cluster = NULL, data = NULL, weights = c("frequency", "sampling"))

hc1_sandwich(cluster = NULL, data = NULL, weights = c("frequency", "sampling"))
```

## Arguments

- cluster:

  `NULL` (the default) or a one-sided formula with a single variable of
  the data, such as `~ id`.

- data:

  `NULL` (the default) or the data frame that holds the cluster
  variable, with one row for each observation of the fit.

- weights:

  `"frequency"` (the default) or `"sampling"`: whether the prior weights
  of the fit are frequencies or sampling weights.

## Value

A `SandwichMethod` object with the properties `type` (`"HC0"` or
`"HC1"`), `cluster`, `data` and `weights`.

## Details

The score \\u_i\\ of observation \\i\\ is the derivative of its
log-likelihood, with no penalty in it, in the whole vector of
coefficients: the fixed coefficients and the smooth coefficients of
every distribution parameter, and the random effects, which are treated
as penalized coefficients. The matrix \\B\\ is \\\sum_i u_i u_i'\\, or
\\\sum_g U_g U_g'\\ with \\U_g = \sum\_{i \in g} u_i\\ over the clusters
\\g\\ when `cluster` is given. The estimator is therefore the robust
counterpart of the frequentist variance of
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md),
conditional on the random effects and on the hyperparameters. It does
not contain the smoothing bias of a penalized term, as the frequentist
variance does not.

`hc0_sandwich()` applies no correction (HC0, and CR0 with a cluster).
`hc1_sandwich()` multiplies by \\n/(n - p)\\ (HC1), and by \\G/(G - 1)
\cdot (n - 1)/(n - p)\\ with a cluster (CR1), where \\G\\ is the number
of clusters. A cluster for each observation reduces CR0 to HC0 and CR1
to HC1. Here \\p\\ is the effective degrees of freedom of the whole
model, random effects included, which is not the number of columns of
the design where a term is penalized.

The cluster variable is read from `data` when that argument is given and
from the data of the call to
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
otherwise, which is evaluated again in the environment of the model
formula. A fit does not keep its data, so the second route fails once
the data frame is gone or has changed.

The meat always uses the score of the likelihood. With
[`ml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
the scores sum to zero at the fit, apart from the penalty. Under
[`reml()`](https://statmodels7.github.io/statmodels7/reference/reml.md)
a coordinate that the marginal criterion estimates, such as the scale of
a Gaussian, has a likelihood score that is not zero at the fit, and the
sandwich includes that offset.

A cluster must contain every group of a
[`modelterms7::random()`](https://statmodels7.github.io/modelterms7/reference/random.html)
term entirely. With no `cluster` the observations are treated as
independent given the random effects.

Models with a structural term
([`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html),
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html))
are rejected, and
[`hac_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hac_sandwich.md)
covers them. A term whose block moves with its coefficients
([`modelterms7::seg()`](https://statmodels7.github.io/modelterms7/reference/seg.html),
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html)
and the smoothed
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
and
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html))
is accepted: its block is the Jacobian of its contribution in the
coefficients and the position of the break-point, read at the fit, so
the score of the position is part of \\u_i\\. The position of a sharp
[`modelterms7::jump()`](https://statmodels7.github.io/modelterms7/reference/jump.html)
or
[`modelterms7::jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.html)
is held by the fit, and its slot has no variance, as in the other
conventions; a block that is still a frozen working linearization is
rejected, because a step has no derivative in its position. A coordinate
that a kinked penalty holds at zero has no variance either.

The prior weights \\w_i\\ of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
multiply the log-likelihood, and the argument `weights` says how the
sandwich reads them. With `weights = "frequency"` (the default) row
\\i\\ stands for \\w_i\\ identical observations. Without a cluster the
copies are separate observations and \\B = \sum_i w_i g_i g_i'\\, with
\\g_i\\ the score of one copy, and \\n = \sum_i w_i\\ in the HC1 factor.
With a cluster the copies stay in the cluster of the row, so \\U_g =
\sum\_{i \in g} w_i g_i\\, \\G\\ is the number of clusters with positive
weight and \\n = \sum_i w_i\\. A row with weight zero contributes
nothing and is not counted. The result equals the unweighted estimator
on the data expanded to \\w_i\\ copies of each row, for integer weights,
and it changes when every weight is multiplied by a constant.

With `weights = "sampling"` the weights are sampling weights, such as
the inverse probabilities of selection of a survey, and each row is one
unit whose score is the weighted score \\w_i g_i\\. Without a cluster
\\B = \sum_i w_i^2 g_i g_i'\\ and \\n\\ in the HC1 factor is the number
of rows with positive weight. With a cluster \\B\\ is the same as under
the frequency reading, and only the \\n\\ of the CR1 factor differs. The
variance does not change when every weight is multiplied by a constant.
The point estimates are the same under both readings, and so are the
variances of the other conventions, which read the weights as
frequencies.

## References

White, H. (1980). A heteroskedasticity-consistent covariance matrix
estimator and a direct test for heteroskedasticity. *Econometrica*,
48(4), 817–838.

MacKinnon, J. G. and White, H. (1985). Some
heteroskedasticity-consistent covariance matrix estimators with improved
finite sample properties. *Journal of Econometrics*, 29(3), 305–325.

Liang, K.-Y. and Zeger, S. L. (1986). Longitudinal data analysis using
generalized linear models. *Biometrika*, 73(1), 13–22.

## See also

[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md),
which takes the object.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(150))
dd$y <- 1 + 2 * dd$x + rnorm(150, sd = 0.2 + dd$x)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
sqrt(diag(vcov(fit, vcov_type = frequentist_vcov())))
#>    mu:(Intercept)              mu:x sigma:(Intercept) 
#>        0.13312662        0.23015913        0.05812382 
sqrt(diag(vcov(fit, vcov_type = frequentist_vcov(
  sandwich = hc0_sandwich()))))
#>    mu:(Intercept)              mu:x sigma:(Intercept) 
#>        0.09366361        0.23964065        0.06933911 
sqrt(diag(vcov(fit, vcov_type = frequentist_vcov(
  sandwich = hc1_sandwich()))))
#>    mu:(Intercept)              mu:x sigma:(Intercept) 
#>        0.09461453        0.24207361        0.07004307 

# observations in clusters
dd$id <- rep(1:30, each = 5)
fit2 <- statmod(y ~ x, distributions7::gaussian1_distrib(), dd)
vcov(fit2, vcov_type = frequentist_vcov(
  sandwich = hc1_sandwich(cluster = ~ id)), param = "mu")
#>                mu:(Intercept)        mu:x
#> mu:(Intercept)     0.01088101 -0.02516280
#> mu:x              -0.02516280  0.07551173
```
