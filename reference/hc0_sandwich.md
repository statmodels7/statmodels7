# Sandwich Variance of the Coefficients

Select the meat of the frequentist variance \\V = V_b B V_b\\, with
\\V_b = (H + S)^{-1}\\ the inverse of the penalized information and
\\B\\ built from the scores of the observations. Pass the object as the
`sandwich` argument of
[`frequentist_vcov()`](https://statmodels7.github.io/statmodels7/reference/frequentist_vcov.md).

## Usage

``` r
hc0_sandwich(cluster = NULL, data = NULL)

hc1_sandwich(cluster = NULL, data = NULL)
```

## Arguments

- cluster:

  `NULL` (the default) or a one-sided formula with a single variable of
  the data, such as `~ id`.

- data:

  `NULL` (the default) or the data frame that holds the cluster
  variable, with one row for each observation of the fit.

## Value

A `SandwichMethod` object with the properties `type` (`"HC0"` or
`"HC1"`), `cluster` and `data`.

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
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)),
with a term whose block moves with its coefficients (the break-point
terms and
[`modelterms7::nl()`](https://statmodels7.github.io/modelterms7/reference/nl.html))
are rejected.
[`hac_sandwich()`](https://statmodels7.github.io/statmodels7/reference/hac_sandwich.md)
covers the structural terms.

The prior weights \\w_i\\ of
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
are frequency weights: row \\i\\ stands for \\w_i\\ identical
observations. Without a cluster the copies are separate observations and
\\B = \sum_i w_i g_i g_i'\\, with \\g_i\\ the score of one copy, and \\n
= \sum_i w_i\\ in the HC1 factor. With a cluster the copies stay in the
cluster of the row, so \\U_g = \sum\_{i \in g} w_i g_i\\, \\G\\ is the
number of clusters with positive weight and \\n = \sum_i w_i\\. A row
with weight zero contributes nothing and is not counted. The result
equals the unweighted estimator on the data expanded to \\w_i\\ copies
of each row, for integer weights.

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
