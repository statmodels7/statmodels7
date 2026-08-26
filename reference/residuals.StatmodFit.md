# The Residuals of a Fitted Model

One residual per observation, comparing it with the whole distribution
the model puts on it rather than with any one of that distribution's
parameters.

## Usage

``` r
# S3 method for class 'StatmodFit'
residuals(
  object,
  type = c("quantile", "pearson", "response"),
  seed = NULL,
  ...
)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- type:

  `"quantile"` (the default), `"pearson"` or `"response"`. Matched with
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

- seed:

  An integer to seed the randomization with, or `NULL` for the ambient
  state. Read only where the distribution function jumps, so it changes
  nothing for a continuous family.

- ...:

  Unused.

## Value

A numeric vector with one entry per observation. Standard normal under a
correct model for `"quantile"`; approximately standardized for
`"pearson"`; on the response's own scale for `"response"`.

## One residual per observation

A residual asks whether an observation is consistent with the law its
row was given, and that law is one object carrying every parameter at
once. So there is one residual per observation, however many parameters
the model develops over covariates.

A per-parameter quantity exists and is a different thing: the
contribution \\\partial \ell_i / \partial \eta\_{ip}\\ says which
equation an observation strains, and the partial residuals of one
equation are what a term's effect is drawn against.

## The quantile residual, and why it is the default

\$\$r_i = \Phi^{-1}(F(y_i; \hat\theta_i))\$\$

Under a correct model \\F(y_i; \theta_i)\\ is exactly uniform, so
\\r_i\\ is exactly standard normal: whatever the family, and whichever
of its parameters are modeled. It privileges no parameter, and its
reference distribution needs no asymptotics (Dunn and Smyth, 1996).

## Where the distribution function jumps

The construction is randomized: \\u_i\\ is drawn uniformly on
\\(F(y_i^-), F(y_i))\\ and the residual is \\\Phi^{-1}(u_i)\\. Exact
again, at the price of being random, so two calls give two answers.

This applies to every discrete family, and at the atom alone to a mixed
one, which is the zero-adjusted wrapper of a continuous parent. `seed`
makes a call reproducible without disturbing the caller's stream; left
`NULL`, the ambient state is used and nothing is set.

## Pearson and response residuals

\\(y_i - \mathbb{E}\[Y_i\]) / \mathrm{sd}(Y_i)\\ and its numerator
alone. Both are defined against the mean, and for a skewed family the
Pearson residual is not standard normal even where the model is right,
so its quantile-quantile plot misleads in exactly the case a
distributional model is for. They are here because they are familiar.

## References

Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and Graphical Statistics* 5(3), 236–244.

## See also

[`fitted.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/fitted.StatmodFit.md)
and
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
for the fitted parameters these are read against.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(200, -2, 2))
d$y <- 1 + 0.8 * d$x + rnorm(200, 0, 0.5)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), d)

# Under a correct model the quantile residuals are standard normal.
r <- residuals(fit)
c(mean = mean(r), sd = stats::sd(r))
#>         mean           sd 
#> 2.090233e-17 1.002509e+00 
stats::shapiro.test(r)$p.value
#> [1] 0.8361219

# For a discrete family the same residual is randomized, so two calls
# differ unless a seed is given.
dp <- data.frame(x = runif(200, -1, 1))
dp$y <- rpois(200, exp(1 + dp$x))
fp <- statmod(y ~ x, distributions7::poisson_distrib(), dp)
identical(residuals(fp, seed = 1), residuals(fp, seed = 1))
#> [1] TRUE
identical(residuals(fp, seed = 1), residuals(fp, seed = 2))
#> [1] FALSE
```
