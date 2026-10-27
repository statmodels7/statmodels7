# The Residuals of a Fitted Model

Returns one residual per observation, which compares the observation
with the whole distribution that the model assigns to it and not with
any single parameter of that distribution.

## Usage

``` r
# S3 method for class 'StatmodFit'
residuals(
  object,
  residual_type = c("quantile", "pearson", "response"),
  seed = NULL,
  ...
)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- residual_type:

  `"quantile"` (the default), `"pearson"` or `"response"`. The argument
  is matched with
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

- seed:

  An integer used to seed the randomization, or `NULL` for the ambient
  state. Read only where the distribution function jumps, so it changes
  nothing for a continuous family.

- ...:

  Unused.

## Value

A numeric vector with one entry per row of the data. The residuals are
standard normal under a correct model for `"quantile"`; approximately
standardized for `"pearson"`; on the response's own scale for
`"response"`.

## One residual per observation

A residual measures whether an observation is consistent with the law
that was given to its row, and that law is one object that carries every
parameter at once. There is therefore one residual per observation,
however many parameters that the model develops over covariates.

A per-parameter quantity is a different object: the contribution
\\\partial \ell_i / \partial \eta\_{ip}\\ identifies which equation an
observation strains, and the partial residuals of one equation are the
quantities against which a term's effect is drawn.

## The quantile residual, and why it is the default

\$\$r_i = \Phi^{-1}(F(y_i; \hat\theta_i))\$\$

Under a correct model \\F(y_i; \theta_i)\\ is uniform where the
distribution function is continuous, so \\r_i\\ is standard normal, for
any such family and whichever of its parameters are modeled. The
residual does not single out any parameter, and its reference
distribution needs no asymptotics (Dunn and Smyth, 1996). Where the
distribution function jumps, the residual is randomized as described
next.

## Where the distribution function jumps

The construction is randomized: \\u_i\\ is drawn uniformly on
\\(F(y_i^-), F(y_i))\\ and the residual is \\\Phi^{-1}(u_i)\\. The
residual is again exactly standard normal under a correct model, but it
is random, so that two calls give two answers.

Where \\F(y_i)\\ exceeds one half, the residual is read from the
survival function, \\\Phi^{-1}(u_i) = -\Phi^{-1}(1 - u_i)\\, and \\1 -
u_i\\ is computed as \\1 - F\\ directly, so a count far in the upper
tail, where \\F\\ rounds to one, has a finite residual. A family whose
distribution function does not take `lower.tail` keeps the lower-tail
reading.

The randomization applies to every discrete family and, for a mixed
family such as the zero-adjusted wrapper of a continuous parent, only at
the declared atoms. `seed` makes a call reproducible without disturbing
the stream of the caller. If `seed` is `NULL`, the ambient state is used
and the seed is not set.

## Pearson and response residuals

The Pearson residual is \\(y_i - \mathbb{E}\[Y_i\]) /
\mathrm{sd}(Y_i)\\, and the response residual is its numerator. Both are
defined against the mean, and for a skewed family the Pearson residual
is not standard normal even where the model is right, so its
quantile-quantile plot can mislead. They are included for comparison
with the residuals that other modeling packages report.

## References

Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and Graphical Statistics* 5(3), 236–244.

## See also

[`fitted.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/fitted.StatmodFit.md)
and
[`predict.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)
for the fitted parameters against which these are read.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(200, -2, 2))
d$y <- 1 + 0.8 * d$x + rnorm(200, 0, 0.5)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), d)

# Under a correct model the quantile residuals are standard normal.
r <- residuals(fit)
c(mean = mean(r), sd = stats::sd(r))
#>          mean            sd 
#> -1.949639e-17  9.974843e-01 
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
