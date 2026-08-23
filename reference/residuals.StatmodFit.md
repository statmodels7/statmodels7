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

  A fitted model.

- type:

  The residual to compute.

- seed:

  An integer to seed the randomization with, or `NULL`. Read only where
  the distribution function jumps.

- ...:

  Unused.

## Value

A numeric vector, one entry per observation.

## Details

A residual asks whether an observation is consistent with the LAW its
row was given, and that law is one object carrying every parameter at
once. So there is one residual per observation and not one per
distribution parameter, however many of them the model develops over
covariates. What IS per parameter is a different quantity: the
contribution \\\partial \ell_i / \partial \eta\_{ip}\\ says which
equation an observation strains, and the partial residuals of one
equation are what a term's effect is drawn against.

The QUANTILE residual is \\r_i = \Phi^{-1}(F(y_i; \hat\theta_i))\\.
Under a correct model \\F(y_i; \theta_i)\\ is exactly uniform, so
\\r_i\\ is exactly standard normal – whatever the family, and whichever
of its parameters are modelled. That is why it is the default here: it
privileges no parameter, and it needs no asymptotics for its reference
distribution to hold (Dunn and Smyth, 1996).

Where the distribution function JUMPS the same construction is
randomized: \\u_i\\ is drawn uniformly on \\(F(y_i^-), F(y_i))\\ and the
residual is \\\Phi^{-1}(u_i)\\. That is exact again, at the price of
being random, so two calls give two answers. It applies to every
discrete family and, at the atom alone, to a mixed one – the
zero-adjusted wrapper of a continuous parent. `seed` makes a call
reproducible without disturbing the caller's stream; left `NULL` the
ambient state is used and nothing is set.

The PEARSON residual is \\(y_i - \mathbb{E}\[Y_i\])/\mathrm{sd}(Y_i)\\
and the RESPONSE residual the numerator alone. Both are defined against
the mean, and for a skewed family the first is not standard normal even
where the model is right, so its quantile-quantile plot misleads in
exactly the case a distributional model is for. They are offered because
they are familiar, not because they answer the question the quantile
residual does.

## References

Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals.
*Journal of Computational and Graphical Statistics* 5(3), 236–244.

## See also

[`fitted.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/fitted.StatmodFit.md),
[`predict.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/predict.StatmodFit.md)

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(80, -2, 2))
d$y <- 1 + 0.8 * d$x + rnorm(80, 0, 0.5)
fit <- statmod(y ~ x, distributions7::gaussian1_distrib(), d)
r <- residuals(fit)
c(mean = mean(r), sd = stats::sd(r))
#>         mean           sd 
#> 5.498294e-17 1.006309e+00 
```
