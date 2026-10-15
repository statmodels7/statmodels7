# A Summary of a Fitted Model

What
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
returns: the blocks of each distribution parameter, the degrees of
freedom, the information criteria and the notes that qualify the
numbers.

## Usage

``` r
StatmodSummary(
  call = NULL,
  distrib_name = character(0),
  n_obs = integer(0),
  weight_total = integer(0),
  tables = list(),
  classes = list(),
  links = character(0),
  edf = NULL,
  structural = NULL,
  loglik = integer(0),
  df = integer(0),
  aic = integer(0),
  bic = integer(0),
  converged = logical(0),
  elapsed = integer(0),
  level = integer(0),
  type = character(0),
  test = "wald",
  max_coef = NULL,
  notes = character(0),
  certificate = NULL
)
```

## Arguments

- call:

  The fit's call.

- distrib_name:

  The distribution's name.

- n_obs:

  The number of rows of the data.

- weight_total:

  The sum of the prior weights, which
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) reports.

- tables:

  A named list, one entry per distribution parameter, each a list of
  block records.

- classes:

  The covariance blocks that two or more pieces share (terms, or
  developments of one term), one block record each. They belong to no
  equation and are printed ahead of the equations.

- links:

  The link of each distribution parameter, a character vector named by
  parameter.

- edf:

  The per-term degrees of freedom.

- structural:

  The table of the structural term's parameters, or `NULL`.

- loglik:

  The conditional log-likelihood at the fit.

- df:

  The effective degrees of freedom in total, including the correction
  for estimated hyperparameters where it was requested.

- aic, bic:

  The information criteria.

- converged:

  Whether every loop stopped on its own rule.

- elapsed:

  The elapsed time, in seconds.

- level:

  The confidence level the intervals were built at.

- type:

  Which variance convention was used.

- test:

  Which statistic the coefficient tables report: `"wald"`, `"lr"`,
  `"score"` or `"gradient"`.

- max_coef:

  The number of coefficient rows that each block prints, or `NULL` where
  the option `statmodels7.summary_max_coef` decides.

- notes:

  A character vector of the notes that qualify the numbers.

- certificate:

  The list returned by
  [`statmod_certificate()`](https://statmodels7.github.io/statmodels7/reference/statmod_certificate.md),
  or `NULL` where it could not be computed.

## Value

An object of class `StatmodSummary`.

## See also

[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)

## Examples

``` r
dd <- data.frame(y = rnorm(30), x = runif(30))
S7::S7_inherits(summary(statmod(y ~ x,
                                distributions7::gaussian1_distrib(), dd)),
                StatmodSummary)
#> [1] TRUE
```
