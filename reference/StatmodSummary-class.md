# A Summary of a Fitted Model

What
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
returns: the blocks of each distribution parameter, the degrees of
freedom, the information criteria and whatever has to be said about how
the numbers should be read.

## Usage

``` r
StatmodSummary(
  call = NULL,
  distrib_name = character(0),
  n_obs = integer(0),
  tables = list(),
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

  The number of observations.

- tables:

  A named list, one entry per distribution parameter, each a list of
  block records.

- edf:

  The per-term degrees of freedom.

- loglik:

  The maximized log-likelihood.

- df:

  The effective degrees of freedom in total.

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

- notes:

  Character vector of things the reader has to know.

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
