# The Specification of a Model, Before It Is Fitted

Everything the formula, the data and the distribution produce: one
equation per distribution parameter, the terms each equation names built
against the data, the response, and the prior weights and offsets.

## Usage

``` r
StatmodSpec(
  formula = NULL,
  distrib = NULL,
  equations = list(),
  terms = list(),
  response = NULL,
  n_obs = integer(0),
  weights = integer(0),
  offsets = list(),
  intercepts = logical(0),
  newdata = NULL,
  structural = list(),
  linpar = list(),
  threads = 1L,
  workers = 1L
)
```

## Arguments

- formula:

  The model formula, as given.

- distrib:

  The distributions7 object.

- equations:

  A named list of one-sided formulas, one per parameter, in the family's
  order.

- terms:

  A named list, one entry per parameter, each a named list of built
  modelterms7 terms.

- response:

  The evaluated left-hand side.

- n_obs:

  The number of observations.

- weights:

  Prior weights, one per observation.

- offsets:

  A named list of offsets, one per parameter or `NULL`.

- intercepts:

  A named logical, whether each equation carried one.

## Value

An object of class `StatmodSpec`.

## See also

[`statmod_spec`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md)

## Examples

``` r
dd <- data.frame(y = rnorm(10), x = runif(10))
S7::S7_inherits(statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd),
                StatmodSpec)
#> [1] TRUE
```
