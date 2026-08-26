# The Specification of a Model, Before It Is Fitted

Everything the formula, the data and the distribution produce, before
any fitting: one equation per distribution parameter, the terms each
equation names built against the data, the response, and the prior
weights and offsets. A fit keeps one of these in its `spec` property, so
it is also what [`summary()`](https://rdrr.io/r/base/summary.html),
[`predict()`](https://rdrr.io/r/stats/predict.html) and the accessors
read.

Build one with
[`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md),
which interprets the formula and validates. This raw constructor takes
the pieces already made.

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

  The model formula, as given, bars included.

- distrib:

  The distributions7 distribution object.

- equations:

  A named list of one-sided formulas, one per parameter, in the family's
  order.

- terms:

  A named list with one entry per parameter, each a named list of built
  modelterms7 terms keyed by the term's call as written.

- response:

  The evaluated left-hand side.

- n_obs:

  The number of observations, a single integer.

- weights:

  Prior weights, a numeric vector of length `n_obs`.

- offsets:

  A named list of offsets, one entry per parameter, `NULL` where an
  equation has none.

- intercepts:

  A named logical, one per parameter: whether that equation's parametric
  block carried an intercept.

## Value

An object of class `StatmodSpec` with one property per argument above,
plus `threads`, `workers`, `linpar`, `newdata` and `structural`, which
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
fills in.

## See also

[`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md),
the constructor to use,
[`statmod_design()`](https://statmodels7.github.io/statmodels7/reference/statmod_design.md)
for what is assembled from one.

## Examples

``` r
dd <- data.frame(y = rnorm(10), x = runif(10))
spec <- statmod_spec(y ~ x, distributions7::gaussian1_distrib(), dd)
S7::S7_inherits(spec, StatmodSpec)
#> [1] TRUE

# One equation per parameter of the family, whatever the formula wrote.
names(spec@equations)
#> [1] "mu"    "sigma"
spec@n_obs
#> [1] 10
```
