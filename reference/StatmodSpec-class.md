# The Specification of a Model, Before It Is Fitted

The specification holds everything that the formula, the data and the
distribution produce before any fitting: one equation per distribution
parameter, the terms that each equation names built against the data,
the response, and the prior weights and offsets. A fit keeps one of
these in its `spec` property, so it is also what
[`summary()`](https://rdrr.io/r/base/summary.html),
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
  held_coef = list(),
  held_hyper = list(),
  threads = 1L,
  cores = 1L
)
```

## Arguments

- formula:

  The model formula as given, including the bars.

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

- newdata:

  The rows that a respecification reads, or `NULL` for the fitting rows.

- structural:

  A named list of the parameters of the structural term, keyed by the
  term; empty until a fit sets them.

- linpar:

  A list, the options of the implicit parametric block as
  [`linpar_options()`](https://statmodels7.github.io/statmodels7/reference/linpar_options.md)
  returns them.

- held_coef:

  A named list, one entry per distribution parameter, of coefficients
  held at a value, each a named numeric vector keyed by coefficient
  name. The list is empty by default.

- held_hyper:

  A named list keyed by distribution parameter, then by the key of a
  term, of hyperparameters held at a value, each a named numeric vector
  on the scale of the hyperparameter. The list is empty by default.

- threads:

  The number of threads, a single integer.

- cores:

  The number of separate R processes for the independent fits of a
  cross-validation and of a kinked path's product grid, a single
  integer. It is distinct from `threads`, the number of threads inside
  one process.

## Value

An object of class `StatmodSpec` with one property per argument above
and the properties described next. `linpar` holds the options of the
implicit parametric block and is set by
[`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md).
`newdata` holds the rows that a respecification reads, `NULL` until
[`statmod_respec()`](https://statmodels7.github.io/statmodels7/reference/statmod_respec.md)
sets it. `structural` holds the parameters of a structural term where a
fit has set them. `held_coef` holds coefficients held at a value, by
name, empty by default. `held_hyper` holds hyperparameters at a value,
by term and name, empty by default. `threads` and `cores` hold the
thread and process counts that
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
sets.

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

# The specification has one equation per parameter of the family, whatever
# the formula wrote.
names(spec@equations)
#> [1] "mu"    "sigma"
spec@n_obs
#> [1] 10
```
