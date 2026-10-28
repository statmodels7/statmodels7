# The Constant at Which a Statically Started Filter Begins

Sets the property `init_value` of every
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
term written with `init = "static"` to the predictor of the model
without the term, read at the first row of each group of the filter.

## Usage

``` r
static_filter_starts(spec, split, data, weights, offsets, distrib, env, linpar)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- split:

  The equations, as
  [`statmod_equations()`](https://statmodels7.github.io/statmodels7/reference/statmod_equations.md)
  returns them.

- data, weights, offsets, distrib, env, linpar:

  The arguments of
  [`statmod_spec()`](https://statmodels7.github.io/statmodels7/reference/statmod_spec.md),
  with `env` the formula's environment.

## Value

The specification, with the constant set on each such term.

## Details

With `init = "static"` the first predictor of every group of the filter
equals this constant whatever the parameters are, so the constant of the
equation stays the level of the filter and does not double as its
starting value. The model without the term is the same formula with each
`gas()` call replaced by the development of its level (the right side of
an `omega ~` subformula) or, without one, by an intercept, so that the
static model carries a level where the filter carries one. It is fitted
once by
[`statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.md)
with the same family, data, weights and offsets, and the constant of
each group is its predictor of the parameter at the group's first row.
With covariates the constant follows them: on the share of front-seat
casualties in Seatbelts, `share ~ law + gas(...)` starts at the
predictor of `share ~ law` at the first month, where the intercept-only
fit would start it at the mean of the whole series and lose 2.07 in
log-likelihood.

The fitted static model travels with the term as `init_source`, so a
prediction that rebuilds the series on other rows reads the constants of
its groups from the same model
([`static_filter_restart()`](https://statmodels7.github.io/statmodels7/reference/static_filter_restart.md)).
