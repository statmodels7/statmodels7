# The Parameters of a Regime Term

Returns the levels of the regimes, the transition matrix and the
stationary distribution at the fitted parameters of the term.

## Usage

``` r
regime_parameters(r)
```

## Arguments

- r:

  An entry of the `regimes` element of
  [`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md).

## Value

A list with `mu` (the levels), `P` (the transition matrix, whose rows
are the current regime) and `delta` (the stationary distribution).
