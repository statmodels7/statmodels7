# The Entry of the Evaluated Model for a Regime Term

Evaluates the predictors of the fit and returns the entry of the regime
term `name`, with the static predictor, the callbacks of the density and
the parameters of the term.

## Usage

``` r
regime_entry(fit, name)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- name:

  The name of the term.

## Value

A list as in the `regimes` element of
[`statmod_eta()`](https://statmodels7.github.io/statmodels7/reference/statmod_eta.md).
