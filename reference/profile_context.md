# The Context of a Profile

Resolves the labels of `which` against the fit, rejects the coefficients
that cannot be held, and collects what every refit of the profile reads:
the specification, the stacked estimates, the positions of the held
coefficients, the settings of the inner fit and the outer criterion.

## Usage

``` r
profile_context(fit, which, hyper)
```

## Arguments

- fit:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

- which:

  The labels.

- hyper:

  `"held"` or `"refit"`.

## Value

A list.
