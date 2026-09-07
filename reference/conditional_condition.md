# The Condition an Unavailable Correction Warns Through

A classed warning, so a caller can catch the one raised when
`vcov(type = "unconditional")` cannot read the hyperparameters' own
uncertainty and returns the conditional variance instead.

## Usage

``` r
conditional_condition(msg)
```

## Arguments

- msg:

  The message.

## Value

A condition of class `statmod_conditional_variance`.

## Details

It carries a class of its own for the same reason the other two do:
[`summary.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/summary.StatmodFit.md)
calls
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
more than once, so without one the same message would reach the reader
several times, and muffling it by position would swallow whatever else
was raised on the way.

## See also

[`hyper_correction()`](https://statmodels7.github.io/statmodels7/reference/hyper_correction.md),
which decides whether it fires, and
[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md),
which raises it.
