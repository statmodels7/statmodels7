# Which Information Matrix a Fit Used

`TRUE` when the fit inverted the expected information, which is what
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
does unless asked otherwise.

## Usage

``` r
fit_expected(object)
```

## Arguments

- object:

  A
  [`StatmodFit`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

## Value

A single logical.

## Details

The default of
[`vcov.StatmodFit`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)
follows this rather than choosing for itself, so that a standard error
comes from the same matrix the fit did, and a caller who wants the other
one asks for it.
