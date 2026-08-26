# Which Information Matrix a Fit Used

`TRUE` when the fit inverted the expected information, as
[`iwls()`](https://statmodels7.github.io/statmodels7/reference/iwls.md)
does unless asked otherwise.

## Usage

``` r
fit_expected(object)
```

## Arguments

- object:

  A
  [`StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/StatmodFit-class.md).

## Value

A single logical.

## Details

[`vcov.StatmodFit()`](https://statmodels7.github.io/statmodels7/reference/vcov.StatmodFit.md)'s
default follows this instead of choosing for itself, so a standard error
comes from the same matrix the fit did. A caller who wants the other one
asks for it.
