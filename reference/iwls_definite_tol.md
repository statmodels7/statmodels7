# The Relative Level Under Which a Penalized Information Is Indefinite

Returns the relative floor with which
[`iwls_fit()`](https://statmodels7.github.io/statmodels7/reference/iwls_fit.md)
tests whether the penalized information of the step is positive
definite: the `repair_floor` of the
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
that was settled, and the `floor` of the solver when the curvature was
named.

## Usage

``` r
iwls_definite_tol(method)
```

## Arguments

- method:

  An
  [`Iwls()`](https://statmodels7.github.io/statmodels7/reference/Iwls-class.md)
  object.

## Value

A single positive number.

## See also

[`pieces_definite()`](https://statmodels7.github.io/statmodels7/reference/pieces_definite.md),
[`auto_info()`](https://statmodels7.github.io/statmodels7/reference/auto_info.md)
