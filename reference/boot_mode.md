# The Kind of Model That the Bootstrap Treats Apart

Returns `"structural"` for a model with a
[`modelterms7::gas()`](https://statmodels7.github.io/modelterms7/reference/gas.html)
or
[`modelterms7::regime()`](https://statmodels7.github.io/modelterms7/reference/regime.html)
term, `"refresh"` for a model with a term whose block moves with its
coefficients, and `NULL` for every other model.

## Usage

``` r
boot_mode(design0)
```

## Arguments

- design0:

  The design of the fit.

## Value

A string or `NULL`.
