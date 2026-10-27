# The Information at the Context's Point

The Information at the Context's Point

## Usage

``` r
ctx_information(ctx, spec, design, coef, hyper, expected, emethod = NULL)
```

## Arguments

- ctx:

  A context, or `NULL`.

- spec, design, coef:

  The fallback arguments, used when `ctx` is `NULL`.

- hyper:

  The hyperparameters, for the context's own check.

- expected:

  Whether the expected information is read.

- emethod:

  The method that obtains the expected information: `NULL` (the closed
  form of the family where it has one, and
  [`distributions7::expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.html)
  otherwise) or an object that inherits from
  [distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html).
  It is read only where the expected information is used.

## Value

A square matrix.
