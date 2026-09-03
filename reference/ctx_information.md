# The Information at the Context's Point

The Information at the Context's Point

## Usage

``` r
ctx_information(ctx, spec, design, coef, hyper, expected, approx = "opg")
```

## Arguments

- ctx:

  A context, or `NULL`.

- spec, design, coef:

  The fallback arguments, used when `ctx` is `NULL`.

- hyper:

  The hyperparameters, for the context's own check.

- expected:

  Whether the expected information is wanted.

- approx:

  The approximation for the expected information.

## Value

A square matrix.
