# The Method With Which a Context Was Built

The Method With Which a Context Was Built

## Usage

``` r
ctx_emethod(ctx)
```

## Arguments

- ctx:

  A context, or `NULL`.

## Value

`NULL`, for the closed form of the family or the outer product of the
scores, or the object that inherits from
[distributions7::ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.html)
with which the context was built. It is `NULL` where `ctx` is `NULL` or
carries none.
