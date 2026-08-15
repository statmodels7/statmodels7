# Refuse a Context That Belongs Somewhere Else

Stops when a context is handed to a consumer reading a different point.

## Usage

``` r
ctx_usable(ctx, coef, hyper)
```

## Arguments

- ctx:

  A context, or `NULL`.

- coef:

  The coefficients the caller is reading.

- hyper:

  The hyperparameters the caller is reading.

## Value

`TRUE` when the context is usable, `FALSE` when it is `NULL`; an error
when it belongs to another point.

## Details

A stale cached information is a silently wrong gradient, which is the
one failure a shared cache can introduce and the one that would be
hardest to find. The comparison is over the coefficients and the
hyperparameters, which is linear in their length and so negligible
beside the \\O(np^2)\\ assembly it guards.
