# Check That a Context Belongs to the Point Read

Signals an error when a context is handed to a consumer that reads a
different point.

## Usage

``` r
ctx_usable(ctx, coef, hyper)
```

## Arguments

- ctx:

  A context, or `NULL`.

- coef:

  The coefficients that the caller is reading.

- hyper:

  The hyperparameters that the caller is reading.

## Value

`TRUE` when the context is usable, `FALSE` when it is `NULL`; an error
when it belongs to another point.

## Details

A cached information that belongs to another point yields a wrong
gradient without any other sign of failure. The comparison is over the
coefficients and the hyperparameters, so its cost is linear in their
length, which is negligible beside the \\O(np^2)\\ assembly that it
guards.
