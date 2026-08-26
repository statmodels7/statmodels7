# The Penalty's Hessian at the Context's Point

The Penalty's Hessian at the Context's Point

## Usage

``` r
ctx_penalty(ctx, spec, design, coef, hyper)
```

## Arguments

- ctx:

  A context, or `NULL`.

- spec, design, coef, hyper:

  The fallback arguments.

## Value

A square matrix.

## Details

The non-finite entries are zeroed here, once. Three callers did it
separately before, and a fourth would have had to remember.
