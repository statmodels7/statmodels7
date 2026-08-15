# The Penalized Matrix and Its Inverse

\\K = H + S\\ with the OBSERVED information, its Cholesky factor and its
inverse, which the gradient and the Hessian both read.

## Usage

``` r
ctx_penalized(ctx, spec, design, coef, hyper)
```

## Arguments

- ctx:

  A context, or `NULL`.

- spec, design, coef, hyper:

  The fallback arguments.

## Value

A list with `K` and `inv`, or `NULL`.

## Details

The inverse of a sparse matrix is dense, so densifying costs nothing
sparsity could have kept. `NULL` is returned where the factorization
fails, which is the answer both callers already gave there.
