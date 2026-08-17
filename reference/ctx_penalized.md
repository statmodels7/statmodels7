# The Penalized Matrix and Its Inverse

\\K = H + S\\ with the criterion's own information, and its inverse,
which the gradient and the Hessian both read.

## Usage

``` r
ctx_penalized(ctx, spec, design, coef, hyper, expected = FALSE)
```

## Arguments

- ctx:

  A context, or `NULL`.

- spec, design, coef, hyper:

  The fallback arguments.

- expected:

  Whether \\H\\ is the expected information, which is what the criterion
  carries under `reml(hessian = "expected")`. It used to be hard-coded
  to the observed one, correctly, because the exact gradient ran on no
  other route; admitting the expected route makes the criterion's
  determinant a different matrix, and reading the wrong one would be a
  gradient of the wrong function.

## Value

A list with `K`, `inv` and `logdet`, or `NULL`.

## Details

**The storage is the matrix's own.** \\H\\ is already sparse wherever
the design is – a grouping indicator, a factor `by`, a `linpar` over
many levels – and it was being densified here only because the penalty's
accumulator is a base matrix. Where the sum is large enough and sparse
enough to be worth it
([`worth_sparse`](https://statmodels7.github.io/statmodels7/reference/worth_sparse.md))
it is kept sparse and factorized as such: measured on a random intercept
over 500 levels, p = 503 at a density of 0.014, the factorization and
its log-determinant cost 0.102 ms against 10.811 ms dense, and the full
inverse 3.280 ms against 25.000 ms, each route timed WITH its own
factorization. End to end that is 1.25x at 500 levels and 2.01x at 1000,
the difference between the operation and the fit being the lesson this
package records three times over: removing the dearer half leaves the
cheaper one. Nothing here asks which term produced the matrix; both
quantities are read off the matrix.

⚠️ Those figures are the ones measured with Matrix's factorization CACHE
defeated.
[`Matrix::Cholesky`](https://rdrr.io/pkg/Matrix/man/Cholesky-methods.html)
stores its result in the matrix's `factors` slot, so a benchmark that
refactorizes the same object measures a cache hit – 0.004 ms rather than
0.102 – and reports the sparse route as three times better than it is. A
fit never gets that hit, the penalized matrix being a new one at every
point.

**The inverse stays dense whatever the factor is.** Its readers take
full matrix products against it –
[`block_leverage`](https://statmodels7.github.io/statmodels7/reference/block_leverage.md)
and the Hessian's pair loop – so the inverse of a sparse matrix being
dense costs nothing sparsity could have kept. What the sparse factor
buys is the cost of producing it.

`NULL` is returned where the matrix is not positive definite, which is
the answer both callers already gave there.

## See also

[`pd_factor`](https://statmodels7.github.io/statmodels7/reference/pd_factor.md)
