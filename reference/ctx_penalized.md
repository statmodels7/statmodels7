# The Penalized Matrix and Its Inverse

Builds \\K = H + S\\ with the criterion's own information, and its
inverse, which the gradient and the Hessian both read.

## Usage

``` r
ctx_penalized(
  ctx,
  spec,
  design,
  coef,
  hyper,
  expected = FALSE,
  laplace = FALSE
)
```

## Arguments

- ctx:

  A context, or `NULL`.

- spec, design, coef, hyper:

  The fallback arguments.

- expected:

  Whether \\H\\ is the expected information, as the criterion carries
  under `reml(information = expected_info())`, or the observed one. The
  two give different matrices and keep separate cache entries.

- laplace:

  Whether the matrix is the criterion's, which leaves out the kinked
  coordinates, or the mode's.

## Value

A list with `K`, `inv` and `logdet`, or `NULL`.

## Details

The storage is the matrix's own. \\H\\ is sparse wherever the design is,
which covers a grouping indicator, a factor `by` and a `linpar` over
many levels. Where the sum \\H + S\\ is large enough and sparse enough
([`worth_sparse()`](https://statmodels7.github.io/statmodels7/reference/worth_sparse.md))
it is kept sparse and factorized as such. The decision depends only on
the matrix and not on which term produced it. The penalized matrix is a
new object at every point, so a factorization stored on an earlier
matrix is never reused.

The inverse stays dense whatever the factor is. Its readers take full
matrix products against it,
[`block_leverage()`](https://statmodels7.github.io/statmodels7/reference/block_leverage.md)
and the Hessian's pair loop among them, so a dense inverse of a sparse
matrix costs nothing that sparsity could have kept. The sparse factor
saves the cost of producing the inverse.

`NULL` is returned where the matrix is not positive definite. For the
mode's matrix (`laplace = FALSE`) the directions that
[`deficient_coords()`](https://statmodels7.github.io/statmodels7/reference/deficient_coords.md)
reports as unidentified are pinned and the factorization is attempted
once more before `NULL` is returned.

Where a penalty has a kink, there are two matrices. The criterion's
determinant leaves out every coordinate that a kinked penalty covers
([`laplace_pinned()`](https://statmodels7.github.io/statmodels7/reference/laplace_pinned.md)),
which is the matrix that `laplace = TRUE` selects. The default is the
mode's own curvature, which leaves out the coordinates that the marginal
criterion holds or does not identify
([`pinned_coords()`](https://statmodels7.github.io/statmodels7/reference/pinned_coords.md)),
the boundary coordinates
([`pin_boundary()`](https://statmodels7.github.io/statmodels7/reference/pin_boundary.md))
and the coordinates that a kinked penalty holds at zero
([`zero_kinked()`](https://statmodels7.github.io/statmodels7/reference/zero_kinked.md)):
those stay at zero when a hyperparameter moves, and a Newton correction
of the mode does not move them either. Where no penalty has a kink the
two are one matrix and one cache entry.

## See also

[`pd_factor()`](https://statmodels7.github.io/statmodels7/reference/pd_factor.md)
