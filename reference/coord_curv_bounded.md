# Whether Every Column Curvature Is Finite and Positive, Decided Without Computing Them

`coord_curv_bounded()` tests whether every column curvature \\v_j =
\sum_i w_i x\_{ij}^2\\ of a block would be finite and positive, from the
block's column norms and the range of the weights, so that
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
need not compute the curvatures of columns that it never visits.
`coord_curv()` computes them on the columns given in `k`.

## Usage

``` r
coord_curv_bounded(X, w, design, p, cols, moves)

coord_curv(X, w, k, threads = 1L)
```

## Arguments

- X:

  The block, dense or `dgCMatrix`.

- w:

  The working weights, finite and positive.

- design:

  The design from which the block was read.

- p:

  The equation, a parameter name.

- cols:

  The block's column positions within the equation.

- moves:

  Whether any term recomputes its block with the coefficients.

- k:

  Integer column positions within `X`.

- threads:

  The thread count that
  [`wxsq()`](https://statmodels7.github.io/statmodels7/reference/xtv.md)
  would be given.

## Value

`coord_curv_bounded()` gives a single logical. `coord_curv()` gives a
numeric vector of `length(k)` curvatures.

## Details

The condition is sufficient and not necessary. The weights are finite
and positive, as
[`coord_working()`](https://statmodels7.github.io/statmodels7/reference/coord_working.md)
has already checked, so a column's curvature is at least \\\min_i w_i\\
\lVert x_j\rVert^2\\ and at most \\\max_i w_i\\ \lVert x_j\rVert^2\\.
Hence where every column norm is positive and \\\min w \min_j \lVert
x_j\rVert^2 \> 10^{-200}\\ and \\\max w \max_j \lVert x_j\rVert^2 \<
10^{200}\\, every curvature stays nonzero and does not underflow or
overflow. Where any of this fails, `FALSE` sends the caller to the full
computation and its own check, so the outcome is the one that the full
vector would have given.

The norms are cached in the design's `eta_memo` environment, keyed by
the equation and the block's columns. A design whose blocks move with
the coefficients carries no such environment, and there the answer is
always `FALSE`.

`coord_curv()` takes the route that
[`wxsq()`](https://statmodels7.github.io/statmodels7/reference/xtv.md)
would take for the whole block: each curvature is its own column's sum,
so it is the same number either way.

## See also

[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md),
[`wxsq()`](https://statmodels7.github.io/statmodels7/reference/xtv.md)
