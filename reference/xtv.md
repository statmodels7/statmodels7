# The Two Vector Products of the Coordinate-Descent Loop

The two per-sweep reads of the compiled coordinate descent, each one
vector of length `p` over an `n x p` design.

`xtv()` computes \\X^\top v\\, one dot product per column. It is the
gradient read: `v` is the running residual, or the residual times the
working weights.

`wxsq()` computes \\\sum_i w_i x\_{ij}^2\\ for each column \\j\\, the
column curvatures of a working model. The `n x p` matrix of squares is
never formed.

Both take the package's threaded kernel above the same measured work
threshold
[`wcrossprod()`](https://statmodels7.github.io/statmodels7/reference/wcrossprod.md)
uses, and otherwise evaluate the expressions they replace,
`as.numeric(crossprod(X, v))` and `as.numeric(crossprod(w, X^2))`.

## Usage

``` r
xtv(X, v, threads = 1L)

wxsq(X, w, threads = 1L)
```

## Arguments

- X:

  A design block, `n x p`.

- v:

  The vector \\X^\top v\\ is taken against, length `n`. Read by `xtv()`
  only.

- threads:

  The thread count, a plain integer as
  [`numericals7::thread_count()`](https://statmodels7.github.io/numericals7/reference/thread_count.html)
  returns it. `1L` takes the sequential route unconditionally.

- w:

  Per-observation weights, length `n`. Read by `wxsq()` only, and
  multiplied onto the squared entry.

## Value

An unnamed numeric vector of length `ncol(X)`: for `xtv()` the column
dot products with `v`, for `wxsq()` the weighted sums of squares.
Neither carries the column names of `X`, since both feed arithmetic
indexed by position.

## Why the answer does not depend on the thread count

Each element of the output is one column's accumulation, owned in full
by one thread and summed in ascending row order. Nothing is split
between threads, so the result is identical bit for bit at any count.

The elementwise product inside each term is formed in the order the
replaced expression forms it: `v` arrives already multiplied by whatever
the caller wanted, and in `wxsq()` the entry is squared before the
weight multiplies it. Against an optimized BLAS, which blocks its own
accumulations, the two routes agree to the rounding of one dot product.

## When the kernel is used

`threads > 1`, `X` a base dense matrix, the vector of length `nrow(X)`,
and \\n \times p\\ at least `2e5`. A sparse `X` takes the Matrix route.

## See also

[`wcrossprod()`](https://statmodels7.github.io/statmodels7/reference/wcrossprod.md)
for the matrix-valued kernel and the measured threshold both share,
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
for the descent that calls these.
