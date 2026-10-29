# The Two Vector Products of the Coordinate-Descent Loop

`xtv()` and `wxsq()` compute the two per-sweep reads of the compiled
coordinate descent; each returns one vector of length `p` over an
`n x p` design.

`xtv()` computes \\X^\top v\\, one dot product per column. It is the
gradient read: `v` is the running residual, or the residual times the
working weights.

`wxsq()` computes \\\sum_i w_i x\_{ij}^2\\ for each column \\j\\, the
column curvatures of a working model. The threaded kernel does not form
the `n x p` matrix of squares, whereas the sequential expression does.

Both take the package's threaded kernel above a work threshold, the one
constant that
[`wcrossprod()`](https://statmodels7.github.io/statmodels7/reference/wcrossprod.md)
also uses, and otherwise evaluate the expressions that they replace,
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

  The vector against which \\X^\top v\\ is taken, length `n`. Only
  `xtv()` reads it.

- threads:

  The thread count, a plain integer as
  [`numericals7::local_threads()`](https://statmodels7.github.io/numericals7/reference/local_threads.html)
  validates it. `1L` takes the sequential route unconditionally.

- w:

  Per-observation weights, length `n`. Only `wxsq()` reads it, and it
  multiplies the squared entry.

## Value

An unnamed numeric vector of length `ncol(X)`: for `xtv()` the column
dot products with `v`, for `wxsq()` the weighted sums of squares. Both
functions omit the column names of `X`, because their results feed
arithmetic indexed by position.

## Why the answer does not depend on the thread count

Each element of the output is one column's accumulation, owned in full
by one thread and summed in ascending row order. The kernel never splits
an accumulation between threads, so the result is identical bit for bit
at any count.

The elementwise product inside each term is formed in the order in which
the replaced expression forms it: `v` arrives already multiplied by
whatever the caller wanted, and in `wxsq()` the entry is squared before
the weight multiplies it. Against an optimized BLAS, which blocks its
own accumulations, the two routes agree to the rounding of one dot
product.

## When the kernel is used

The kernel is used when `threads > 1`, `X` is a base dense matrix, the
vector has length `nrow(X)`, and \\n \times p\\ is at least `2e5`. A
sparse `X` takes the Matrix route.

## See also

[`wcrossprod()`](https://statmodels7.github.io/statmodels7/reference/wcrossprod.md)
for the matrix-valued kernel and the work threshold,
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
for the descent that calls these.
