# The Two Vector Products of the Coordinate-Descent Loop

`xtv()` is `as.numeric(crossprod(X, v))`, one dot product per column;
`wxsq()` is `as.numeric(crossprod(w, X^2))`, the column curvatures of a
working model, computed without materializing the `n x p` square. Both
run through the package's threaded kernels above the same measured work
gate as
[`wcrossprod`](https://statmodels7.github.io/statmodels7/reference/wcrossprod.md)
and through the exact expressions they replace everywhere else.

## Usage

``` r
xtv(X, v, threads = 1L)

wxsq(X, w, threads = 1L)
```

## Arguments

- X:

  A design block, `n x p`.

- v, w:

  Per-observation vectors of length `n`.

- threads:

  The thread count, a plain integer.

## Value

A numeric vector of length `ncol(X)`.

## Details

Each output element is accumulated in full by one thread in ascending
row order, with any elementwise product rounded exactly as the replaced
expression rounds it (`v` arrives precomputed; the square is taken
before the weight multiplies it), so the result does not depend on the
thread count, bit for bit. These are the per-sweep reads the plan's
section 0quinquies classifies as reductions: they are parallel over the
OUTPUT elements, never over a split of one accumulation.
