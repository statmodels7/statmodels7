# Which Coordinates a Path Point Has to Visit

The sequential strong rule: a coordinate whose gradient at the previous
point of the path is below \\2s_k - s\_{k-1}\\ is left out of the
sweeps.

## Usage

``` r
coord_screen(X, w, z, beta, s_now, s_prev, threads = 1L)
```

## Arguments

- X:

  The block's own columns, `n x p`, dense or `dgCMatrix`.

- w:

  The working weights, length `n`.

- z:

  The working response, length `n`.

- beta:

  The block's coefficients at the previous point of the path, length
  `p`.

- s_now:

  The size of the kink at this point, a single number.

- s_prev:

  The size of the kink at the previous point, or `NULL` when there is no
  previous point.

- threads:

  The thread count the gradient read may use, a plain integer.

## Value

An integer vector of one-based column indices to visit, in ascending
order. Every column when `s_prev` is `NULL` or either kink size is not
usable. A coordinate already away from zero is always kept, whatever its
gradient, so the rule can only ever add coordinates to the active set.
Never empty: where the test discards everything, the column with the
largest gradient is kept.

## Details

The rule rests on the gradient moving no faster than the threshold does.
That is an assumption and not a bound, so the rule screens without
proving, and the caller checks what it discarded. With no previous point
there is nothing to screen against and every coordinate is visited.

## See also

[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
