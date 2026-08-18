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

  The block's columns.

- w:

  The working weights.

- z:

  The working response.

- beta:

  The coefficients at the previous point.

- s_now:

  The size of the kink here.

- s_prev:

  The size of the kink at the previous point, or `NULL`.

- threads:

  The thread count the gradient read may use.

## Value

An integer vector of column indices, never empty.

## Details

The rule rests on the gradient moving no faster than the threshold,
which is an assumption rather than a bound, so it screens rather than
proves and the caller checks what it discarded. With no previous point
there is nothing to screen against and every coordinate is visited.

## See also

[`coord_fit`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
