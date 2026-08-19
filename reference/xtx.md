# The Unweighted Cross Product of a Square-Root Design

`crossprod(A)` through the threaded kernel where the count and the
measured work gate allow it, and through `crossprod(A)` itself
everywhere else. What
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)'s
subset route assembles at every scoring iteration.

## Usage

``` r
xtx(A, threads = 1L)
```

## Arguments

- A:

  A dense design block, `n x p`.

- threads:

  The thread count, a plain integer.

## Value

A `p x p` matrix.

## Details

Element \\(j, k)\\ is one dot product accumulated in full by one thread
in ascending row order, and elements \\(j, k)\\ and \\(k, j)\\ are the
same products summed in the same order, so the result is exactly
symmetric and the kernel does not depend on the thread count, bit for
bit. The reference `dsyrk` behind `crossprod` accumulates in the same
order; an optimized BLAS does not, and agrees to the rounding of one dot
product.
