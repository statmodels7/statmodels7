# The Weighted Cross Product of the Assembly

Computes \\A^\top \mathrm{diag}(w) B\\, the weighted cross product the
information matrix and the outer Hessian are assembled from. Above a
work threshold, and only when more than one thread is requested, the
work goes to the package's own threaded kernel. Everywhere else the
function evaluates `crossprod(A * w, B)`, which is the expression the
kernel replaces.

## Usage

``` r
wcrossprod(A, w, B, threads = 1L)
```

## Arguments

- A, B:

  Design blocks with the same number of rows, `n x pa` and `n x pb`. If
  either is not a base dense matrix, the sequential route is taken and
  the result is whatever
  [`crossprod()`](https://rdrr.io/r/base/crossprod.html) gives for that
  pair.

- w:

  Per-observation weights, a numeric vector of length `n`. A weight of
  any sign is accepted; nothing here assumes the working weights are
  positive.

- threads:

  The thread count, a plain integer as
  [`numericals7::thread_count()`](https://statmodels7.github.io/numericals7/reference/thread_count.html)
  returns it. `1L`, the default, takes the sequential route
  unconditionally.

## Value

A `pa x pb` numeric matrix, with the column names of `A` and `B` as its
dimnames when either block carries them.

## Why the answer does not depend on the thread count

The kernel splits the work over the **elements of the output**: one
thread owns one entry \\(j,k)\\ and accumulates its whole dot product,
in ascending row order, with the weight multiplied onto `A`'s entry
before the product is formed. No accumulation is ever split between
threads, so the sum is the same sequence of additions at any count and
the result is identical bit for bit at any thread count.

That guarantee is about the kernel alone. It says nothing about the two
routes agreeing with each other: engaging the kernel replaces a BLAS
call, and an optimized BLAS (OpenBLAS, Accelerate) blocks its
accumulations into a different order. The two routes therefore agree to
the last bit on R's reference BLAS and to the rounding of a single dot
product elsewhere.

## When the kernel is used

All of the following must hold: `threads > 1`; `A` and `B` are base
dense matrices with the same number of rows; `w` has one entry per row;
and the work \\n \times p_A \times p_B\\ is at least `2e5`
multiply-adds. Otherwise `crossprod(A * w, B)` is evaluated, which
signals its own error for blocks of different heights. A sparse design
keeps its Matrix route, where a dense kernel would gain nothing: the
cost there is set by the number of stored nonzeros, and the threshold
reads the dense shape.

The threshold is a constant in the source, with no argument to set it.
It lies above the work at which opening a parallel region costs what it
saves.

## Dimnames

The kernel returns a bare matrix, so the column names of `A` and `B` are
put back as the row and column names of the result when either has them.
[`crossprod()`](https://rdrr.io/r/base/crossprod.html) does this itself,
so both routes name their output alike.

## See also

[`xtv()`](https://statmodels7.github.io/statmodels7/reference/xtv.md),
[`wxsq()`](https://statmodels7.github.io/statmodels7/reference/xtv.md)
and
[`xtx()`](https://statmodels7.github.io/statmodels7/reference/xtx.md),
the other threaded assembly kernels. They use the same constant, \\2
\times 10^5\\, on their own measures of work: \\n p\\ for
[`xtv()`](https://statmodels7.github.io/statmodels7/reference/xtv.md)
and
[`wxsq()`](https://statmodels7.github.io/statmodels7/reference/xtv.md),
and \\n p^2\\ for
[`xtx()`](https://statmodels7.github.io/statmodels7/reference/xtx.md).
