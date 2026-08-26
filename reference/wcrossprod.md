# The Weighted Cross Product of the Assembly

Computes \\A^\top \mathrm{diag}(w) B\\, the weighted cross product the
information matrix and the outer Hessian are assembled from. Above a
measured work threshold, and only when the caller asked for more than
one thread, the work goes to the package's own threaded kernel.
Everywhere else it evaluates `crossprod(A * w, B)`, which is the
expression the kernel replaces.

## Usage

``` r
wcrossprod(A, w, B, threads = 1L)
```

## Arguments

- A, B:

  Dense design blocks with the same number of rows, `n x pa` and
  `n x pb`. Either may be sparse, in which case the Matrix route is
  taken and the result is whatever
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
the result is identical bit for bit at 1, 2 or 24 threads.

That guarantee is about the kernel alone. It says nothing about the two
routes agreeing with each other: engaging the kernel replaces a BLAS
call, and an optimized BLAS (OpenBLAS, Accelerate) blocks its
accumulations into a different order. The two routes therefore agree to
the last bit on R's reference BLAS and to the rounding of a single dot
product elsewhere.

## When the kernel is used

All four conditions must hold: `threads > 1`, `A` and `B` are base dense
matrices, `w` has one entry per row, and the work \\n \times p_A \times
p_B\\ is at least `2e5` multiply-adds. A sparse design keeps its Matrix
route, where a dense kernel would gain nothing: the cost there is set by
the number of stored nonzeros, and the threshold above reads the dense
shape.

The threshold is a constant in the source, with no argument to set it.
It was measured: the crossover where opening a parallel region costs
what it saves sits near \\9 \times 10^4\\ multiply-adds on the
development machine (0.94x at \\8 \times 10^4\\, 1.7x at \\10^5\\, 10x
at \\2 \times 10^6\\, 19x to 20x at the shapes a real fit assembles),
and the gate is set above it, where a misjudgment in either direction
costs almost nothing.

## Dimnames

The kernel returns a bare matrix, so the column names of `A` and `B` are
put back as the row and column names of the result when either has them.
[`crossprod()`](https://rdrr.io/r/base/crossprod.html) does this itself,
so both routes name their output alike.

## See also

[`xtv()`](https://statmodels7.github.io/statmodels7/reference/xtv.md)
and
[`xtx()`](https://statmodels7.github.io/statmodels7/reference/xtx.md),
the other two threaded assembly kernels, which share this threshold.
