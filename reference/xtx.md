# The Unweighted Cross Product of a Square-Root Design

Computes \\A^\top A\\ for a square-root design, through the package's
threaded kernel above a work threshold and through `crossprod(A)`
everywhere else. The argument is a design that has already been scaled
by the square root of the working weights, so the function takes no
weight vector.

## Usage

``` r
xtx(A, threads = 1L)
```

## Arguments

- A:

  A dense design block, `n x p`, already scaled by the square root of
  the working weights. A sparse `A` takes the Matrix route and returns
  whatever [`crossprod()`](https://rdrr.io/r/base/crossprod.html) gives
  for it.

- threads:

  The thread count, a plain integer as
  [`numericals7::thread_count()`](https://statmodels7.github.io/numericals7/reference/thread_count.html)
  returns it. `1L` takes the sequential route unconditionally.

## Value

A `p x p` symmetric numeric matrix, with the column names of `A` as both
dimnames when `A` has them.

## Symmetry, and independence of the thread count

Only the upper triangle is accumulated: one thread owns column \\k\\ and
walks \\j = 1, \ldots, k\\, and each entry that it computes is written
to both \\(j,k)\\ and \\(k,j)\\. The two halves are therefore the same
number, rounded once, so the result is exactly symmetric and needs no
symmetrization. Each pair belongs to exactly one thread, so the mirrored
write is disjoint. The kernel splits no accumulation, so the answer is
identical bit for bit at any thread count.

The reference `dsyrk` behind
[`crossprod()`](https://rdrr.io/r/base/crossprod.html) also halves the
work in this way and accumulates in the same ascending row order, so the
two routes agree to the last bit there. An optimized BLAS blocks its
accumulations and agrees to the rounding of one dot product.

Only the upper triangle is computed, as `dsyrk` does, because computing
all \\p^2\\ entries gave no gain.

## When the kernel is used

The kernel is used when `threads > 1`, `A` is a base dense matrix, and
\\n \times p^2\\ is at least `2e5`. This kernel writes \\p^2\\ outputs
where
[`xtv()`](https://statmodels7.github.io/statmodels7/reference/xtv.md)
writes \\p\\, so a design wide enough to reach the threshold here can be
too narrow to reach it there.

## See also

[`wcrossprod()`](https://statmodels7.github.io/statmodels7/reference/wcrossprod.md)
for the weighted two-block form,
[`xtv()`](https://statmodels7.github.io/statmodels7/reference/xtv.md)
for the vector-valued kernels.
