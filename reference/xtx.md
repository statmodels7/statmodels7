# The Unweighted Cross Product of a Square-Root Design

Computes \\A^\top A\\ for a square-root design, through the package's
threaded kernel above the measured work threshold and through
`crossprod(A)` everywhere else. This is what
[`fit_smooth()`](https://statmodels7.github.io/statmodels7/reference/fit_smooth.md)'s
subset route assembles at every scoring iteration: the design has
already been scaled by the square root of the working weights, so no
weight vector enters.

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
walks \\j = 1, \ldots, k\\, and each entry it computes is written to
both \\(j,k)\\ and \\(k,j)\\. The two halves are therefore the same
number, one rounding shared, so the result is exactly symmetric with
nothing to symmetrize away. Each pair belongs to exactly one thread, so
the mirrored write is disjoint, and no accumulation is split, so the
answer is identical bit for bit at any thread count.

That halving is what the reference `dsyrk` behind
[`crossprod()`](https://rdrr.io/r/base/crossprod.html) does too, and it
accumulates in the same ascending row order, so the two routes agree to
the last bit there. An optimized BLAS blocks its accumulations and
agrees to the rounding of one dot product.

A full \\p^2\\ version was measured and is not what ships: it was no
better than `dsyrk` at eight threads, so the second half of the work
bought nothing.

## When the kernel is used

`threads > 1`, `A` a base dense matrix, and \\n \times p^2\\ at least
`2e5`. Note the \\p^2\\: this kernel writes \\p^2\\ outputs where
[`xtv()`](https://statmodels7.github.io/statmodels7/reference/xtv.md)
writes \\p\\, so a design wide enough to reach the threshold here can be
too narrow to reach it there.

## See also

[`wcrossprod()`](https://statmodels7.github.io/statmodels7/reference/wcrossprod.md)
for the weighted two-block form,
[`xtv()`](https://statmodels7.github.io/statmodels7/reference/xtv.md)
for the vector-valued kernels.
