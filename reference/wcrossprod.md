# The Weighted Cross Product of the Assembly

`t(A) %*% (w * B)`, through the package's threaded kernel where the
thread count and the measured work threshold allow it, and through
`crossprod(A * w, B)` – exactly the expression it replaces – everywhere
else.

## Usage

``` r
wcrossprod(A, w, B, threads = 1L)
```

## Arguments

- A, B:

  Design blocks, `n x pa` and `n x pb`.

- w:

  The per-observation weights, length `n`.

- threads:

  The thread count, a plain integer.

## Value

A `pa x pb` matrix.

## Details

The kernel decomposes over the ELEMENTS of the output, each dot product
accumulated in full by one thread in ascending row order with the weight
multiplied onto `A`'s entry first, which is the same arithmetic in the
same order as the sequential product: the result does not depend on the
thread count, bit for bit. Only base dense matrices with a full-length
weight are eligible; a sparse design keeps its Matrix route, where a
threaded dense kernel would do nothing (piano_parallel.txt: the
threshold is on the work, p and the density, not on n alone).

The threshold is internal and measured, not an argument: the crossover
where the cost of opening the region meets its gain sits near `9e4`
multiply-adds on this machine's grid (0.94x at 8e4, 1.7x at 1e5, 10x at
2e6, 19-20x at the profile shapes), and the gate is set a factor of two
above it, where the asymmetry argument says an error either way costs
almost nothing.
