# Which Readings a Chart's Coordinate Enters

Reads a readable block's Jacobian as a statement about dependence:
`TRUE` where a quantity depends on a coordinate, `FALSE` where the entry
is the floating-point image of a structural zero.

## Usage

``` r
jacobian_depends(J)
```

## Arguments

- J:

  The Jacobian of a readable block, quantities by coordinates.

## Value

A logical matrix of the same shape, `TRUE` where the quantity in that
row depends on the coordinate in that column. A row that is zero
throughout depends on nothing and comes back all `FALSE`.

## Details

Two callers ask this question and both used to ask it as `== 0`, which
is an exact test on a quantity the chart arrives at by arithmetic. ⚠️
**Measured, it is wrong.** On a covariance a class carries, whose chart
is `parameters7::dr_prod(2)`, a standard deviation is a function of its
own log coordinate and of nothing else, so its entry in the correlation
column is structurally zero – and at one fitted point that entry came
back `1.338e-23` against a row whose own size is `0.5056`, a relative
`2.6e-23`, while at another it came back exactly `0`. Read exactly, the
first blanks the standard deviation's standard error and its interval
for a dependence that is not there.

The scale a Jacobian entry means anything against is the row's own
largest entry, that row being a gradient, and the tolerance is
`sqrt(.Machine$double.eps)` – the point past which a differentiated
quantity cannot be told from the rounding of what it was formed from,
which is the constant this toolkit derives for that question everywhere
else. The two errors it arbitrates are wildly asymmetric, which is why
the generous side is the right one: a genuine dependence of relative
size below `sqrt(eps)` read as zero costs a variance contribution of
relative order `1e-16`, while a cancellation artifact read as a
dependence costs a quantity its standard error outright.

## See also

[`readable_hyper_rows()`](https://statmodels7.github.io/statmodels7/reference/readable_hyper_rows.md),
[`summary_blocks()`](https://statmodels7.github.io/statmodels7/reference/summary_blocks.md)
