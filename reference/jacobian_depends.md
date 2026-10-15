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

The reading replaces the test `== 0`, an exact comparison on a quantity
that the chart produces by arithmetic. On a covariance that a class
carries, whose chart is `parameters7::dr_prod(2)`, a standard deviation
is a function of its own log coordinate only, so its entry in the
correlation column is structurally zero. The entry that is computed can
be a small nonzero number from cancellation, and an exact comparison
then reads it as a dependence, which removes the standard error and the
interval of the standard deviation for a dependence that does not exist.

A Jacobian entry is judged against the largest entry of its own row, a
row being a gradient, and the tolerance is `sqrt(.Machine$double.eps)`,
the point beyond which a differentiated quantity cannot be distinguished
from the rounding of the quantities it was formed from. The two possible
errors are asymmetric, so the tolerance errs on the generous side: a
genuine dependence of relative size below `sqrt(eps)` that is read as
zero changes a variance contribution by a relative amount of order
`1e-16`, whereas a cancellation artifact read as a dependence removes
the standard error of a quantity.

## See also

[`readable_hyper_rows()`](https://statmodels7.github.io/statmodels7/reference/readable_hyper_rows.md),
[`summary_blocks()`](https://statmodels7.github.io/statmodels7/reference/summary_blocks.md)
