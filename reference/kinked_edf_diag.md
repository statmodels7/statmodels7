# The Count of Each Coordinate Where a Block Has a Kink

The diagonal of \\(H\_{AA} + S\_{AA})^{-1}H\_{AA}\\ over the coordinates
\\A\\ away from a kink, and zero elsewhere: the count
[`statmod_pe()`](https://statmodels7.github.io/statmodels7/reference/statmod_pe.md)
prices a point of a path with, one coordinate at a time.

## Usage

``` r
kinked_edf_diag(
  spec,
  coef,
  design,
  hyper,
  expected,
  approx,
  aliased = integer(0)
)
```

## Arguments

- spec, coef, design, hyper:

  The specification, the coefficients, the design and the
  hyperparameters.

- expected, approx:

  Which information.

- aliased:

  Stacked positions the fit did not estimate.

## Value

A numeric vector over the stacked coefficients, or `NULL`.

## Details

The total is
[`statmod_pe()`](https://statmodels7.github.io/statmodels7/reference/statmod_pe.md)'s
\\\tau\\ at the same coefficients: where the penalty's Hessian vanishes
on \\A\\ each active coordinate counts one, which is the number of
non-zero coefficients of a lasso, and otherwise the trace is read
through the same Cholesky factor. A coordinate in the concave zone of
SCAD or MCP has a negative penalty curvature and counts more than one,
which is Stein's count for those penalties.

It answers only for a model with a kinked block and no structural term,
and declines (`NULL`) where the active columns are not of full rank,
where
[`statmod_pe()`](https://statmodels7.github.io/statmodels7/reference/statmod_pe.md)
reads a rank in place of a count, or where the factor does not exist;
[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
then reads its other rules.

## See also

[`statmod_pe()`](https://statmodels7.github.io/statmodels7/reference/statmod_pe.md),
[`statmod_edf()`](https://statmodels7.github.io/statmodels7/reference/statmod_edf.md)
