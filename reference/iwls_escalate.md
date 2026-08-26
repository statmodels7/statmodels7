# Raise the Levenberg Damping

The next \\\lambda\\ to try after a step that failed or moved nothing.

## Usage

``` r
iwls_escalate(damp, pieces)
```

## Arguments

- damp:

  The current damping, a single non-negative number. `0` starts the
  schedule.

- pieces:

  What
  [`iwls_pieces()`](https://statmodels7.github.io/statmodels7/reference/iwls_pieces.md)
  built, read only for the scale.

## Value

The next damping to try, a single positive number.

## Details

The first value is \\10^{-8}\\ of the curvature's own scale, as
[`iwls_scale()`](https://statmodels7.github.io/statmodels7/reference/iwls_scale.md)
measures it. That is negligible beside a well-curved coordinate and
already large beside one whose information has fallen by eight orders,
which is the disparity the damping exists to correct.

Each further escalation multiplies by a hundred, so eight of them span
sixteen orders and reach a damping that dominates the largest diagonal.
Measured on the case this was built for, the offending coordinate needed
\\\lambda \approx 100\\ against a largest diagonal of 2328, which the
sixth escalation passes.

## See also

[`iwls_scale()`](https://statmodels7.github.io/statmodels7/reference/iwls_scale.md)
for the scale,
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
for what the damping does to the system.
