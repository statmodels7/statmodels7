# The Curvature a SCAD or MCP Penalty Is Scaled By

Writes into every SCAD and MCP penalty of the model the curvature of the
likelihood in each of its coordinates, read at the given coefficients.

## Usage

``` r
statmod_curv(spec, design, coef, expected, approx)
```

## Arguments

- spec, design:

  The specification and its design.

- coef:

  The coefficients, split by parameter.

- expected, approx:

  Which information.

## Value

The specification with the curvature written into the terms, or `NULL`
where no penalty was written, so that the caller rebuilds its blocks
only when a penalty changed.

## Details

The curvature of coordinate \\j\\ is the one the coordinate descent
steps with: \\c_j = \sum_i w_i (x\_{ij} - \bar x_j)^2\\, the working
weights being the ones
[`coord_working()`](https://statmodels7.github.io/statmodels7/reference/coord_working.md)
returns for the coordinate's own equation and \\\bar x_j\\ the column's
weighted mean where the equation's intercept is profiled out of the
descent
([`coord_centers()`](https://statmodels7.github.io/statmodels7/reference/coord_centers.md)),
and \\\sum_i w_i x\_{ij}^2\\ where it is not. Reading the uncentered
diagonal of the information instead, as 0.159.0 did, gave a curvature
the centered step does not have: on
[`MASS::UScrime`](https://rdrr.io/pkg/MASS/man/UScrime.html) with raw
predictors, whose means reach 33 standard deviations, the ratio of the
two reached 1138 and the step condition failed on every fit.

It is called at every pass of
[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
and once more at the fitted coefficients, which makes the scaling
SELF-CONSISTENT (decided 2026-09-29): at the fit, \\c_j\\ is the
curvature at the fit. Measured on 60 simulated sparse regressions it
gives 44 false positives against the 58 of a curvature read once at the
start, and it cannot fail the step condition,
[`coord_table_penalty()`](https://statmodels7.github.io/statmodels7/reference/coord_table_penalty.md)
writing the current step's curvature into the table. Under a diagonal
map \\u = D\beta\\ the curvature is divided by \\d_j^2\\, the penalty
being written on \\u\\. A penalty under any other map, or reached
through a sub-term, is left as it is.
