# The Curvature by Which a SCAD or MCP Penalty Is Scaled

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

The specification in which the curvature is written into the terms, or
`NULL` where no penalty was written, so that the caller rebuilds its
blocks only when a penalty changed.

## Details

The curvature of coordinate \\j\\ is the one that the coordinate descent
uses in its steps. It is \\c_j = \sum_i w_i (x\_{ij} - \bar x_j)^2\\
where the equation's intercept is profiled out of the descent
([`coord_centers()`](https://statmodels7.github.io/statmodels7/reference/coord_centers.md)),
with \\w_i\\ the working weights that
[`coord_working()`](https://statmodels7.github.io/statmodels7/reference/coord_working.md)
returns for the coordinate's own equation and \\\bar x_j\\ the column's
weighted mean, and it is \\c_j = \sum_i w_i x\_{ij}^2\\ where the
intercept is not. The uncentered diagonal of the information is not
used, because it differs from the curvature of the centered step when
the column means are large relative to their standard deviations.

It is called at the starting coefficients, at every pass of
[`statmod_alternate()`](https://statmodels7.github.io/statmodels7/reference/statmod_alternate.md)
and once more at the fitted coefficients, so that the scaling is
self-consistent: at the fit, \\c_j\\ is the curvature at the fit. The
step condition of the descent is met because
[`coord_table_penalty()`](https://statmodels7.github.io/statmodels7/reference/coord_table_penalty.md)
writes the current step's curvature into the table. Under a diagonal map
\\u = D\beta\\ the curvature is divided by \\d_j^2\\, because the
penalty is written on \\u\\. A penalty under any other map, or reached
through a sub-term, is left as it is.
