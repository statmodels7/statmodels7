# The Step That a Coordinate Descent Would Take on a Block

Returns \\1/v_j\\, with \\v_j = \sum_i w_i x\_{ij}^2\\, one per column,
at the coefficients in hand.

## Usage

``` r
path_steps(spec, design, block, beta, split)
```

## Arguments

- spec:

  A
  [`StatmodSpec()`](https://statmodels7.github.io/statmodels7/reference/StatmodSpec-class.md).

- design:

  The design.

- block:

  One entry of `statmod_blocks()$sparse`.

- beta:

  The current coefficients.

- split:

  The objective's own splitter, which puts a stacked coefficient vector
  back into one piece per distribution parameter.

## Value

A numeric vector, one step per column, or `NULL`.

## Details

This is the step that the sweeps of
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
use, and it determines whether a shape parameter is admissible: SCAD's
proximal operator requires the step to be below \\a - 1\\ and MCP's
below \\\gamma\\, where the step is multiplied by the squared entries of
the diagonal map that standardization applies. The lower limit of a
shape is therefore a property of the data, and not only the constant 2
or 1 over which the penalty is defined.

The working weights are those of
[`coord_working()`](https://statmodels7.github.io/statmodels7/reference/coord_working.md)
for the observed information at the coefficients in hand. The result is
`NULL` for a structural block, for a block without columns, and wherever
the working weights or the squared column sums are not finite and
positive; the caller reads `NULL` as the penalty's own bound and nothing
more. The function also returns `NULL` when any error occurs inside the
computation, since the grid built from the steps is a starting grid and
may be approximate.

## See also

[`shape_floor()`](https://statmodels7.github.io/statmodels7/reference/shape_floor.md),
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
