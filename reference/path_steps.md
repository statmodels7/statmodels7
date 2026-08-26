# The Step a Coordinate Descent Would Take on a Block

\\1/v_j\\ with \\v_j = \sum_i w_i x\_{ij}^2\\, one per column, at the
coefficients in hand.

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

This is the step the sweeps of
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
use, and it is what decides whether a shape parameter is admissible:
SCAD's proximal operator needs \\t \< a - 1\\ and MCP's \\t \< \gamma\\,
tightened to \\t d^2\\ under the diagonal map standardization applies.
So the useful lower limit of a shape is a property of the data, not the
constant 2 or 1 the penalty is defined above, and it binds at ordinary
settings: measured, a standardized penalty on a column of spread 20 at
\\n = 200\\ needs \\a \> 3\\, and a Poisson block whose fitted means are
near \\10^{-3}\\ needs \\a \> 11\\.

Everything is guarded, and `NULL`, which the caller reads as "the
penalty's own bound and nothing more", is the answer wherever the
working weights are not usable. A starting grid may be approximate; what
it may not do is fail.

## See also

[`shape_floor()`](https://statmodels7.github.io/statmodels7/reference/shape_floor.md),
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
