# Which Coordinates a Boundary Has Frozen

The positions whose curvature is not finite, which is what a parameter
at the clamp its link keeps it strictly inside leaves behind.

## Usage

``` r
boundary_coords(K)
```

## Arguments

- K:

  A square matrix, the information or the penalized information.

## Value

An integer vector of positions, empty where there are none.

## Details

The DIAGONAL and not the column. A frozen coordinate makes its whole row
non-finite, cross terms included, so a column test marks its neighbours
too: measured on a Student t whose \\\nu\\ had reached `double.xmax`,
testing columns held \\\sigma\\ along with \\\nu\\ and left the fit
exactly where it had been.

## See also

[`pin_boundary`](https://statmodels7.github.io/statmodels7/reference/pin_boundary.md),
[`iwls_solve`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
