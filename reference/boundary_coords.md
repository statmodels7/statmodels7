# Which Coordinates a Boundary Has Frozen

The positions whose curvature is not finite, as a parameter sitting at
the clamp its link keeps it strictly inside leaves behind.

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

The test is on the **diagonal**. A frozen coordinate makes its whole row
non-finite, cross terms included, so a test over whole columns marks its
neighbors too. Measured on a Student t whose \\\nu\\ had reached
`double.xmax`, a column test held \\\sigma\\ along with \\\nu\\ and left
the fit exactly where it had been.

## See also

[`pin_boundary()`](https://statmodels7.github.io/statmodels7/reference/pin_boundary.md),
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
