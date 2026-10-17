# Which Coordinates a Boundary Has Frozen

Returns the positions whose curvature is not finite, as left by a
parameter sitting at the clamp of its link, which keeps it strictly
inside its bounds.

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

The test is on the diagonal. A frozen coordinate makes its whole row
non-finite, including the cross terms, so a test over whole columns
marks its neighbors too.

## See also

[`pin_boundary()`](https://statmodels7.github.io/statmodels7/reference/pin_boundary.md),
[`iwls_solve()`](https://statmodels7.github.io/statmodels7/reference/iwls_solve.md)
