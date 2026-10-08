# The Size of the Kink in Each Coordinate

Returns, for every coordinate of a kinked penalty, the jump of its
derivative at the kink, which is the threshold a coordinate's gradient
is compared with to decide whether it can stay at zero.

## Usage

``` r
coord_kinks(pen, theta, eps = 1e-04)
```

## Arguments

- pen:

  A kinked penalty.

- theta:

  Its hyperparameters.

- eps:

  The step of the one-sided derivatives.

## Value

A numeric vector with one entry per coordinate of `pen`, all zero where
the penalty reports no finite kink. Its first entry is
[`kink_scale()`](https://statmodels7.github.io/statmodels7/reference/kink_scale.md)'s
answer.

## Details

[`kink_scale()`](https://statmodels7.github.io/statmodels7/reference/kink_scale.md)
reads the first coordinate only, which is what a path needs to place its
grid. A coordinate descent needs every coordinate's: under a diagonal
map \\D\\, which is what `standardize` writes, coordinate \\j\\'s kink
is \\\lambda\lvert d_j\rvert\\. The jump is measured as in
[`kink_scale()`](https://statmodels7.github.io/statmodels7/reference/kink_scale.md),
by a Richardson extrapolation of the one-sided derivatives, on every
coordinate at once. A scaled SCAD or MCP has its kink at \\\lambda\lvert
d_j\rvert\\ whatever its curvature, the scaling leaving the slope at the
origin where it was.

## See also

[`kink_scale()`](https://statmodels7.github.io/statmodels7/reference/kink_scale.md),
[`coord_screen()`](https://statmodels7.github.io/statmodels7/reference/coord_screen.md),
[`coord_fit()`](https://statmodels7.github.io/statmodels7/reference/coord_fit.md)
