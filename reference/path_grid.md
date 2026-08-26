# The Values a Path Visits Over a Bounded Hyperparameter

An equally spaced grid strictly inside the hyperparameter's own
interval.

## Usage

``` r
path_grid(pen, theta, name, n_values = 25L, steps = NULL)
```

## Arguments

- pen:

  A penalties7 penalty.

- theta:

  The hyperparameters in force.

- name:

  Which hyperparameter the path varies.

- n_values:

  How many points.

- steps:

  What
  [`path_steps()`](https://statmodels7.github.io/statmodels7/reference/path_steps.md)
  returned, or `NULL`.

## Value

A numeric vector of values for `name`.

## Details

[`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
walks the size of the kink, from the value that empties the block down,
and a bounded hyperparameter cannot reach that end: the elastic net's
kink is \\\lambda\alpha\\, so at a given \\\lambda\\ no admissible
\\\alpha\\ empties the block and every point of such a path is dropped.
The interval itself is what a bounded shape is swept over instead,
\\\alpha\\ between the ridge and the lasso and SCAD's and MCP's shape
between their limits, and the sweeps being cyclic, one coordinate at a
time, the kink still moves through \\\lambda\\.

The endpoints are excluded because the bounds are open: the elastic net
at \\\alpha = 0\\ has no kink at all, and the path would be scoring a
penalty of another kind.

A shape parameter is swept above the smallest value at which the block
can be fitted, which
[`shape_floor()`](https://statmodels7.github.io/statmodels7/reference/shape_floor.md)
derives from the proximal condition at the steps the block's coordinate
descent will take, rather than above the constant the penalty is defined
over. The two coincide on an ordinary well-conditioned block and differ
where the steps are long: with a standardized penalty on a column of
spread 20 the floor is 3 where SCAD's own bound is 2, so a quarter of
the old grid named shapes no fit could reach.

## See also

[`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md),
[`path_bounded()`](https://statmodels7.github.io/statmodels7/reference/path_bounded.md),
[`shape_floor()`](https://statmodels7.github.io/statmodels7/reference/shape_floor.md)
