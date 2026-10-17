# The Values That a Path Visits Over a Hyperparameter Not Swept by Kink Size

Returns the grid of a hyperparameter that
[`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md)
does not reach: equally spaced points strictly inside its interval where
the interval is bounded, and a geometric grid above
[`shape_floor()`](https://statmodels7.github.io/statmodels7/reference/shape_floor.md)
where it has no upper bound.

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
kink is \\\lambda\alpha\\, so at a given \\\lambda\\ the block stays
nonempty for every admissible \\\alpha\\ and every point of such a path
is dropped. A bounded hyperparameter is therefore swept over its
interval, and \\\alpha\\ runs between the ridge and the lasso. The kink
still moves through \\\lambda\\, which is swept on its own axis.

The endpoints are excluded because the bounds are open: the elastic net
at \\\alpha = 0\\ has no kink at all, and the path would be scoring a
penalty of another kind.

A hyperparameter with no upper bound that does not move the kink is the
shape of SCAD or MCP. Its grid is the floor plus
`exp(seq(log(0.25), log(25), length.out = n_values))`, which covers the
values in common use (3.7 for the SCAD of Fan and Li, 2001, and 3 for
the MCP of Zhang, 2010) and does not reach the floor itself, where the
operator is set-valued or the penalty is not defined.

The floor is the smallest value at which the block can be fitted, which
[`shape_floor()`](https://statmodels7.github.io/statmodels7/reference/shape_floor.md)
derives from the proximal condition at the steps that the block's
coordinate descent takes, and not the constant over which the penalty is
defined. The two coincide on an ordinary well-conditioned block and
differ where the steps are long.

## See also

[`path_values()`](https://statmodels7.github.io/statmodels7/reference/path_values.md),
[`path_bounded()`](https://statmodels7.github.io/statmodels7/reference/path_bounded.md),
[`shape_floor()`](https://statmodels7.github.io/statmodels7/reference/shape_floor.md)
