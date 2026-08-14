# Is a Hyperparameter Swept by the Size of Its Kink?

TRUE where the geometric path of path_values() reaches it: the
hyperparameter scales the kink and has no upper bound, so the value that
empties the block is admissible.

## Usage

``` r
path_by_kink(pen, theta, name)
```

## Arguments

- pen:

  A penalties7 penalty.

- theta:

  The hyperparameters in force.

- name:

  Which hyperparameter.

## Value

A single logical.

## Details

The two conditions fail in different ways and both have to hold. The
elastic net's alpha scales the kink and is bounded by one, so no
admissible value of it empties the block at a given lambda; the shape of
SCAD and MCP has no upper bound and does not move the kink at all, so
the solve has nothing to solve. Either way the sweep is
[`path_grid`](https://statmodels7.github.io/statmodels7/reference/path_grid.md).

## See also

[`path_values`](https://statmodels7.github.io/statmodels7/reference/path_values.md),
[`path_grid`](https://statmodels7.github.io/statmodels7/reference/path_grid.md)
