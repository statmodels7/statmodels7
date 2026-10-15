# The Smallest Admissible Value of a Shape Parameter

The lower end of the range a shape may be swept over: the penalty's own
bound, raised to wherever its proximal operator starts to exist at the
steps the block's coordinate descent will take.

## Usage

``` r
shape_floor(pen, theta, name, steps = NULL)
```

## Arguments

- pen:

  A penalties7 penalty.

- theta:

  The hyperparameters in force.

- name:

  Which one is the shape.

- steps:

  What
  [`path_steps()`](https://statmodels7.github.io/statmodels7/reference/path_steps.md)
  returned, or `NULL`.

## Value

A single number: the penalty's lower bound where `steps` is `NULL` or
not finite, or where no admissible value is found below the upper bound;
otherwise the smallest value, to bisection accuracy, at which the
proximal table exists.

## Details

The limit is found by bisection on whether
[`penalties7::penalty_prox_spec()`](https://statmodels7.github.io/penalties7/reference/penalty_prox_spec.html)
returns a table for the penalty at the block's steps, so that a family
added later is covered and neither the \\a - 1\\ of SCAD nor the
\\\gamma\\ of MCP appears here. A grid starting just above the constant
over which the penalty is defined would otherwise contain points at
which the block cannot be fitted with those data by the only route that
a kinked penalty has.

## See also

[`path_steps()`](https://statmodels7.github.io/statmodels7/reference/path_steps.md),
[`path_grid()`](https://statmodels7.github.io/statmodels7/reference/path_grid.md)
