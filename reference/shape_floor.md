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

A single number strictly above the penalty's lower bound.

## Details

The limit is derived from the condition, with nothing written down. The
question is put to the penalty, namely whether it produces a table at
this step, and bisected, so a family added later is covered and neither
the \\a - 1\\ of SCAD nor the \\\gamma\\ of MCP appears here. A grid
starting just above the constant the penalty is defined over would
otherwise contain points at which that block, with those data, cannot be
fitted by the only route a kinked penalty has.

## See also

[`path_steps()`](https://statmodels7.github.io/statmodels7/reference/path_steps.md),
[`path_grid()`](https://statmodels7.github.io/statmodels7/reference/path_grid.md)
